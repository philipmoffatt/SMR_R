library(yaml)
library(tidyverse)
library(sf)
library(terra)
library(FedData)
library(nhdplusTools)
library(whitebox)

## ---- load project config ----
cfg               <- yaml::read_yaml("config.yaml")
setwd(cfg$root_dir)
wshed_label       <- cfg$watershed$label
huc_filter_id     <- cfg$watershed$huc12_filter_id
basin_mask_value  <- cfg$watershed$basin_mask_value
wshed_id          <- cfg$watershed$wshed_id
res_vol           <- cfg$watershed$res_vol
res_coeff         <- cfg$watershed$res_coeff
depth_a           <- cfg$soil$depth_a_cm
residual_A        <- cfg$soil$residual_a
residual_B        <- cfg$soil$residual_b
stream_thresh     <- cfg$dem$stream_threshold
isobasin_thresh   <- cfg$dem$isobasin_threshold
smooth_filter     <- cfg$dem$smooth_filter
dem_res           <- cfg$dem$resolution_m
raw_path          <- cfg$data$raw_dir
processed_path    <- cfg$data$processed_dir
wshed_lat         <- cfg$watershed$lat
wshed_lon         <- cfg$watershed$lon
## ---------------------------

# ---- HUC12 boundary ----
huc_template_path <- file.path(raw_path, "template", "huc12.shp")
needs_huc <- !file.exists(huc_template_path) ||
             nrow(sf::st_read(huc_template_path, quiet = TRUE)) == 0
if (needs_huc) {
  message("Downloading HUC12 boundary for (", wshed_lat, ", ", wshed_lon, ") ...")
  point    <- sf::st_sfc(sf::st_point(c(wshed_lon, wshed_lat)), crs = 4326)
  huc12_dl <- nhdplusTools::get_huc(AOI = point, type = "huc12")
  if (nrow(huc12_dl) == 0) stop("nhdplusTools returned 0 features — check lat/lon in config.yaml.")
  huc12_dl <- huc12_dl[, "geometry"]
  huc12_dl$Id <- huc_filter_id
  dir.create(file.path(raw_path, "template"), recursive = TRUE, showWarnings = FALSE)
  sf::st_write(huc12_dl, huc_template_path, delete_layer = TRUE)
  message("Saved to ", huc_template_path)
}

PRW     <- sf::st_read(huc_template_path, quiet = TRUE)
MFC     <- PRW %>% dplyr::filter(Id == huc_filter_id)
MFC_buf <- PRW %>% dplyr::filter(Id == huc_filter_id) %>% sf::st_buffer(300)

# ---- Download data (FedData v4 API — no raw.dir / extraction.dir) ----
message("Downloading NED ...")
NED <- FedData::get_ned(template = MFC_buf, label = wshed_label, res = "1")

message("Downloading SSURGO ...")
ssurgo <- FedData::get_ssurgo(template = MFC_buf, label = wshed_label)

message("Downloading NLCD ...")
NLCD <- FedData::get_nlcd(template = MFC_buf, label = wshed_label,
                           year = 2019, dataset = "landcover")

# ---- Reproject NED to UTM Zone 10N at 30 m ----
NED_utm  <- terra::project(NED, "EPSG:32610", res = dem_res)

ned_dir  <- file.path(processed_path, "NED")
if (!dir.exists(ned_dir)) dir.create(ned_dir, recursive = TRUE)
dem_path <- file.path(ned_dir, paste0(wshed_label, "_NED_UTM.tif"))
terra::writeRaster(NED_utm, dem_path, overwrite = TRUE)
dem <- dem_path

# ---- WhiteBox DEM processing ----
whitebox::wbt_init()

outpath  <- file.path(processed_path, "dem")
dir_path <- file.path(outpath, wshed_label)
if (!dir.exists(outpath))  dir.create(outpath)
if (!dir.exists(dir_path)) dir.create(dir_path)

wbt_feature_preserving_smoothing(
  dem    = dem,
  output = file.path(dir_path, "dem_smoothed.tif"),
  filter = smooth_filter)

wbt_breach_depressions(
  dem    = file.path(dir_path, "dem_smoothed.tif"),
  output = file.path(dir_path, "dem_breached.tif"))

wbt_d8_pointer(
  dem    = file.path(dir_path, "dem_breached.tif"),
  output = file.path(dir_path, "dem_d8_pointer.tif"))

wbt_basins(
  d8_pntr = file.path(dir_path, "dem_d8_pointer.tif"),
  output  = file.path(dir_path, "dem_basins.tif"))

basins <- terra::rast(file.path(dir_path, "dem_basins.tif"))
plot(basins, main = "Basin IDs — confirm basin_mask_value in config.yaml")
message("Check basin plot. Current basin_mask_value = ", basin_mask_value,
        ". Update config.yaml if needed, then re-run from here.")

# Build watershed mask (1 = keep, NA = outside)
mfc_mask_01       <- basins == basin_mask_value
mfc_mask_01[is.na(mfc_mask_01)] <- 0   # replace NA with 0 before writing
terra::writeRaster(mfc_mask_01, file.path(dir_path, "dem_mfc_mask.tif"), overwrite = TRUE)

# terra mask helper: keep cells where mfc_mask_01 == 1
mask_to_wshed <- function(x) terra::mask(x, mfc_mask_01, maskvalues = 0)

wbt_d8_flow_accumulation(
  input  = file.path(dir_path, "dem_breached.tif"),
  output = file.path(dir_path, "d8_flow_accum.tif"))

wbt_extract_streams(
  flow_accum = file.path(dir_path, "d8_flow_accum.tif"),
  output     = file.path(dir_path, "dem_streams.tif"),
  threshold  = stream_thresh)

wbt_isobasins(
  dem    = file.path(dir_path, "dem_breached.tif"),
  output = file.path(dir_path, paste0("iso_basin_", isobasin_thresh, ".tif")),
  size   = isobasin_thresh)

# ---- Export DEM-derived layers as masked ASC ----
asc_path <- file.path(processed_path, "ASC")
mfc_path <- file.path(processed_path, "ASC", wshed_label)
if (!dir.exists(asc_path)) dir.create(asc_path)
if (!dir.exists(mfc_path)) dir.create(mfc_path)

dem_files <- list.files(dir_path, pattern = "\\.tif$", full.names = TRUE)
purrr::walk(dem_files, function(h) {
  out_file <- file.path(mfc_path, sub("\\.tif$", ".asc", basename(h)))
  temp     <- mask_to_wshed(terra::rast(h))
  terra::writeRaster(temp, filename = out_file, overwrite = TRUE)
})

# ---- Soils ----
majority_comps <- ssurgo$tabular$component %>%
  dplyr::group_by(mukey) %>%
  dplyr::slice(which.max(comppct.r))

hoz_col <- c("hzname","hzdept.r","hzdepb.r","hzthk.r","om.r",
             "dbovendry.r","dbthirdbar.r","ksat.r",
             "wsatiated.r","wthirdbar.r","wfifteenbar.r","cokey")

hoz_data <- ssurgo$tabular$chorizon[
  ssurgo$tabular$chorizon$cokey %in% majority_comps$cokey, ] %>%
  dplyr::select(dplyr::all_of(hoz_col)) %>%
  dplyr::mutate(hzthk.r = hzdepb.r - hzdept.r)

res_data <- ssurgo$tabular$corestrictions[
  ssurgo$tabular$corestrictions$cokey %in% majority_comps$cokey, ] %>%
  dplyr::select(c("reskind","resdept.l","reshard","cokey"))

hoz_comp <- c("om.r","dbovendry.r","dbthirdbar.r","ksat.r",
              "wsatiated.r","wthirdbar.r","wfifteenbar.r")

depth_weight <- function(data, cols, layer) {
  data <- data %>% dplyr::mutate(
    hzthk.r = hzdepb.r - hzdept.r,
    weights  = hzthk.r / sum(hzthk.r))
  (data %>% dplyr::select(dplyr::all_of(cols)) * data$weights) %>%
    dplyr::summarise_all(sum) %>%
    dplyr::mutate(depth = max(data$hzdepb.r), layer = layer, cokey = unique(data$cokey))
}

layer <- list()
for (i in unique(hoz_data$cokey)) {
  temp_hoz  <- hoz_data %>% dplyr::filter(cokey == i)
  temp_res  <- res_data %>% dplyr::filter(cokey == i)
  max_depth <- min(temp_res$resdept.l, max(temp_hoz$hzdepb.r))
  depth_a_i <- min(depth_a, max_depth)

  layer_A_temp <- temp_hoz %>% dplyr::filter(hzdept.r < depth_a_i) %>% dplyr::arrange(hzdept.r)
  layer_A_temp$hzdepb.r[nrow(layer_A_temp)] <- depth_a_i
  layer_A_out  <- depth_weight(layer_A_temp, hoz_comp, "A")

  if (max_depth > depth_a) {
    layer_B_temp <- temp_hoz %>%
      dplyr::filter(hzdepb.r > depth_a_i, hzdept.r <= max_depth) %>% dplyr::arrange(hzdept.r)
    layer_B_temp$hzdept.r[1] <- depth_a_i
    layer_B_temp$hzdepb.r[nrow(layer_B_temp)] <- max_depth
    layer_B_out <- depth_weight(layer_B_temp, hoz_comp, "B")
    layer <- rbind(layer, layer_A_out, layer_B_out)
  } else {
    layer <- rbind(layer, layer_A_out)
  }
}

MFC_spatial <- dplyr::left_join(
  ssurgo$spatial %>% dplyr::rename("mukey" = "MUKEY"),
  majority_comps %>% dplyr::select("compname","mukey","cokey") %>% dplyr::mutate(mukey = as.character(mukey))
) %>% dplyr::full_join(layer)

soil_cols <- names(MFC_spatial)[8:15]
for (i in unique(MFC_spatial$layer)) {
  for (j in soil_cols) {
    temp <- terra::rasterize(
      terra::project(terra::vect(MFC_spatial %>% dplyr::filter(layer == i)), "EPSG:32610"),
      NED_utm,
      field = j)
    temp_mfc <- mask_to_wshed(temp)
    terra::writeRaster(temp_mfc,
                       filename = file.path(mfc_path, paste0(j, "_", i, ".asc")),
                       overwrite = TRUE)
  }
}

# ---- Land cover ----
NLCD_reclass_ma <- matrix(c(
  11,1, 21,2, 22,2, 23,2, 24,2,
  31,2, 42,3, 52,4, 71,5, 81,5,
  82,6, 90,3, 95,5), ncol = 2, byrow = TRUE)

NLCD_reclass <- terra::classify(NLCD, NLCD_reclass_ma)
NLCD_utm     <- terra::project(NLCD_reclass, NED_utm, method = "near")
temp_mfc     <- mask_to_wshed(NLCD_utm)
terra::writeRaster(temp_mfc,
                   filename = file.path(mfc_path, "NLCD_simple.asc"),
                   overwrite = TRUE)

# ---- Unit conversions and renaming for Perl ----
imitate_smr_setup <- file.path(processed_path, "imitate_smr_setup")
copy_from_path    <- mfc_path
if (!dir.exists(imitate_smr_setup)) dir.create(imitate_smr_setup)

area_cells  <- sum(terra::values(mfc_mask_01) == 1, na.rm = TRUE)
wshed_frame <- data.frame(wshed_id, area_cells, res_vol, res_coeff)
write.table(wshed_frame, sep = " ", row.names = FALSE, col.names = FALSE,
            file = file.path(imitate_smr_setup, "wshed_res_properties.ini"))

copy_from <- c("dem_breached.asc","dem_streams.asc","dem_mfc_mask.asc","NLCD_simple.asc")
copy_to   <- c("el.asc","strms_30m.asc","watershed.asc","landuse.asc")
for (x in seq_along(copy_from)) {
  src <- file.path(copy_from_path, copy_from[x])
  dst <- file.path(imitate_smr_setup, copy_to[x])
  if (!file.exists(dst) && file.exists(src)) file.copy(src, dst)
}

# moisture content: % -> fraction
pct_files <- c("wfifteenbar.r_A.asc","wfifteenbar.r_B.asc",
               "wthirdbar.r_A.asc","wthirdbar.r_B.asc",
               "wsatiated.r_A.asc","wsatiated.r_B.asc")
mc_files  <- c("wiltpt_mc_A.asc","wiltpt_mc_B.asc",
               "fieldcap_mc_A.asc","fieldcap_mc_B.asc",
               "sat_mc_A.asc","sat_mc_B.asc")
for (x in seq_along(pct_files)) {
  dst <- file.path(imitate_smr_setup, mc_files[x])
  if (!file.exists(dst))
    terra::writeRaster(terra::rast(file.path(copy_from_path, pct_files[x])) / 100, dst)
}

# Ksat: um/s -> cm/day (*8.64)
for (s in c("A","B")) {
  dst <- file.path(imitate_smr_setup, paste0("Ksat_matrix_", s, ".asc"))
  if (!file.exists(dst))
    terra::writeRaster(terra::rast(file.path(copy_from_path, paste0("ksat.r_", s, ".asc"))) * 8.64, dst)
}

# Soil depths (stream cells capped at 1 cm for layer A)
soil_depth_fun <- function(streams_v, depth_v) ifelse(streams_v > 0, 1, depth_v)
for (s in c("A","B")) {
  dst <- file.path(imitate_smr_setup, paste0("soil_depth_", s, ".asc"))
  if (!file.exists(dst)) {
    depth_r  <- terra::rast(file.path(copy_from_path, paste0("depth_", s, ".asc")))
    streams  <- terra::rast(file.path(imitate_smr_setup, "strms_30m.asc"))
    depth_r[is.na(depth_r)]  <- 0
    streams[is.na(streams)]  <- 0
    soil_d   <- terra::lapp(c(streams, depth_r), fun = soil_depth_fun)
    soil_d[soil_d == 0] <- NA
    terra::writeRaster(soil_d, dst)
  }
}

r <- function(f) terra::rast(file.path(imitate_smr_setup, f))  # shorthand

combined_depth <- r("soil_depth_A.asc") + r("soil_depth_B.asc")
terra::writeRaster(combined_depth, file.path(imitate_smr_setup, "soil_depth.asc"), overwrite = TRUE)

sat_amt <- r("sat_mc_A.asc") * r("soil_depth_A.asc") +
           r("sat_mc_B.asc") * r("soil_depth_B.asc")
terra::writeRaster(sat_amt, file.path(imitate_smr_setup, "sat_amt.asc"), overwrite = TRUE)

wiltpt_amt <- r("wiltpt_mc_A.asc") * r("soil_depth_A.asc") +
              r("wiltpt_mc_B.asc") * r("soil_depth_B.asc")
terra::writeRaster(wiltpt_amt, file.path(imitate_smr_setup, "wiltpt_amt.asc"), overwrite = TRUE)

fieldcap_amt_A <- r("fieldcap_mc_A.asc") * r("soil_depth_A.asc") - r("soil_depth_A.asc") * residual_A
fieldcap_amt_B <- r("fieldcap_mc_B.asc") * r("soil_depth_B.asc") - r("soil_depth_B.asc") * residual_B
fieldcap_amt   <- fieldcap_amt_A + fieldcap_amt_B
terra::writeRaster(fieldcap_amt_A, file.path(imitate_smr_setup, "fieldcap_amt_A.asc"), overwrite = TRUE)
terra::writeRaster(fieldcap_amt_B, file.path(imitate_smr_setup, "fieldcap_amt_B.asc"), overwrite = TRUE)
terra::writeRaster(fieldcap_amt,   file.path(imitate_smr_setup, "fieldcap_amt.asc"),   overwrite = TRUE)

Ksubsurface <- r("Ksat_matrix_A.asc") / 500
terra::writeRaster(Ksubsurface, file.path(imitate_smr_setup, "Ksubsurface.asc"), overwrite = TRUE)

ETreduction_mc <- r("fieldcap_amt.asc") * 0.8 / r("soil_depth.asc")
terra::writeRaster(ETreduction_mc, file.path(imitate_smr_setup, "ETreduction_mc.asc"), overwrite = TRUE)

Ksat_mpore_A <- r("Ksat_matrix_A.asc") * 10
Ksat_mpore_B <- r("Ksat_matrix_B.asc") / 2
terra::writeRaster(Ksat_mpore_A, file.path(imitate_smr_setup, "Ksat_mpore_A.asc"), overwrite = TRUE)
terra::writeRaster(Ksat_mpore_B, file.path(imitate_smr_setup, "Ksat_mpore_B.asc"), overwrite = TRUE)

sat_mc_A <- r("sat_mc_A.asc"); sd_A <- r("soil_depth_A.asc")
sat_mc_B <- r("sat_mc_B.asc"); sd_B <- r("soil_depth_B.asc")
Kfc_A <- r("Ksat_matrix_A.asc") * exp((-13.0 / sat_mc_A) * (sat_mc_A - fieldcap_amt_A / sd_A))
Kfc_B <- r("Ksat_matrix_B.asc") * exp((-13.0 / sat_mc_B) * (sat_mc_B - fieldcap_amt_B / sd_B))
terra::writeRaster(Kfc_A, file.path(imitate_smr_setup, "Kfc_A.asc"), overwrite = TRUE)
terra::writeRaster(Kfc_B, file.path(imitate_smr_setup, "Kfc_B.asc"), overwrite = TRUE)

# ---- Flow direction percent maps (needed by Perl lateral flow routing) ----
el <- r("el.asc")

diff_n  <-  el - (terra::shift(el,  0,  30) %>% terra::extend(., el) %>% terra::crop(., el))
diff_ne <- (el - (terra::shift(el,  30, 30) %>% terra::extend(., el) %>% terra::crop(., el))) * 0.707
diff_e  <-  el - (terra::shift(el, 30,  0)  %>% terra::extend(., el) %>% terra::crop(., el))
diff_se <- (el - (terra::shift(el, 30,-30) %>% terra::extend(., el) %>% terra::crop(., el))) * 0.707
diff_s  <-  el - (terra::shift(el,  0,-30)  %>% terra::extend(., el) %>% terra::crop(., el))
diff_sw <- (el - (terra::shift(el,-30,-30) %>% terra::extend(., el) %>% terra::crop(., el))) * 0.707
diff_w  <-  el - (terra::shift(el,-30,  0)  %>% terra::extend(., el) %>% terra::crop(., el))
diff_nw <- (el - (terra::shift(el,-30, 30) %>% terra::extend(., el) %>% terra::crop(., el))) * 0.707

slope_delta <- c(diff_n, diff_ne, diff_e, diff_se, diff_s, diff_sw, diff_w, diff_nw)
names(slope_delta) <- c("n","ne","e","se","s","sw","w","nw")
slope_delta[slope_delta < 0] <- 0

flowunits <- sum(slope_delta)
pct_flow  <- slope_delta / flowunits
names(pct_flow) <- c("n","ne","e","se","s","sw","w","nw")

for (d in names(pct_flow)) {
  terra::writeRaster(pct_flow[[d]],
                     filename = file.path(imitate_smr_setup, paste0(d, ".asc")),
                     overwrite = TRUE)
}
terra::writeRaster(flowunits,
                   filename = file.path(imitate_smr_setup, "flowunits.asc"),
                   overwrite = TRUE)

message("Stage 1 complete. Outputs in: ", imitate_smr_setup)
message("Check basin_mask_value in config.yaml, then run wxData.R")
