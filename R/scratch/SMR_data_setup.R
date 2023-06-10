library(tidyverse)
library(mapview)
library(raster)
library(FedData)
library(sf)
library(terra)
library(whitebox)

## define MFC ----
raw_path = "./raw_data"

## boolean for whether to drop pullman
drop_pullman = TRUE

# grab dummy data with a covering extent: PRW raster. 
PRW <- st_read(file.path(raw_path, "template", "huc12.shp"))
  #   mapview(PRW)

# get MFC and buffer by 300 m or 10 cells to ensure we 
  #   capture the eatershed divide 
MFC <- PRW %>% filter(Id == "79")

MFC_buf <- PRW %>% filter(Id == "79") %>% st_buffer(., 300)
mapview(MFC)

mfc_no_pullman <- 
  st_read(
    file.path(
      raw_path, 
      "template", "mfc_no_pullman", "layers", "globalwatershed.shp")) %>%
  st_buffer(., 300) %>%
  st_transform(crs = st_crs(MFC_buf)) 
mapview(MFC_buf) + mapview(mfc_no_pullman)

mfc_no_pullman <- st_intersection(mfc_no_pullman, MFC_buf)

#if(drop_pullman == TRUE) {
#  MFC_buf <- st_intersection(MFC_buf, mfc_no_pullman)
#} else {
#  MFC_buf <- MFC_buf
#}

### pulling data 
NED <- FedData::get_ned(
  template = MFC_buf,
  res = 1,
  raw.dir = "./raw_data/NED",
  extraction.dir = "./processed_data/NED",
  label = "MFC_buf", 
  force.redo = F)

ssurgo <- FedData::get_ssurgo(
  template = MFC_buf,
  raw.dir = "./raw_data/SSURGO",
  extraction.dir = "./processed_data/SSURGO",
  label = "MFC", 
  force.redo = F)

NLCD <- FedData::get_nlcd(
  template = MFC_buf,
  year = 2019,
  dataset = "landcover",
  # raw.dir = "./raw_data/NLCD",
  # extraction.dir = "./procesed_data/NLCD",
  label = "MFC", 
  force.redo = F)


### white box DEM processing -----------------
## https://www.whiteboxgeo.com/manual/wbt_book/available_tools/hydrological_analysis.html#Basins
library(whitebox)
whitebox::wbt_init()
  #map temp 


NED_utm = projectRaster(from = NED,
                       crs = crs(MFC), 
                       res = 30, 
                       filename = gsub(NED$layer@file@name,
                                       pattern = "_1",
                                       replacement = "_UTM"),
                       overwrite = T)
dem <- NED_utm$layer@file@name

# build DEM directory 
outpath="./processed_data/dem"
if(!dir.exists(outpath)) dir.create(outpath)

# Build wshed specific out path directory 
dir_path <- file.path(outpath,"MFC")
  # make directory per site 
if(!dir.exists(dir_path)) dir.create(dir_path)
  

wbt_feature_preserving_smoothing(
    dem = dem,
    output = paste0(dir_path, "/dem_smoothed.tif"),
    filter = 9)

wbt_breach_depressions(dem = paste0(dir_path, "/dem_smoothed.tif"),
                         output = paste0(dir_path,"/dem_breached.tif"))

wbt_d8_pointer(dem = paste0(dir_path, "/dem_breached.tif"),
               output = paste0(dir_path,"/dem_d8_pointer.tif"))

# d8 <- raster(paste0(dir_path,"/dem_d8_pointer.tif"))
# plot(d8)

#preset basin size builds MFC - need to check pour point defaults 
wbt_basins(d8_pntr =paste0(dir_path,"/dem_d8_pointer.tif"),
           output = paste0(dir_path,"/dem_basins.tif")
          )

basins <- raster(paste0(dir_path,"/dem_basins.tif"))
plot(basins)

# build a boolean mask of MFC 
# this needs to be genrallized. 
mfc_mask = basins == 84
mfc_mask[is.na(mfc_mask[])] = 0
 #plot(mfc_mask)

#if(drop_pullman == TRUE) {
  # crop mfc_mask to extent of mfc_no_pullman
#  mfc_mask <- crop(mfc_mask, extent(mfc_no_pullman)) %>%
#    mask(., mfc_no_pullman_buffered)
#} else {
  # Leave mft_mask unchanged
#  mft_mask <- mft_mask
#}

mapview(mfc_mask)


raster::writeRaster(mfc_rast, paste0(dir_path,"/dem_mfc_mask.tif"), overwrite=T)


wbt_d8_flow_accumulation(input = paste0(dir_path,"/dem_breached.tif"),
                                  output = paste0(dir_path,"/d8_flow_accum.tif"))

wbt_extract_streams(flow_accum = paste0(dir_path,"/d8_flow_accum.tif"),
                    output = paste0(dir_path,"/dem_streams.tif"),
                    threshold = 1000) #

## explore the stream data 
# streams <- raster(paste0(dir_path,"/dem_streams.tif"))
# test <- (streams*mfc_rast)
# # mapview(mfc_rast)
# mapview(mfc_rast)+mapview(test)


# build small basins for portential parrellel computing. 
wbt_isobasins(dem = paste0(dir_path,"/dem_breached.tif"),
              output = paste0(dir_path,"/iso_basin_1000.tif"), 
              size = 1000)

### set to same threshold as streams and we can run all basins in parallel fast on kamiak maybe
iso = raster(paste0(dir_path,"/iso_basin_1000.tif"))
# test2 <- (mask(iso,MFC))
# 
# mapview(test2) + mapview(test)


# Prep terrain maps, NEED dem(pit filled), ASPECT, FLOWDISTIONS(check needed unit),
  # STREAMS, region, and watersheds.
# all mask to smallest extent and output as .asc
# my guess it will be fastest to run each of the small watershed totally independently. 
# need a standard directory stucture and we can can stack up the input data and split it up. 
 
dem_files <- list.files(dir_path, pattern = "tif",full.names = T)

asc_path = dir_path %>% gsub(., pattern = "dem/MFC", replacement = "ASC")
mfc_path = dir_path %>% gsub(., pattern = "dem/", replacement = "ASC/")

if(!dir.exists(asc_path)){dir.create(asc_path)}
if(!dir.exists(mfc_path)){dir.create(mfc_path)}


# Pass out all ASC
# we might need to clip these to the watershed before to reduce calcualtion time
# 

# mfc_mask
# test <- mask(NED_utm, mfc_mask, maskvalue = 1, incervse = T)

map(dem_files, function(h){
  print(h)
  temp <- raster(h) %>% mask(., mfc_mask, maskvalue = 1, inverse = T)
  writeRaster(temp, filename = gsub(
    gsub(h,pattern = "dem/", replacement = "ASC/"),
    pattern = ".tif", replacement = ".asc"), 
    overwrite=T)
})

## loop for sub watersheds -- do this later 

iso_mask = mask(iso, mfc_mask, maskvalue = 1, inverse = T)

# for(i in unique(iso_mask)){
#   temp_path = mfc_path %>% gsub(., pattern = "MFC", replacement = paste0("sub_",i))
#   if(!dir.exists(temp_path)){dir.create(temp_path)}
#   
#   map(dem_files[1:(length(dem_files)-2)], function(h){
#     # print(h)
#     temp <- raster(h) %>% mask(., mfc_mask, maskvalue = 1, inverse = T) %>% 
#       mask(., iso_mask==i, maskvalue = 1, inverse = T)
#     
#     #trim_NA, check shared extent and that it makes sense for post processing 
#     
#     print(
#       gsub(
#         gsub(h,pattern = "dem/MFC", replacement = paste0("ASC/sub_",i,"")),
#         pattern = ".tif", replacement = ".asc")
#       )
#       
#     temp <- temp %>% trim()
#     
#    #save out 
#    writeRaster(temp, filename = gsub(
#       gsub(h,pattern = "dem/MFC", replacement = paste0("ASC/sub_",i,"")),
#       pattern = ".tif", replacement = ".asc"), 
#       overwrite=T)
#   })
#   
# }




### soils-----------------------------------------

# how to treat ssurgo components majority or component weighted?
majority_comps = ssurgo$tabular$component %>% group_by(mukey) %>% 
  slice(which.max(comppct.r))

# check no dropped or duplicated mamp unit keys (mukeys)
if(!all((unique(ssurgo$tabular$component$mukey) %in% unique(majority_comps$mukey)))){
  warning(paste0("MUKEY [", 
                (unique(ssurgo$tabular$component$mukey)[
                  !(unique(ssurgo$tabular$component$mukey) %in% 
                     unique(majority_comps$mukey))]), 
                "] is not represented by a major component \n"))}

if(all(duplicated(majority_comps$mukey))){
  warning(paste0("MUKEY [", 
                 majority_comps$mukey[duplicated(majority_comps$mukey)], 
                 "] is duplicated\n"))
}

# ssurgo$tabular$chorizon$dbthirdbar.r %>% names

hoz_col <- c("hzname","hzdept.r", "hzdepb.r", "hzthk.r", "om.r",
             "dbovendry.r","dbthirdbar.r", "dbthirdbar.r","ksat.r",
             "wsatiated.r", "wthirdbar.r", "wfifteenbar.r", "cokey")
# 
# can pull soil texture percentages if needed for pedotransfer functions
# could pull H and L for calibration range on all soil parameters though often missing 

#pull majority component of each mukey 
hoz_data <- ssurgo$tabular$chorizon[
  ssurgo$tabular$chorizon$cokey %in% majority_comps$cokey,] %>% 
  # as.data.frame() %>% 
  dplyr::select(all_of(hoz_col))

# rebuild Hoz thickness
hoz_data <- hoz_data %>% mutate(hzthk.r = hzdepb.r - hzdept.r)

# pull restriction data for MFC 
res_data <- ssurgo$tabular$corestrictions[
  ssurgo$tabular$corestrictions$cokey %in% majority_comps$cokey,] %>% 
  dplyr::select(c("reskind","resdept.l","reshard","cokey"))
   # wbt_slope(dem = paste0(dir_path,"/dem_smoothed.tif"),
    #             output = paste0(dir_path,"/slope.tif"),
    #             units = "degrees")
    #   
# join hoz_data and res_data and treat depth (first guse weighted average 
# layer A 0-30cm or max/res if very shallow, 
# layer B =  30 - max or restriction), need to decide on what to do if B = 0 

# check if res depth is less than depth a 
# res_data$resdept.l %>% table
# # 25  41  46  50  51 100 102 
# # 2   1   2   2   6   2   4 

depth_a = 30 

#define parameters to depth weight
hoz_comp <- c("om.r","dbovendry.r","dbthirdbar.r","ksat.r",
"wsatiated.r","wthirdbar.r","wfifteenbar.r")

# direction weight weighting 
depth_weight <- function(data, cols, layer){
  data = data %>% mutate(
    hzthk.r = hzdepb.r - hzdept.r, 
    weights = hzthk.r/sum(hzthk.r, na.rm = T)
  )
  
  data_sum <- (data %>% 
                 dplyr::select(all_of(cols))*data$weights) %>% 
    summarise_all(function(H){mean(H, na.rm = TRUE)}) %>% 
    mutate(depth = max(data$hzdepb.r), 
           layer = layer, 
           cokey = data$cokey %>% unique())
}

#empty list for each build 
layer = list()

### build the list of soil characterisitcs in MF in layers A and B
for(i in unique(hoz_data$cokey)){
  layer_A_temp_2 <- list()
  layer_B_temp_2 <- list()
  
  temp_hoz = hoz_data %>% filter(cokey == i)
  temp_res = res_data %>% filter(cokey == i)
  
  # ! if there is no restriction in the component 
  # find the wighed average 
  max_depth = min(temp_res$resdept.l, max(temp_hoz$hzdepb.r))
    # we might make the max depth == 152cm which is the common SSURGO max depth.  
  depth_a = min(30, max_depth)  
          # select layers with topdepth < depth_a
  layer_A_temp = temp_hoz %>% filter(hzdept.r < depth_a) %>% 
    arrange(hzdept.r)
  layer_A_temp$hzdepb.r[nrow(layer_A_temp)] <- depth_a
  
  layer_A_temp[grepl(layer_A_temp$hzname, pattern = "O"),]$wsatiated.r = layer_A_temp[grepl(layer_A_temp$hzname, pattern = "O"),]$wthirdbar.r*1.1
  # need a better scalar for this...
  
  layer_A_temp_2 = depth_weight(data =layer_A_temp, cols = hoz_comp,
                                layer = "A")
  
  if(max_depth>30){
    layer_B_temp = temp_hoz %>% 
      filter(hzdepb.r > depth_a)%>% 
      filter(hzdept.r <= max_depth)%>% 
      arrange(hzdept.r)
    
    layer_B_temp$hzdept.r[1] <- depth_a
    layer_B_temp$hzdepb.r[nrow(layer_B_temp)] <- max_depth
    
    layer_B_temp= layer_B_temp %>% filter(!is.na(ksat.r))
    
    layer_B_temp_2 = depth_weight(layer_B_temp, cols = hoz_comp, layer = "B")
  }else if(max_depth <= 30){
    layer_B_temp_2 = layer_A_temp_2 %>% mutate(depth = 0.1, layer ="B")
    
  }
  rm(layer_A_temp)
  rm(layer_B_temp)
    
    layer = rbind(layer, layer_A_temp_2, layer_B_temp_2)
    rm(layer_A_temp_2, layer_B_temp_2)

      # build weights by thickness, weight bottom layer. (simple or harmonic weights?)
  }

  
## join by each layer (1 cokey = 1 layer + attributes)
# MFC_spatial_A = left_join(ssurgo$spatial %>% 
#                           rename("mukey" = "MUKEY"), 
#                         majority_comps %>% 
#                           select("compname", "mukey" , 'cokey') %>% 
#                           mutate(mukey = as.character(mukey))) %>% 
#   left_join(., layer %>% filter(layer=="A"))
# 
# MFC_spatial_B = left_join(ssurgo$spatial %>% 
#                             rename("mukey" = "MUKEY"), 
#                           majority_comps %>% 
#                             select("compname", "mukey" , 'cokey') %>% 
#                             mutate(mukey = as.character(mukey))) %>% 
#   left_join(., layer %>% filter(layer=="B"))


## join by each layer (1 cokey = 2 layers + attributes)
MFC_spatial = left_join(ssurgo$spatial %>% 
                            rename("mukey" = "MUKEY"), 
                          majority_comps %>% 
                            dplyr::select("compname", "mukey" , 'cokey') %>% 
                            mutate(mukey = as.character(mukey))) %>% 
  full_join(., layer)


# build out rasters for each MFC and ISO watershed, 
# **Check extent between hydro and soil asc

for(i in unique(MFC_spatial$layer)){
  for(j in names(MFC_spatial)[8:15]){
    temp = terra::rasterize(MFC_spatial %>% 
                              filter(layer==i) %>% 
                              vect() %>% 
                              project(., crs(NED_utm)), 
                            NED_utm %>% rast,
                            field = j)
    
    
    #build MFC layers 
    temp_mfc <-raster(temp) %>% mask(., mfc_mask, maskvalue = 1, inverse = T) 
   
     writeRaster(temp_mfc, filename = paste0(
       "./processed_data/ASC/MFC/",j,"_",i,".asc"), 
      overwrite=T)
  }
}

### land cover ----- 
# reclassify to row_crop, forest, grass, other
# unique(NLCD)
# [1] 11 21 22 23 24 31 42 52 71 81 82 90 95
NLCD_reclass_ma =matrix(data = c(
  11, 1, #water
  21, 2, #urban/other
  22, 2, #urban/other
  23, 2, #urban/other
  24, 2, #urban/other
  31, 2, #rock/barren/other
  42, 3, #forest
  52, 4, #shrub
  71, 5, #grass
  81, 5, #grass/pasture
  82, 6, # row crop
  90, 3, #woody wetlands
  95, 5), # grassy wetland
  ncol =2, byrow = T)

NLCD_reclass <- reclassify(NLCD, NLCD_reclass_ma) %>% projectRaster(., NED_utm, method = "ngb")
temp_mfc <- NLCD_reclass %>% mask(., mfc_mask, maskvalue = 1, inverse = T) 

writeRaster(temp_mfc, filename = paste0(
  "./processed_data/ASC/MFC/NLCD_simple.asc"), 
  overwrite=T)

## clip all thefiles at once 
## loop for sub watersheds -- do this later 
#  load all MFC asc 

mfc_asc_files = list.files(path = "./processed_data/ASC/MFC/", full.names = T)
tst <- rast(mfc_asc_files)
iso_mask = mask(iso, mfc_mask, maskvalue = 1, inverse = T)

#for(k in unique(iso_mask)){
#    temp_path = mfc_path %>% gsub(., pattern = "MFC", replacement = paste0("sub_",k))
#    if(!dir.exists(temp_path)){dir.create(temp_path)}

#      temp <- tst %>% crop(., rast(iso_mask)) %>% 
#        mask(., rast(iso_mask)==k, maskvalue = 1, inverse = T) %>% trim
#      print(temp_path)
#      writeRaster(temp_mfc, filename = paste0(temp_path,"/",j,"_",i,".asc"),
#        overwrite=T)
  
#  }

# re imagine the building of iso build full set then load as a stack and mask all at once
#

print("The asc file package complete")


### ----- Unit Conversions and Renaming for Perl Script ----- ###

## all maps need to be cropped to mfc_no_pullman
mapview(mfc_no_pullman)

# setting up out path to folder for unit converted maps and watershed properties file
imitate_smr_setup = file.path("./processed_data", "imitate_smr_setup")
copy_from_path = file.path ("./processed_data/ASC/MFC")

if(!dir.exists(imitate_smr_setup)) { 
  dir.create(imitate_smr_setup)
}

# builds out the tif to set mapset CRS using the watershed mask
file.copy(from = "./processed_data/dem/MFC/dem_mfc_mask.tif",
          to = "./processed_data/imitate_smr_setup/watershed.tif")

watershed_unmasked <- raster("./processed_data/dem/MFC/dem_mfc_mask.tif")

if(drop_pullman == TRUE) {
  watershed <- crop(watershed_unmasked, extent(mfc_no_pullman)) %>%
    mask(., mfc_no_pullman)
  watershed <- crop(watershed_unmasked, extent(MFC)) %>%
    mask(., MFC)
} else {
  watershed <- watershed_unmasked
}

watershed[watershed==0] <- NA
mapview(watershed)

raster::writeRaster(watershed, "./processed_data/imitate_smr_setup/watershed.tif", overwrite=T)

# set values for watershed properties
wshed_id = 79
area_cells = freq(watershed, value = 1)
res_vol = 110
res_coeff = 0.05

# build frame from values
wshed_frame = data.frame(wshed_id, area_cells, res_vol, res_coeff)

# write out .ini file to imitation folder
write.table(wshed_frame, sep = " ", row.names = F, col.names = F, file = paste0(imitate_smr_setup, "/wshed_res_properties.ini"))

# rename and write out the temp_mfc raster as 'landuse'
if(drop_pullman == TRUE) {
  temp_mfc <- crop(temp_mfc, extent(mfc_no_pullman)) %>%
    mask(., mfc_no_pullman)
} else {
  temp_mfc <- temp_mfc
}

writeRaster(temp_mfc, './processed_data/imitate_smr_setup/landuse.asc', datatype = "INT4S", overwrite = T, )

# copying over and renaming rasters so they fit into SMR perl
copy_from = c("/dem_breached.tif", "/dem_streams.tif")
copy_to = c("/el.asc", "/strms_30m.asc")

dem_path <- "./processed_data/dem/MFC"

# rasters that just need a name change
for (x in 1:length(copy_from)) {
    path_to_file <- paste0(dem_path, copy_from[x])
   
    #print(paste0("raster copied from: ", path_to_file))
    #print(copy_from[x])
    converted <- raster(path_to_file)
    
    if(drop_pullman == TRUE) {
      converted <- crop(converted, mfc_no_pullman) %>%
        mask(., mfc_no_pullman)
      converted <- crop(converted, MFC) %>%
        mask(., MFC)
    } else {
      converted <- converted
    }
    
    plot(converted, main=copy_from[x])
    
    copy_to_path <- paste0(imitate_smr_setup, copy_to[x])
    print(copy_to_path)
    
    #file.copy(paste0(copy_from_path, copy_from[x]), paste0(imitate_smr_setup, copy_to[x]))
    print(paste0('writing raster: ', copy_to[x]))
    writeRaster(converted, filename = copy_to_path, overwrite = T)
  #}
}

# stream fix
streams <- raster(path_to_file) 
streams[is.na(streams)] <- 0

crs(streams) <- crs(mfc_no_pullman)

streams <- crop(streams, extent(mfc_no_pullman)) %>% 
  mask(., mfc_no_pullman)
streams <- crop(streams, extent(MFC)) %>% 
  mask(., MFC)

writeRaster(streams, copy_to_path, overwrite=T)

# converting from % to moisture content (/100)
percentages = c("/wfifteenbar.r_A.asc", "/wfifteenbar.r_B.asc",
                "/wthirdbar.r_A.asc", "/wthirdbar.r_B.asc",
                "/wsatiated.r_A.asc", "/wsatiated.r_B.asc"
                )

moisture_contents = c("/wiltpt_mc_A.asc", "/wiltpt_mc_B.asc",
                        "/fieldcap_mc_A.asc", "/fieldcap_mc_B.asc",
                        "/sat_mc_A.asc", "/sat_mc_B.asc"
                        )

for (x in 1:length(percentages)) {
  #if (!file.exists(paste0(imitate_smr_setup, moisture_contents[x]))) {
    converted <- raster(paste0(copy_from_path, percentages[x])) / 100

    if(drop_pullman == TRUE) {
      converted <- crop(converted, extent(mfc_no_pullman)) %>%
        mask(., mfc_no_pullman)
    } else {
      converted <- converted
    }
    plot(converted)
    
    print(paste0('writing raster: ', moisture_contents[x]))
    raster::writeRaster(converted, paste0(imitate_smr_setup, moisture_contents[x]), overwrite = T)
  #}
}

# converting from um/second to cm/day (*8.64)
um_per_second = c("/ksat.r_A.asc", "/ksat.r_B.asc")
cm_per_day = c("/Ksat_matrix_A.asc", "/Ksat_matrix_B.asc")

for (x in 1:length(um_per_second)) {
  #if (!file.exists(paste0(imitate_smr_setup, cm_per_day[x]))) {
    converted <- raster(paste0(copy_from_path, um_per_second[x])) * 8.64
    
    if(drop_pullman == TRUE) {
      converted <- crop(converted, extent(mfc_no_pullman)) %>%
        mask(., mfc_no_pullman)
    } else {
      converted <- converted
    }
    plot(converted)
    
    print(paste0('writing raster: ', cm_per_day[x]))
    
    raster::writeRaster(converted, paste0(imitate_smr_setup, cm_per_day[x]), overwrite = T)
 # }
}

# calculate soil depths based based on stream cells
initial_depth_names = c("/depth_A.asc", "/depth_B.asc")
output_depth_names = c("/soil_depth_A.asc", "/soil_depth_B.asc")

soil_depth_fun <- function(x, y) {
  z <- ifelse(x>0, 1, y) 
  return(z)
  }

for (x in 1:length(initial_depth_names)) {
  #if (!file.exists(paste0(imitate_smr_setup, output_depth_names[x]))) {
    
    horizon <- ifelse(grepl("A", initial_depth_names[x]), "A", ifelse(grepl("B", initial_depth_names[x]), "B", NA))
    print(horizon)
    landuse <- raster(file.path(imitate_smr_setup, "landuse.asc"))
    plot(landuse)
    
    depth_rast <- raster(paste0(copy_from_path, initial_depth_names[x]))
    plot(depth_rast)
    streams <- raster(paste0(imitate_smr_setup, "/strms_30m.asc"))
    
    #depth_rast[is.na(depth_rast)] <- 0  # Set NA values to 0
    #streams[!is.na(mfc_mask)] <- 0  # Set NA values to 0
    
    #soil_depth <- overlay(streams, depth_rast, fun=soil_depth_fun)
    #soil_depth[soil_depth == 0] <- NA

    if(drop_pullman == TRUE) {
      depth_rast <- crop(depth_rast, extent(mfc_no_pullman)) %>%
        mask(., mfc_no_pullman)
      streams <- crop(streams, extent(mfc_no_pullman)) %>%
        mask(., mfc_no_pullman)
    }
    
    plot(streams)
    
    if (horizon=="A") {
      depth_rast[streams>0] <- 1
      depth_rast[landuse==2] <- 2
      plot(depth_rast)
    } else if (horizon=="B") {
      depth_rast[streams>0] <- 19
      depth_rast[landuse==2] <- 20
      plot(depth_rast)
    } else {
      depth_rast = depth_rast
    }
    

    print(paste0("writing out streams raster with drop_pullman being: ", drop_pullman))
    raster::writeRaster(streams, file.path(imitate_smr_setup, "strms_30m.asc"), overwrite=T)
    
    print(paste0("writing out renamed raster: ", output_depth_names[x]))
    raster::writeRaster(depth_rast, paste0(imitate_smr_setup, output_depth_names[x]), overwrite=T)
  #}
}

# fixing stream values
streams <- raster(paste0(imitate_smr_setup, "/strms_30m.asc"))
#streams[is.na(streams)] <- 0

#mfc_no_pullman <- crop(mfc_mask, mfc_no_pullman) --> don't think i need/want this
#mfc_no_pullman[mfc_no_pullman == 0] <- NA

if(drop_pullman == TRUE) {
  streams <- crop(streams, extent(mfc_no_pullman)) %>%
    mask(., mfc_no_pullman)
  #streams[is.na(mfc_no_pullman)] <- NA
  converted <- streams
} else {
  streams <- streams
}

mapview(streams)
writeRaster(streams, file.path(imitate_smr_setup, "strms_30m.asc"), overwrite=T)

# combined soil depth (A+B)
combined_depth <- raster(paste0(imitate_smr_setup, "/soil_depth_A.asc")) + 
                  raster(paste0(imitate_smr_setup, "/soil_depth_B.asc"))

raster::writeRaster(
  combined_depth,
  paste0(imitate_smr_setup, "/soil_depth.asc"),
  overwrite = T
  )

# calculate the combined saturation amount
sat_amt <- raster(paste0(imitate_smr_setup, "/sat_mc_A.asc")) * 
           raster(paste0(imitate_smr_setup, "/soil_depth_A.asc")) + 
           raster(paste0(imitate_smr_setup, "/sat_mc_B.asc")) * 
           raster(paste0(imitate_smr_setup, "/soil_depth_B.asc"))

raster::writeRaster(
  sat_amt,
  paste0(imitate_smr_setup, "/sat_amt.asc"),
  overwrite = T
  )

# calculate the combined wilting point amount
wiltpt_amt <- 
  raster(paste0(imitate_smr_setup, "/wiltpt_mc_A.asc")) * 
  raster(paste0(imitate_smr_setup, "/soil_depth_A.asc")) + 
  raster(paste0(imitate_smr_setup, "/wiltpt_mc_B.asc")) * 
  raster(paste0(imitate_smr_setup, "/soil_depth_B.asc"))

raster::writeRaster(
  wiltpt_amt,
  paste0(imitate_smr_setup, "/wiltpt_amt.asc"),
  overwrite = T
)

# setting residual floats for calculating fieldcap_amts
residual_A <- 0.02
residual_B <- 0.02

fieldcap_amt_A <- 
  raster(paste0(imitate_smr_setup, "/fieldcap_mc_A.asc")) * 
  raster(paste0(imitate_smr_setup, "/soil_depth_A.asc")) -
  raster(paste0(imitate_smr_setup, "/soil_depth_A.asc")) *
  residual_A

fieldcap_amt_B <- 
  raster(paste0(imitate_smr_setup, "/fieldcap_mc_B.asc")) * 
  raster(paste0(imitate_smr_setup, "/soil_depth_B.asc")) -
  raster(paste0(imitate_smr_setup, "/soil_depth_B.asc")) *
  residual_B

fieldcap_amt <- fieldcap_amt_A + fieldcap_amt_B

raster::writeRaster(fieldcap_amt_A, paste0(imitate_smr_setup, "/fieldcap_amt_A.asc"), overwrite = T)
raster::writeRaster(fieldcap_amt_B, paste0(imitate_smr_setup, "/fieldcap_amt_B.asc"), overwrite = T)
raster::writeRaster(fieldcap_amt, paste0(imitate_smr_setup, "/fieldcap_amt.asc"), overwrite = T)

# make Ksubsurface from Ksat_matrix_A
Ksubsurface <- 
  raster(paste0(imitate_smr_setup, "/Ksat_matrix_A.asc")) / 500

raster::writeRaster(Ksubsurface, paste0(imitate_smr_setup, "/Ksubsurface.asc"), overwrite = T)

# make ETreduction from fieldcap_amt and total soil_depth
ETreduction_mc <-
  raster(paste0(imitate_smr_setup, "/fieldcap_amt.asc")) *
  0.8 /
  raster(paste0(imitate_smr_setup, "/soil_depth.asc"))

raster::writeRaster(ETreduction_mc, paste0(imitate_smr_setup, "/ETreduction_mc.asc"), overwrite = T)

# make Ksat_mpores from Ksat_matrices
Ksat_mpore_A <- 
  raster(paste0(imitate_smr_setup, "/Ksat_matrix_A.asc")) *
  10

Ksat_mpore_B <- 
  raster(paste0(imitate_smr_setup, "/Ksat_matrix_B.asc")) /
  2

raster::writeRaster(Ksat_mpore_A, paste0(imitate_smr_setup, "/Ksat_mpore_A.asc"), overwrite = T)
raster::writeRaster(Ksat_mpore_B, paste0(imitate_smr_setup, "/Ksat_mpore_B.asc"), overwrite = T)

# make Kfc_A and Kfc_B with: Ksat_matrix_A/B, sat_mc_A/B, fieldcap_amt_A/B,
# and soil_depth_A/B --> using Bresler's formula

Kfc_A <-
  raster(paste0(imitate_smr_setup, "/Ksat_matrix_A.asc")) *
  exp(
    (-13.0 / raster(paste0(imitate_smr_setup, "/sat_mc_A.asc"))) * 
    (raster(paste0(imitate_smr_setup, "/sat_mc_A.asc")) - 
       fieldcap_amt_A/raster(paste0(imitate_smr_setup, "/soil_depth_A.asc")))
    )

Kfc_B <-
  raster(paste0(imitate_smr_setup, "/Ksat_matrix_B.asc")) *
  exp(
    (-13.0 / raster(paste0(imitate_smr_setup, "/sat_mc_B.asc"))) * 
      (raster(paste0(imitate_smr_setup, "/sat_mc_B.asc")) - 
         fieldcap_amt_B/raster(paste0(imitate_smr_setup, "/soil_depth_B.asc")))
  )

raster::writeRaster(Kfc_A, paste0(imitate_smr_setup, "/Kfc_A.asc"), overwrite = T)
raster::writeRaster(Kfc_B, paste0(imitate_smr_setup, "/Kfc_B.asc"), overwrite = T)


# FLOW DIRECTION PROGRAM
#    A "flow unit" is an elevation difference of one unit for an adjacent
#    cell.  An elevation difference of one unit for a diagonal cell would
#    be 1/(square root of 2)=0.707  flow units.

el <- rast("./processed_data/imitate_smr_setup/el.asc")

if(drop_pullman == TRUE) {
  el <- raster(el)
  el <- crop(el, extent(mfc_no_pullman)) %>%
    mask(., mfc_no_pullman)
  el <- rast(el)
}

#el <- rast(el)

# Calculate the differences with the 8-neighbors cells 

## check direction is not flipped ####
diff_n <- el - (terra::shift(el,0, 30) %>% terra::extend(.,el) %>% terra::crop(.,el))
diff_ne <- (el - (terra::shift(el, 30, 30) %>% terra::extend(.,el) %>% terra::crop(.,el)))* 0.707
diff_e <- el - (terra::shift(el, 30, 0) %>% terra::extend(.,el) %>% terra::crop(.,el))
diff_se <- (el - (terra::shift(el, 30, -30) %>% terra::extend(.,el) %>% terra::crop(.,el)))* 0.707
diff_s <- el - (terra::shift(el, 0,-30) %>% terra::extend(.,el) %>% terra::crop(.,el))
diff_sw <- (el - ((terra::shift(el, -30, -30) %>% terra::extend(.,el) %>% terra::crop(.,el))))* 0.707
diff_w <- el - (terra::shift(el, -30, 0) %>% terra::extend(.,el) %>% terra::crop(.,el))
diff_nw <- (el - (terra::shift(el, -30, 30) %>% terra::extend(.,el) %>% terra::crop(.,el)))* 0.707


# Take the maximum with 0.0
# 
# diff_ne[diff_ne < 0] <- 0
# diff_e[diff_e < 0] <- 0
# diff_se[diff_se < 0] <- 0
# diff_s[diff_s < 0] <- 0
# diff_sw[diff_sw < 0] <- 0
# diff_w[diff_w < 0] <- 0
# diff_nw[diff_nw < 0] <- 0
# 

slope_delta <- rast(list(diff_n, diff_ne, diff_e, diff_se, diff_s, diff_sw, diff_w, diff_nw))
names(slope_delta) <- c("diff_n", "diff_ne", "diff_e", "diff_se", "diff_s", "diff_sw", "diff_w", "diff_nw")
slope_delta[slope_delta < 0] <- 0


# Sum the differences
flowunits <- sum(slope_delta)

#    The following maps are the percent of the cell's neighbors in a direction
#    that will flow to the cell.  For example, "north", is the percent of the
#    lateral flow out of the cell to the north that will flow to the cell.

percent_flow = slope_delta/flowunits
names(percent_flow) <- c("north", "northeast", "east", "southeast", "south", "southwest", "west", "northwest")


plot(percent_flow)

raster::writeRaster(flowunits, paste0(imitate_smr_setup, "/flowunits.asc"), overwrite=T)

for (direction in names(percent_flow)) {
  file_name <- paste0(imitate_smr_setup, "/", direction, ".asc")
  
  raster::writeRaster(percent_flow[[direction]], file_name, overwrite = T)
}








