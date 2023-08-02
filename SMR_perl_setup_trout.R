library(tidyverse)
library(mapview)
library(raster)
library(FedData)
library(sf)
library(terra)
library(whitebox)
library(rasterVis)
library(sp)
library(gridExtra)

# issues 7/26/23
# some na values in the surgo data under compname:  Water and Rock land
# an artifact in the Northeastn portion of the watershed related to soils data
# need to verify the O horzion fix and handle others
# low conductivity, very low Kfc values
# does the iso basins effect modeling outputs?
# error for map() lin 216

# obtain watershed shapefile


## define path to shape file ----
raw_path = "./raw_data"


# obtain crs from the HUC12 shape file
columbia <- st_read(file.path(raw_path, 'WBD_17_HU2_shape/Shape/WBDHU12.shp'))

# mapview(columbia) # used the mapview to find the object id for the trout basin
trout_HUC12 <- columbia %>% filter(ObjectID == "4751")

# # get the shapes from the streamstats downloads. convienient for the pour points.
# barrett <- st_read(file.path(raw_path,'barrett/globalwatershed.shp'))
# b <- mapview(barrett, color = "blue")
# mires <- st_read(file.path(raw_path,'mires/globalwatershed.shp'))
# m <- mapview(mires, color = 'green')
# trout <- st_read(file.path(raw_path,'trout/globalwatershed.shp'))
# t <- mapview(trout, color = 'red')
# # the three creek discharges to trout lake
# b+m+t

# transform the crs for the crappy stream stats shape from the columbia 
trout <- 
  st_read(
    file.path(
      raw_path, 
      'trout/globalwatershed.shp')) %>%
  #st_buffer(., 300) %>%
  st_transform(crs = st_crs(trout_HUC12)) 

mapview(trout)

# Adding the buffer causes the extent to include the lake which translates to weird streams in the whitebox process.

### pulling data 
NED <- FedData::get_ned(
  template = trout,
  res = 1,
  extraction.dir = "./processed_data/NED",
  label = "trout", 
  force.redo = F)

ssurgo <- FedData::get_ssurgo(
  template = trout,
  raw.dir = "./raw_data/SSURGO",
  extraction.dir = "./processed_data/SSURGO",
  label = "trout", 
  force.redo = F)

NLCD <- FedData::get_nlcd(
  template = trout,
  year = 2019,
  dataset = "landcover",
  # raw.dir = "./raw_data/NLCD",
  extraction.dir = "./processed_data/NLCD",
  label = "trout", 
  force.redo = F)


### white box DEM processing -----------------
whitebox::wbt_init()

NED_utm = projectRaster(from = NED,
                        crs = 'EPSG:32611', #
                        res = 30, 
                        filename = gsub(raster::filename(NED),
                                        pattern = "_1",
                                        replacement = "_UTM"),
                        overwrite = TRUE)

dem <- NED_utm@file@name

# looking at the projection
# NED_utm_plt <- raster(dem)
# gplot(NED_utm_plt) +
#   geom_tile(aes(fill = value))


# build DEM directory 
outpath="./processed_data/dem"
if(!dir.exists(outpath)) dir.create(outpath)

# Build wshed specific out path directory 
dir_path <- file.path(outpath,"trout")
  # make directory per site 
if(!dir.exists(dir_path)) dir.create(dir_path)
  

wbt_feature_preserving_smoothing(
    dem = dem,
    output = paste0(dir_path, "/dem_smoothed.tif"),
    filter = 9)

# dem_smooth <- raster("./processed_data/dem/trout/dem_smoothed.tif")
# gplot(dem_smooth) +
#   geom_tile(aes(fill = value))

wbt_breach_depressions(dem = paste0(dir_path, "/dem_smoothed.tif"),
                         output = paste0(dir_path,"/dem_breached.tif"))

wbt_fill_depressions(dem = paste0(dir_path, "/dem_smoothed.tif"),
                       output = paste0(dir_path,"/dem_fill.tif"))

wbt_d8_pointer(dem = paste0(dir_path, "/dem_breached.tif"),
               output = paste0(dir_path,"/dem_d8_pointer.tif"))

d8 <- raster(paste0(dir_path,"/dem_d8_pointer.tif"))
plot(d8)

#preset basin size builds trout - need to check pour point defaults 
wbt_basins(d8_pntr =paste0(dir_path,"/dem_d8_pointer.tif"),
           output = paste0(dir_path,"/dem_basins.tif")
          )

basins <- raster(paste0(dir_path,"/dem_basins.tif"))
plot(basins)

# build a boolean mask of trout 
# this needs to be genralized. for now, examine basins plot and guess/cheack until you find a value that results in the correct basin shape
trout_mask = basins == 58 # guess
trout_mask[is.na(trout_mask[])] = 0 # grey out all other basins
plot(trout_mask) # check

# trout_mask
raster::writeRaster(trout_mask, paste0(dir_path,"/dem_trout_mask.tif"), overwrite=T)


wbt_d8_flow_accumulation(input = paste0(dir_path,"/dem_breached.tif"),
                                  output = paste0(dir_path,"/d8_flow_accum.tif"))

wbt_extract_streams(flow_accum = paste0(dir_path,"/d8_flow_accum.tif"),
                    output = paste0(dir_path,"/dem_streams.tif"),
                    threshold = 1000) #

## explore the stream data 
streams <- raster(paste0(dir_path,"/dem_streams.tif"))
test <- (streams*trout_mask)
mapview(streams)
strm <- raster(paste0(dir_path,"/dem_streams.tif"))
mapview(log(strm))

# build small basins for potential parellel computing. 
wbt_isobasins(dem = paste0(dir_path,"/dem_breached.tif"),
              output = paste0(dir_path,"/iso_basin_1000.tif"), 
              size = 30000)

### set to same threshold as streams and we can run all basins in parallel fast on kamiak maybe
iso = raster(paste0(dir_path,"/iso_basin_1000.tif"))
mapview(iso) + mapview(trout_mask)


# Prep terrain maps, NEED dem(pit filled), ASPECT, FLOWDISTIONS(check needed unit),
# STREAMS, region, and watersheds.
# all mask to smallest extent and output as .asc
#
dem_files <- list.files(dir_path, pattern = "tif",full.names = T)
#
asc_path = dir_path %>% gsub(., pattern = "dem/trout", replacement = "ASC")
trout_path = dir_path %>% gsub(., pattern = "dem/", replacement = "ASC/")
#
if(!dir.exists(asc_path)){dir.create(asc_path)}
if(!dir.exists(trout_path)){dir.create(trout_path)}
#
#
# # Pass out all ASC
#
# # trout_mask
test <- mask(NED_utm, trout_mask, maskvalue = 1, incervse = T)
#
map(dem_files, function(h){
  print(h)
  temp <- raster(h) %>% mask(., trout_mask, maskvalue = 1, inverse = T)
  writeRaster(temp, filename = gsub(
    gsub(h,pattern = "dem/", replacement = "ASC/"),
    pattern = ".tif", replacement = ".asc"),
    overwrite=T)
})

# testing <- raster("./processed_data/dem/trout/d8_flow_accum.tif")
# res(testing)
# testing <- projectRaster(testing, crs = crs(trout_mask), res = res(trout_mask))
# testing <- crop(testing, trout_mask)
# testing <- mask(testing, trout_mask, maskvalue = 1, inverse = TRUE)

# Function to process and write each raster file
process_and_write_raster <- function(h) {
  print(h)
  # Read the raster
  temp <- raster(h)
  # Reproject the raster to match the projection and resolution of the mask
  temp <- projectRaster(temp, crs = crs(trout_mask), res = res(trout_mask))
  ############################### new line 7/27/23 plm##
  temp <- crop(temp, trout_mask)
  # Mask the raster using the trout_mask
  temp <- mask(temp, trout_mask, maskvalue = 1, inverse = TRUE)
  # Write the masked raster to a new file
  writeRaster(temp, filename = gsub(
    gsub(h, pattern = "dem/", replacement = "ASC/"),
    pattern = ".tif", replacement = ".asc"),
    overwrite = TRUE
  )
}

# Apply the function to each raster file
map(dem_files, process_and_write_raster)

# loop for sub watersheds -- do this later
#
iso_mask = mask(iso, trout_mask, maskvalue = 1, inverse = T)
#

# testing the function below
# unique(iso_mask)
# testing <- raster("./processed_data/dem/trout/dem_streams.tif") %>%
#   mask(., trout_mask, maskvalue = 1, inverse = T) %>%
#   mask(., iso_mask==89, maskvalue = 1, inverse = T)
# mapview(testing)

 for(i in unique(iso_mask)){
   temp_path = trout_path %>% gsub(., pattern = "trout", replacement = paste0("sub_",i))
   if(!dir.exists(temp_path)){dir.create(temp_path)}

   map(dem_files[1:(length(dem_files)-2)], function(h){
     print(h)
     temp <- raster(h) %>% mask(., trout_mask, maskvalue = 1, inverse = T) %>%
       mask(., iso_mask==i, maskvalue = 1, inverse = T)

     #trim_NA, check shared extent and that it makes sense for post processing

     print(
       gsub(
         gsub(h,pattern = "dem/trout", replacement = paste0("ASC/sub_",i,"")),
         pattern = ".tif", replacement = ".asc")
       )

     temp <- temp %>% trim()

    #save out
    writeRaster(temp, filename = gsub(
       gsub(h,pattern = "dem/trout", replacement = paste0("ASC/sub_",i,"")),
       pattern = ".tif", replacement = ".asc"),
       overwrite=T)
   })

}

### soils-----------------------------------------

# weighting soil characteristics by majority group
majority_comps = ssurgo$tabular$component %>% group_by(mukey) %>% 
  slice(which.max(comppct.r))

# check no dropped or duplicated map unit keys (mukeys)
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

# for definitions on the column names in surgo see below
# https://data.nal.usda.gov/system/files/SSURGO_Metadata_-_Table_Column_Descriptions.pdf#:~:text=Column%20Physical%20Name%3A%20rvindicatorColumn%20Label%3A%20RV%3F%20A%20yes%2Fno,%28set%20of%20values%29%20is%20representative%20for%20the%20component.

# can pull soil texture percentages if needed for pedotransfer functions
# could pull H and L for calibration range on all soil parameters though often missing 

#pull majority component of each mukey 
hoz_data <- ssurgo$tabular$chorizon[
  ssurgo$tabular$chorizon$cokey %in% majority_comps$cokey,] %>% 
  # as.data.frame() %>% 
  dplyr::select(all_of(hoz_col))

# rebuild Hoz thickness
hoz_data <- hoz_data %>% mutate(hzthk.r = hzdepb.r - hzdept.r)

# pull restriction data for trout 
res_data <- ssurgo$tabular$corestrictions[
  ssurgo$tabular$corestrictions$cokey %in% majority_comps$cokey,] %>% 
  dplyr::select(c("reskind","resdept.l","reshard","cokey"))
   # wbt_slope(dem = paste0(dir_path,"/dem_smoothed.tif"),
    #             output = paste0(dir_path,"/slope.tif"),
    #             units = "degrees")
    #   
# join hoz_data and res_data and treat depth (first use weighted average 
# layer A = 0-30cm or max/res if very shallow, 
# layer B =  30 - max or restriction), if B = 0, set to 0.1 cm

# check if res depth is less than depth a 
res_data$resdept.l %>% table
# 0  15  20  51  56 102 
# 1   2   3  30   1   4 
res_data$resdept.l[21] <- 1 # changed depth to restriction to 1 from 0

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

### build the list of soil characterisitcs in trout in layers A and B
for(i in unique(hoz_data$cokey)){
  layer_A_temp_2 <- list()
  layer_B_temp_2 <- list()
  
  temp_hoz = hoz_data %>% filter(cokey == i)
  temp_res = res_data %>% filter(cokey == i)
  
  # ! if there is no restriction in the component 
  # find the weighted average 
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
# trout_spatial_A = left_join(ssurgo$spatial %>% 
#                           rename("mukey" = "MUKEY"), 
#                         majority_comps %>% 
#                           select("compname", "mukey" , 'cokey') %>% 
#                           mutate(mukey = as.character(mukey))) %>% 
#   left_join(., layer %>% filter(layer=="A"))
# 
# trout_spatial_B = left_join(ssurgo$spatial %>% 
#                             rename("mukey" = "MUKEY"), 
#                           majority_comps %>% 
#                             select("compname", "mukey" , 'cokey') %>% 
#                             mutate(mukey = as.character(mukey))) %>% 
#   left_join(., layer %>% filter(layer=="B"))


## join by each layer (1 cokey = 2 layers + attributes)
trout_spatial = left_join(ssurgo$spatial %>% 
                            rename("mukey" = "MUKEY"), 
                          majority_comps %>% 
                            dplyr::select("compname", "mukey" , 'cokey') %>% 
                            mutate(mukey = as.character(mukey))) %>% 
  full_join(., layer)

# QA check on land/type ####
trout_spatial_cln <- trout_spatial %>% # soil table without water and rock
  filter(compname != 'Water', compname != 'Rock land')

# values for both layers A and B of land identified as Rock outcrop in the soil survey. 
om_rock <- 0.01
dbovendry_rock <- 0.01  # weight of water lost from moist soil from oven drying
dbthirdbar_rock <- 0.01 # The oven dry weight of the less than 2 mm soil material per unit volume of soil at a water tension of 1/3 bar.
ksat_rock <- 0.00001 # The amount of water that would move vertically through a unit area of saturated soil in unit time under unit hydraulic gradient.
wsatiated_rock <- 1 # The volumetric content of soil water retained at a tension of 15 bars (1500 kPa)
wthirdbar_rock <-  0.5 # water content near field capacity
wfifteenbar_rock <- 0.5 # water content at wilt point
# depth <- assigned values of 1 cm for A and 0.1 for B

trout_spatial_rock <- trout_spatial %>%
  filter(compname == 'Rock land') %>%
  mutate(
    om.r = om_rock,
    dbovendry.r = dbovendry_rock,
    dbthirdbar.r = dbthirdbar_rock,
    ksat.r = ksat_rock,
    wsatiated.r = wsatiated_rock,
    wthirdbar.r = wthirdbar_rock,
    wfifteenbar.r = wfifteenbar_rock
  )

# values for both layers A and B of land identified as water in the soil survey. Spatially, two ponds were located in the middle of Borosaprists soils, so the soil characteristics were transfered to the water cells to replace the NAs. 
om_water <- 20
dbovendry_water <- 0.3  # weight of water lost from moist soil from oven drying
dbthirdbar_water <- 0.21 # The oven dry weight of the less than 2 mm soil material per unit volume of soil at a water tension of 1/3 bar.
ksat_water <- 5 # The amount of water that would move vertically through a unit area of saturated soil in unit time under unit hydraulic gradient.
wsatiated_water <- 30 # The estimated volumetric soil water content at or near zero bar tension, expressed as a percentage of the whole soil.
wthirdbar_water <-  26 # water content near field capacity
wfifteenbar_water <- 11 # water content at wilt point
depth_A <- 2
depth_B <- 206

trout_spatial_water <- trout_spatial %>%
  filter(compname == 'Water') %>%
  mutate(
    om.r = om_water,
    dbovendry.r = dbovendry_water,
    dbthirdbar.r = dbthirdbar_water,
    ksat.r = ksat_water,
    wsatiated.r = wsatiated_water,
    wthirdbar.r = wthirdbar_water,
    wfifteenbar.r = wfifteenbar_water,
    layer = c('A', 'B'),
    depth = c(depth_A, depth_B)
    )
# Using rbind to combine the data frames
trout_spatial_clnd <- rbind(trout_spatial_cln, trout_spatial_rock, trout_spatial_water)


# build out rasters for each trout and ISO watershed, ####
# **Check extent between hydro and soil asc

for(i in unique(trout_spatial_clnd$layer)){
  for(j in names(trout_spatial_clnd)[8:15]){
    temp = terra::rasterize(trout_spatial_clnd %>% 
                              filter(layer==i) %>% 
                              vect() %>% 
                              project(., crs(NED_utm)), 
                            NED_utm %>% rast,
                            field = j)
    
    
    #build trout layers 
    temp_trout <-raster(temp) %>% mask(., trout_mask, maskvalue = 1, inverse = T) 
   
     writeRaster(temp_trout, filename = paste0(
       "./processed_data/ASC/trout/",j,"_",i,".asc"), 
      overwrite=T)
  }
}

### land cover ----- 
# reclassify to group grasses, forests, barren cover types. SMR uses land cover to simulate canopy interception and storage
# unique(NLCD)
# [1] 11 21 22 23 24 41 42 43 52 71 81 90 95
# https://www.mrlc.gov/data/legends/national-land-cover-database-class-legend-and-description
NLCD_reclass_ma =matrix(data = c(
  11, 1, #water
  21, 2, #urban/other
  22, 2, #urban/other
  23, 2, #urban/other
  24, 2, #urban/other
  41, 3, # deciduous forest
  42, 3, # evergreen forest
  43, 3, # mixed forest
  52, 4, #shrub
  71, 5, #grass
  81, 5, #grass/pasture
  90, 3, #woody wetlands
  95, 5), # grassy wetland
  ncol =2, byrow = T)

NLCD_reclass <- reclassify(NLCD, NLCD_reclass_ma) %>% projectRaster(., NED_utm, method = "ngb")
temp_trout <- NLCD_reclass %>% mask(., trout_mask, maskvalue = 1, inverse = T) 

writeRaster(temp_trout, filename = paste0(
  "./processed_data/ASC/trout/NLCD_simple.asc"), 
  overwrite=T)

## clip all thefiles at once 
## loop for sub watersheds -- do this later 
#  load all trout asc 

trout_asc_files = list.files(path = "./processed_data/ASC/trout/", full.names = T)
tst <- rast(trout_asc_files)
iso_mask = mask(iso, trout_mask, maskvalue = 1, inverse = T)

#for(k in unique(iso_mask)){
#    temp_path = trout_path %>% gsub(., pattern = "trout", replacement = paste0("sub_",k))
#    if(!dir.exists(temp_path)){dir.create(temp_path)}

#      temp <- tst %>% crop(., rast(iso_mask)) %>% 
#        mask(., rast(iso_mask)==k, maskvalue = 1, inverse = T) %>% trim
#      print(temp_path)
#      writeRaster(temp_trout, filename = paste0(temp_path,"/",j,"_",i,".asc"),
#        overwrite=T)
  
#  }

# re imagine the building of iso build full set then load as a stack and mask all at once
#

print("The asc file package complete")

### This set of functions will use the maps that are pulled using R packages --------------
##   in the SMR_data_setup.R script and it will convert them to match the format
##   (units, naming, and extent) of the PERL SMR script.

## Variable Setups:
imitate_smr_path <- "./processed_data/imitate_smr_setup"
shape_path <- "./raw_data/trout/globalwatershed.shp"
#WBD_17_HU2_Shape/Shape/WBDHU12.shp"
basins_path <- "./processed_data/dem/trout/dem_basins.tif"
watershed_id <- 58 # unique id from line 134 SMR_data_setup_trout
dem_breach_path <- "./processed_data/dem/trout/dem_breached.tif"
wshed_mask_path <- "./processed_data/imitate_smr_setup/watershed.tif"
streams_raw_path <- "./processed_data/dem/trout/dem_streams.tif"
nlcd_simple_path <- "./processed_data/ASC/trout/NLCD_simple.asc"
mfc_folder_path <- "./processed_data/ASC/trout"
initial_depth_names <- c("depth_A.asc", "depth_B.asc")
output_depth_names <- c("soil_depth_A.tif", "soil_depth_B.tif")
landuse_path <- file.path(imitate_smr_path, "landuse.tif")
stream_path <- file.path(imitate_smr_path, "strms_30m.tif")


### Watershed Mask(s):
##   This section reads in the mask, crop it if necessary, and buffer it if
##   necessary. Important to do this first and have it be consistent so all 
##   future maps align. 
read_watershed_mask <- function(shapefile_path, basins, watershed_id, buffer=FALSE, buffer_amount = 0) {
  
  watershed <- st_read(shapefile_path)
  basins <- raster(basins)
  basins[basins != watershed_id] <- NA
  
  watershed <- st_transform(watershed, crs = crs(basins))
  
  watershed_mask <- mask(basins, watershed)
  
  if (buffer) {
    watershed_mask <- buffer(wshed_mask, buffer_amount)
  }
  
  watershed_mask[watershed_mask == watershed_id] <- 1
  
  plot(watershed_mask, main = "Watershed Mask")
  
  
  return(watershed_mask)
}

write_watershed_mask <- function(mask, output_folder, filename="watershed") {
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  output_file <- file.path(output_folder, paste0(filename, ".tif"))
  
  writeRaster(mask, output_file, overwrite = TRUE)
  
  return(output_file)
}


wshed_mask <- read_watershed_mask(shape_path, 
                                  basins_path, 
                                  watershed_id, 
                                  buffer=FALSE,
                                  buffer_amount=0)

write_watershed_mask(mask = wshed_mask, 
                     output_folder = imitate_smr_path, 
                     filename = "watershed")

# set values for watershed properties file that initiates smr
watershed <- raster("./processed_data/imitate_smr_setup/watershed.tif")
watershed[watershed==0] <- NA
mapview(watershed)

wshed_id = 	4751 # huc 12 ObjectID
area_cells = freq(watershed, value = 1) # number of cells in the watershed mask
res_vol = 50 # the initial volume (cm) assumed in the linear reservoir
res_coeff = 0.05 # the recession rate from the linear reservoir to baseflow

# build frame from values
wshed_frame = data.frame(wshed_id, area_cells, res_vol, res_coeff)

# write out .ini file to imitation folder and scratch folder
write.table(wshed_frame, sep = " ", row.names = F, col.names = F, file = paste0(imitate_smr_path, "/wshed_res_properties.ini"))
write.table(wshed_frame, sep = " ", row.names = F, col.names = F, file = "./R/scratch/wshed_res_properties.ini") # this is for initial perl script development

### Elevation:
##   This section uses the DEM pulled using the FedData package in 
##   SMR_data_setup.R, crops it to the watershed mask, and renames it to match
##   The PERL SMR format. 
crop_dem_to_watershed <- function(dem_path, watershed_mask_path, output_folder, file_name="el.tif") {
  
  # Read DEM and watershed mask
  dem <- raster(dem_path)
  watershed_mask <- raster(watershed_mask_path)
  
  # Project DEM to the CRS of the watershed mask
  dem <- projectRaster(dem, crs = crs(watershed_mask))
  
  # Crop DEM to the extent of the watershed mask
  cropped_dem <- crop(dem, extent(watershed_mask))
  
  # Mask the DEM using the watershed mask
  masked_dem <- mask(cropped_dem, watershed_mask)
  
  # Generate the output file path
  output_file <- file.path(output_folder, file_name)
  
  # Write the cropped and masked DEM to the output file
  writeRaster(masked_dem, output_file, overwrite = TRUE)
  
  # Return the output file path
  return(output_file)
}


crop_dem_to_watershed(dem_path = dem_breach_path, 
                      watershed_mask_path = wshed_mask_path,
                      output_folder = imitate_smr_path, 
                      file_name = "el.tif")

### Flow Units:
##   This section uses the watershed mask and the digital elevation model to
##   produce a set of maps with flow units (cardinal directions).
calculate_flow_direction <- function(dem_path, output_folder, 
                                     produce_plot = TRUE) {
  
  el <- terra::rast(dem_path)
  
  diff_n <- el - (terra::shift(el, 0, 30) %>% terra::extend(., el) %>% terra::crop(., el))
  diff_ne <- (el - (terra::shift(el, 30, 30) %>% terra::extend(., el) %>% terra::crop(., el))) * 0.707
  diff_e <- el - (terra::shift(el, 30, 0) %>% terra::extend(., el) %>% terra::crop(., el))
  diff_se <- (el - (terra::shift(el, 30, -30) %>% terra::extend(., el) %>% terra::crop(., el))) * 0.707
  diff_s <- el - (terra::shift(el, 0, -30) %>% terra::extend(., el) %>% terra::crop(., el))
  diff_sw <- (el - (terra::shift(el, -30, -30) %>% terra::extend(., el) %>% terra::crop(., el))) * 0.707
  diff_w <- el - (terra::shift(el, -30, 0) %>% terra::extend(., el) %>% terra::crop(., el))
  diff_nw <- (el - (terra::shift(el, -30, 30) %>% terra::extend(., el) %>% terra::crop(., el))) * 0.707
  
  slope_delta <- terra::rast(list(diff_n, diff_ne, diff_e, diff_se, diff_s, diff_sw, diff_w, diff_nw))
  names(slope_delta) <- c("north", "northeast", "east", "southeast", "south", "southwest", "west", "northwest")
  slope_delta[slope_delta < 0] <- 0
  
  flowunits <- sum(slope_delta)
  
  percent_flow <- slope_delta / flowunits
  names(percent_flow) <- c("north", "northeast", "east", "southeast", "south", "southwest", "west", "northwest")
  
  for (direction in names(percent_flow)) {
    file_name <- paste0(output_folder, "/", direction, ".tif")
    terra::writeRaster(percent_flow[[direction]], file_name, overwrite = TRUE)
  }
  
  if (produce_plot) {
    for (direction in names(percent_flow)) {
      plot(percent_flow[[direction]], main = direction)
    }
  }
  
  flowunits_file <- paste0(output_folder, "/flowunits.tif")
  terra::writeRaster(flowunits, flowunits_file, overwrite = TRUE)
  
  if (produce_plot) {
    plot(flowunits, main = "Flow Units")
  }
  
  return(flowunits_file)
}

calculate_flow_direction(dem_path = file.path(imitate_smr_path, "el.tif"), 
                         output_folder = imitate_smr_path,
                         produce_plot = TRUE)


### Streams:
##   This section ensures that the streams map (derived in SMR_data_setup.R) 
##   is cropped to the extent of the watershed mask. It will also make sure that
##   non-stream cells are values of 0 and not NA. This will protect other maps
##   from decaying in the PERL SMR script. 
process_streams <- function(raw_streams_path, watershed_mask_path, output_folder, output_filename, produce_plots = TRUE) {
  raw_streams <- raster(raw_streams_path)
  watershed_mask <- raster(watershed_mask_path)
  
  raw_streams <- projectRaster(raw_streams, crs = crs(watershed_mask))
  raw_streams <- crop(raw_streams, watershed_mask)
  raw_streams[is.na(raw_streams)] <- 0
  raw_streams[is.na(watershed_mask)] <- NA
  
  output_path <- file.path(output_folder, output_filename)
  writeRaster(raw_streams, output_path, overwrite = TRUE)
  
  if (produce_plots) {
    plot(raw_streams, main = "Processed Stream Raster")
  }
}

process_streams(raw_streams_path = streams_raw_path, 
                watershed_mask_path = wshed_mask_path,
                output_folder = imitate_smr_path, 
                output_filename = "strms_30m.tif",
                produce_plots = TRUE
)

### Land Use:
##   This section uses the land use map produced in SMR_data_setup.R
##   (land use is pulled from NLDC with FedData package). It crops the land use
##   map to the current watershed mask and ensures that the values are written
##   out as integers and not floats (integers are needed for the Boolean logic
##   used in PERL SMR). 
process_land_use <- function(landuse_path, output_folder, output_filename, produce_plots = TRUE, watershed_mask_path = NULL) {
  raster_data <- raster(landuse_path)
  
  if (!is.null(watershed_mask_path)) {
    watershed_mask <- raster(watershed_mask_path)
    raster_data <- mask(raster_data, watershed_mask)
    
  }
  
  output_path <- file.path(output_folder, output_filename)
  
  writeRaster(raster_data, output_path, datatype = "INT4S", overwrite = TRUE)
  
  if (produce_plots) {
    plot(raster_data, main = "Processed Raster")
  }
  
}

process_land_use(landuse_path = nlcd_simple_path, 
                 output_folder = imitate_smr_path,
                 output_filename = "landuse.tif",
                 produce_plots = TRUE,
                 watershed_mask_path = wshed_mask_path)

### Percentages to Moisture Contents:
percent_to_mc <- function(percentage_path, watershed_mask_path, output_folder, output_filename, produce_plots=TRUE) {
  percentage_raster <- raster(percentage_path)
  
  watershed_mask <- raster(watershed_mask_path)
  
  cropped_raster <- mask(percentage_raster, watershed_mask)
  
  converted_raster <- cropped_raster / 100
  
  output_path <- file.path(output_folder, output_filename)
  
  writeRaster(converted_raster, output_path, overwrite = TRUE)
  
  if (produce_plots) {
    plot_title <- basename(output_filename)
    plot(converted_raster, main = plot_title)
  }
  
  return(output_path)
}


percent_to_mc(percentage_path = file.path(mfc_folder_path, "/wfifteenbar.r_A.asc"),
              watershed_mask_path = wshed_mask_path,
              output_folder = imitate_smr_path,
              output_filename = "/wiltpt_mc_A.tif",
              produce_plots = TRUE)

percent_to_mc(percentage_path = file.path(mfc_folder_path, "/wfifteenbar.r_B.asc"),
              watershed_mask_path = wshed_mask_path,
              output_folder = imitate_smr_path,
              output_filename = "/wiltpt_mc_B.tif",
              produce_plots = TRUE)

percent_to_mc(percentage_path = file.path(mfc_folder_path, "/wthirdbar.r_A.asc"),
              watershed_mask_path = wshed_mask_path,
              output_folder = imitate_smr_path,
              output_filename = "/fieldcap_mc_A.tif",
              produce_plots = TRUE)

percent_to_mc(percentage_path = file.path(mfc_folder_path, "/wthirdbar.r_B.asc"),
              watershed_mask_path = wshed_mask_path,
              output_folder = imitate_smr_path,
              output_filename = "/fieldcap_mc_B.tif",
              produce_plots = TRUE)

percent_to_mc(percentage_path = file.path(mfc_folder_path, "/wsatiated.r_A.asc"),
              watershed_mask_path = wshed_mask_path,
              output_folder = imitate_smr_path,
              output_filename = "/sat_mc_A.tif",
              produce_plots = TRUE)

percent_to_mc(percentage_path = file.path(mfc_folder_path, "/wsatiated.r_B.asc"),
              watershed_mask_path = wshed_mask_path,
              output_folder = imitate_smr_path,
              output_filename = "/sat_mc_B.tif",
              produce_plots = TRUE)

### Conductivity Units:
##   This section uses the maps related to conductivity (produced in 
##   SMR_data_setup.R). The main tasks are to rename the maps, convert units 
##   where it is necessary, and crop the maps to the current watershed mask.
um_per_second_to_cm_per_day <- function(um_per_second_path, watershed_mask_path, output_folder, output_filename, produce_plots = TRUE) {
  # Read the um per second raster
  um_per_second_raster <- raster(um_per_second_path)
  
  # Read the watershed mask raster
  watershed_mask <- raster(watershed_mask_path)
  
  # Mask the um per second raster using the watershed mask
  masked_raster <- mask(um_per_second_raster, watershed_mask)
  
  # Convert units from um/second to cm/day (*8.64)
  converted_raster <- masked_raster * 8.64
  
  # Generate the output file path
  output_path <- file.path(output_folder, output_filename)
  
  # Write the converted raster to the output file
  writeRaster(converted_raster, output_path, overwrite = TRUE)
  
  # Create a plot of the converted raster if produce_plots is TRUE
  if (produce_plots) {
    plot_title <- basename(output_filename)
    plot(converted_raster, main = plot_title)
  }
  
  # Return the output file path
  return(output_path)
}

um_per_second_to_cm_per_day(
  um_per_second_path = file.path(mfc_folder_path, "/ksat.r_A.asc"),
  watershed_mask_path = wshed_mask_path,
  output_folder = imitate_smr_path, output_filename = "/Ksat_matrix_A.tif",
  produce_plots = TRUE)

um_per_second_to_cm_per_day(
  um_per_second_path = file.path(mfc_folder_path, "/ksat.r_B.asc"),
  watershed_mask_path = wshed_mask_path,
  output_folder = imitate_smr_path, output_filename = "/Ksat_matrix_B.tif",
  produce_plots = TRUE)

### Soil Depth:
##   This section renames the soil depth rasters produced in SMR_data_setup.R,
##   produces a combined depth raster, and crops the soil depth rasters to the 
##   current watershed mask.
calculate_soil_depth <- function(initial_depth_root_path, initial_depth_names, output_depth_names, landuse_path, stream_path, watershed_mask_path, output_dir, produce_plots = TRUE) {
  watershed_mask <- raster(watershed_mask_path)
  
  for (x in 1:length(initial_depth_names)) {
    horizon <- ifelse(grepl("A", initial_depth_names[x]), "A", ifelse(grepl("B", initial_depth_names[x]), "B", NA))
    print(paste("Soil Horizon:", horizon))
    
    landuse <- raster(landuse_path)
    depth_rast <- raster(file.path(initial_depth_root_path, initial_depth_names[x]))
    streams <- raster(stream_path)
    
    # Mask soil depth raster with watershed mask
    crs(depth_rast) <- crs(watershed_mask)
    depth_rast <- mask(depth_rast, watershed_mask)
    
    if (produce_plots) {
      plot(depth_rast, main = paste("Soil Depth in", horizon, "Horizon (no stream effect)"))
    }
    
    if (horizon == "A") {
      depth_rast_no_streams <- depth_rast
      depth_rast[streams > 0] <- 1
      depth_rast[landuse == 2] <- 2
    } else if (horizon == "B") {
      depth_rast_no_streams <- depth_rast
      depth_rast[streams > 0] <- 19
      depth_rast[landuse == 2] <- 20
    } else {
      depth_rast <- depth_rast
    }
    
    if (produce_plots) {
      plot(depth_rast, main = paste("Soil Depth in", horizon, "Horizon (with stream effect)"))
    }
    
    output_path <- file.path(output_dir, output_depth_names[x])
    
    print(paste("Writing out renamed raster:", output_path))
    raster::writeRaster(depth_rast, output_path, overwrite = TRUE)
  }
  
  # Combine soil depth rasters (A + B)
  combined_depth <- raster::stack(file.path(output_dir, output_depth_names))
  combined_depth <- sum(combined_depth, na.rm=FALSE)
  
  combined_output_path <- file.path(output_dir, "soil_depth.tif")
  print(paste("Writing out combined soil depth raster:", combined_output_path))
  raster::writeRaster(combined_depth, combined_output_path, overwrite = TRUE)
}

calculate_soil_depth(initial_depth_root_path = mfc_folder_path, 
                     initial_depth_names, 
                     output_depth_names, 
                     landuse_path, 
                     stream_path, 
                     output_dir = imitate_smr_path,
                     watershed_mask_path = wshed_mask_path)


## Final Modifications:
##   There are modifications or additions that can be made only once the
##   moisture content and soil depth maps are finished. This section produces
##   these final maps, names them correctly, and crops them to the watershed
##   mask.
combined_saturation_amount <- function(root_folder_path, sat_mc_A_name,
                                       soil_depth_A_name, sat_mc_B_name,
                                       soil_depth_B_name, sat_combined_name,
                                       watershed_mask_path, produce_plots = TRUE) {
  
  sat_mc_A <- raster(file.path(root_folder_path, sat_mc_A_name))
  soil_depth_A <- raster(file.path(root_folder_path, soil_depth_A_name))
  sat_mc_B <- raster(file.path(root_folder_path, sat_mc_B_name))
  soil_depth_B <- raster(file.path(root_folder_path, soil_depth_B_name))
  
  watershed_mask <- raster(watershed_mask_path)
  
  # Assign CRS of output rasters to match watershed mask
  crs(sat_mc_A) <- crs(watershed_mask)
  crs(soil_depth_A) <- crs(watershed_mask)
  crs(sat_mc_B) <- crs(watershed_mask)
  crs(soil_depth_B) <- crs(watershed_mask)
  
  # Mask rasters with watershed mask
  sat_mc_A <- mask(sat_mc_A, watershed_mask)
  soil_depth_A <- mask(soil_depth_A, watershed_mask)
  sat_mc_B <- mask(sat_mc_B, watershed_mask)
  soil_depth_B <- mask(soil_depth_B, watershed_mask)
  
  combined_sat_amt <- (sat_mc_A * soil_depth_A) + (sat_mc_B * soil_depth_B)
  
  output_path <- file.path(root_folder_path, sat_combined_name)
  raster::writeRaster(combined_sat_amt, output_path, overwrite = TRUE)
  
  if (produce_plots) {
    plot(combined_sat_amt, main = "Combined Saturation Amount")
  }
}

combined_saturation_amount(
  root_folder_path = imitate_smr_path,
  sat_mc_A_name = "sat_mc_A.tif",
  sat_mc_B_name = "sat_mc_B.tif",
  soil_depth_A_name = "soil_depth_A.tif",
  soil_depth_B_name = "soil_depth_B.tif",
  sat_combined_name = "sat_amt.tif",
  watershed_mask_path = wshed_mask_path
)

combined_wiltpt_amount <- function(root_folder_path, wiltpt_mc_A_name,
                                   wiltpt_mc_B_name, soil_depth_A_name,
                                   soil_depth_B_name, wiltpt_combined_name,
                                   watershed_mask_path, produce_plots = TRUE) {
  
  wiltpt_mc_A <- raster(file.path(root_folder_path, wiltpt_mc_A_name))
  soil_depth_A <- raster(file.path(root_folder_path, soil_depth_A_name))
  wiltpt_mc_B <- raster(file.path(root_folder_path, wiltpt_mc_B_name))
  soil_depth_B <- raster(file.path(root_folder_path, soil_depth_B_name))
  
  watershed_mask <- raster(watershed_mask_path)
  
  # Assign CRS of output rasters to match watershed mask
  crs(wiltpt_mc_A) <- crs(watershed_mask)
  crs(soil_depth_A) <- crs(watershed_mask)
  crs(wiltpt_mc_B) <- crs(watershed_mask)
  crs(soil_depth_B) <- crs(watershed_mask)
  
  # Mask rasters with watershed mask
  wiltpt_mc_A <- mask(wiltpt_mc_A, watershed_mask)
  soil_depth_A <- mask(soil_depth_A, watershed_mask)
  wiltpt_mc_B <- mask(wiltpt_mc_B, watershed_mask)
  soil_depth_B <- mask(soil_depth_B, watershed_mask)
  
  combined_wilt_amt <- (wiltpt_mc_A * soil_depth_A) + (wiltpt_mc_B * soil_depth_B)
  
  output_path <- file.path(root_folder_path, wiltpt_combined_name)
  raster::writeRaster(combined_wilt_amt, output_path, overwrite = TRUE)
  
  if (produce_plots) {
    plot(combined_wilt_amt, main = "Combined Wilt Point Amount")
  }
}

combined_wiltpt_amount(
  root_folder_path = imitate_smr_path,
  wiltpt_mc_A_name = "wiltpt_mc_A.tif",
  wiltpt_mc_B_name = "wiltpt_mc_B.tif",
  soil_depth_A_name = "soil_depth_A.tif",
  soil_depth_B_name = "soil_depth_B.tif",
  wiltpt_combined_name = "wiltpt_amt.tif",
  watershed_mask_path = wshed_mask_path
)

calculate_fieldcap_amount <- function(root_folder_path, fieldcap_mc_A_name,
                                      soil_depth_A_name, fieldcap_mc_B_name,
                                      soil_depth_B_name, output_dir,
                                      residual_A, residual_B,
                                      watershed_mask_path,
                                      produce_plots = TRUE) {
  
  watershed_mask <- raster(watershed_mask_path)
  
  fieldcap_mc_A <- raster(file.path(root_folder_path, fieldcap_mc_A_name))
  soil_depth_A <- raster(file.path(root_folder_path, soil_depth_A_name))
  fieldcap_mc_B <- raster(file.path(root_folder_path, fieldcap_mc_B_name))
  soil_depth_B <- raster(file.path(root_folder_path, soil_depth_B_name))
  
  fieldcap_amt_A <- (fieldcap_mc_A * soil_depth_A) - (soil_depth_A * residual_A)
  fieldcap_amt_B <- (fieldcap_mc_B * soil_depth_B) - (soil_depth_B * residual_B)
  
  crs(fieldcap_amt_A) <- crs(watershed_mask)
  crs(fieldcap_amt_B) <- crs(watershed_mask)
  
  fieldcap_amt_A <- mask(fieldcap_amt_A, watershed_mask)
  fieldcap_amt_B <- mask(fieldcap_amt_B, watershed_mask)
  
  
  fieldcap_amt <- fieldcap_amt_A + fieldcap_amt_B
  
  output_path_A <- file.path(output_dir, "fieldcap_amt_A.tif")
  output_path_B <- file.path(output_dir, "fieldcap_amt_B.tif")
  output_path_combined <- file.path(output_dir, "fieldcap_amt.tif")
  
  raster::writeRaster(fieldcap_amt_A, output_path_A, overwrite = TRUE)
  raster::writeRaster(fieldcap_amt_B, output_path_B, overwrite = TRUE)
  raster::writeRaster(fieldcap_amt, output_path_combined, overwrite = TRUE)
  
  if (produce_plots) {
    plot(fieldcap_amt_A, main = "Field Capacity Amount A")
    plot(fieldcap_amt_B, main = "Field Capacity Amount B")
    plot(fieldcap_amt, main = "Combined Field Capacity Amount")
  }
}

calculate_fieldcap_amount(
  root_folder_path = imitate_smr_path,
  fieldcap_mc_A_name = "fieldcap_mc_A.tif",
  fieldcap_mc_B_name = "fieldcap_mc_B.tif",
  soil_depth_A_name = "soil_depth_A.tif",
  soil_depth_B_name = "soil_depth_B.tif",
  output_dir = imitate_smr_path,
  residual_A = 0.02,
  residual_B = 0.02,
  watershed_mask_path = wshed_mask_path
)

# The vertical conductivity of the restrictive layer is set to 0.1 cm/day here. This controls the rate of percolation from storage into the restrictive layer.
calculate_Ksubsurface <- function(root_folder_path, Ksat_matrix_A_name, output_name,
                                  watershed_mask_path, produce_plots = TRUE) {
  
  watershed_mask <- raster(watershed_mask_path)
  
  Ksat_matrix_A <- raster(file.path(root_folder_path, Ksat_matrix_A_name))
  Ksat_matrix_A <- mask(Ksat_matrix_A, watershed_mask)
  crs(Ksat_matrix_A) <- crs(watershed_mask)
  
  
  Ksubsurface <- Ksat_matrix_A / Ksat_matrix_A # a value of 1
  Ksubsurface <- Ksubsurface / 10              # a final value of 0.1 cm/day, pretty high
  
  output_path <- file.path(root_folder_path, output_name)
  raster::writeRaster(Ksubsurface, output_path, overwrite = TRUE)
  
  if (produce_plots) {
    plot(Ksubsurface, main = "Ksubsurface")
  }
}

calculate_Ksubsurface(
  root_folder_path = imitate_smr_path,
  Ksat_matrix_A_name = "Ksat_matrix_A.tif",
  output_name = "Ksubsurface.tif",
  watershed_mask_path = wshed_mask_path,
  produce_plots = TRUE
)

calculate_ETreduction_mc <- function(root_folder_path, fieldcap_amt_name, soil_depth_name,
                                     output_name, watershed_mask_path, produce_plots = TRUE) {
  
  watershed_mask <- raster(watershed_mask_path)
  
  fieldcap_amt <- raster(file.path(root_folder_path, fieldcap_amt_name))
  crs(fieldcap_amt) <- crs(watershed_mask)
  fieldcap_amt <- mask(fieldcap_amt, watershed_mask)
  
  soil_depth <- raster(file.path(root_folder_path, soil_depth_name))
  soil_depth <- mask(soil_depth, watershed_mask)
  crs(soil_depth) <- crs(watershed_mask)
  
  ETreduction_mc <- fieldcap_amt * 0.8 / soil_depth
  
  output_path <- file.path(root_folder_path, output_name)
  raster::writeRaster(ETreduction_mc, output_path, overwrite = TRUE)
  
  if (produce_plots) {
    plot(ETreduction_mc, main = "ETreduction_mc")
  }
}

calculate_ETreduction_mc(
  root_folder_path = imitate_smr_path,
  fieldcap_amt_name = "fieldcap_amt.tif",
  soil_depth_name = "soil_depth.tif",
  output_name = "ETreduction_mc.tif",
  watershed_mask_path = wshed_mask_path,
  produce_plots = TRUE
)

calculate_Ksat_mpores <- function(root_folder_path, Ksat_matrix_A_name, Ksat_matrix_B_name,
                                  output_A_name, output_B_name, watershed_mask_path, produce_plots = TRUE) {
  
  watershed_mask <- raster(watershed_mask_path)
  
  Ksat_matrix_A <- raster(file.path(root_folder_path, Ksat_matrix_A_name))
  crs(Ksat_matrix_A) <- crs(watershed_mask)
  Ksat_matrix_A <- crop(Ksat_matrix_A, watershed_mask)
  Ksat_matrix_A <- mask(Ksat_matrix_A, watershed_mask)
  
  Ksat_matrix_B <- raster(file.path(root_folder_path, Ksat_matrix_B_name))
  crs(Ksat_matrix_B) <- crs(watershed_mask)
  Ksat_matrix_B <- crop(Ksat_matrix_B, watershed_mask)
  Ksat_matrix_B <- mask(Ksat_matrix_B, watershed_mask)
  
  Ksat_mpore_A <- Ksat_matrix_A * 10
  Ksat_mpore_B <- Ksat_matrix_B * 2. # this was "/" for some reason?
  
  output_path_A <- file.path(root_folder_path, output_A_name)
  output_path_B <- file.path(root_folder_path, output_B_name)
  
  raster::writeRaster(Ksat_mpore_A, output_path_A, overwrite = TRUE)
  raster::writeRaster(Ksat_mpore_B, output_path_B, overwrite = TRUE)
  
  if (produce_plots) {
    plot(Ksat_mpore_A, main = "Ksat_mpore_A")
    plot(Ksat_mpore_B, main = "Ksat_mpore_B")
  }
}

calculate_Ksat_mpores(
  root_folder_path = imitate_smr_path,
  Ksat_matrix_A_name = "Ksat_matrix_A.tif",
  Ksat_matrix_B_name = "Ksat_matrix_B.tif",
  output_A_name = "Ksat_mpore_A.tif",
  output_B_name = "Ksat_mpore_B.tif",
  watershed_mask_path = wshed_mask_path,
  produce_plots = TRUE
)

calculate_Kfc <- function(root_folder_path, Ksat_matrix_A_name, Ksat_matrix_B_name,
                          sat_mc_A_name, sat_mc_B_name, fieldcap_amt_A_name, fieldcap_amt_B_name,
                          soil_depth_A_name, soil_depth_B_name, output_A_name, output_B_name,
                          watershed_mask_path, produce_plots = TRUE) {
  
  watershed_mask <- raster(watershed_mask_path)
  
  Ksat_matrix_A <- raster(file.path(root_folder_path, Ksat_matrix_A_name))
  crs(Ksat_matrix_A) <- crs(watershed_mask)
  Ksat_matrix_A <- crop(Ksat_matrix_A, watershed_mask)
  Ksat_matrix_A <- mask(Ksat_matrix_A, watershed_mask)
  
  Ksat_matrix_B <- raster(file.path(root_folder_path, Ksat_matrix_B_name))
  crs(Ksat_matrix_B) <- crs(watershed_mask)
  Ksat_matrix_B <- crop(Ksat_matrix_B, watershed_mask)
  Ksat_matrix_B <- mask(Ksat_matrix_B, watershed_mask)
  
  sat_mc_A <- raster(file.path(root_folder_path, sat_mc_A_name))
  sat_mc_B <- raster(file.path(root_folder_path, sat_mc_B_name))
  fieldcap_amt_A <- raster(file.path(root_folder_path, fieldcap_amt_A_name))
  fieldcap_amt_B <- raster(file.path(root_folder_path, fieldcap_amt_B_name))
  soil_depth_A <- raster(file.path(root_folder_path, soil_depth_A_name))
  soil_depth_B <- raster(file.path(root_folder_path, soil_depth_B_name))
  
  Kfc_A <- Ksat_matrix_A * exp((-13.0 / sat_mc_A) *
                                 (sat_mc_A - fieldcap_amt_A / soil_depth_A))
  
  Kfc_B <- Ksat_matrix_B * exp((-13.0 / sat_mc_B) *
                                 (sat_mc_B - fieldcap_amt_B / soil_depth_B))
  
  output_path_A <- file.path(root_folder_path, output_A_name)
  output_path_B <- file.path(root_folder_path, output_B_name)
  
  raster::writeRaster(Kfc_A, output_path_A, overwrite = TRUE)
  raster::writeRaster(Kfc_B, output_path_B, overwrite = TRUE)
  
  if (produce_plots) {
    plot(Kfc_A, main = "Kfc_A")
    plot(Kfc_B, main = "Kfc_B")
  }
}

calculate_Kfc(
  root_folder_path = imitate_smr_path,
  Ksat_matrix_A_name = "Ksat_matrix_A.tif",
  Ksat_matrix_B_name = "Ksat_matrix_B.tif",
  sat_mc_A_name = "sat_mc_A.tif",
  sat_mc_B_name = "sat_mc_B.tif",
  fieldcap_amt_A_name = "fieldcap_amt_A.tif",
  fieldcap_amt_B_name = "fieldcap_amt_B.tif",
  soil_depth_A_name = "soil_depth_A.tif",
  soil_depth_B_name = "soil_depth_B.tif",
  output_A_name = "Kfc_A.tif",
  output_B_name = "Kfc_B.tif",
  watershed_mask_path = wshed_mask_path,
  produce_plots = TRUE
)

