library(tidyverse)
library(mapview)
library(raster)
library(FedData)
library(sf)

## define MFC ----
raw_path = "./raw_data"

# grab dummy data with a covering extent: PRW raster. 
PRW <- st_read(file.path(raw_path, "template", "huc12.shp"))
  #   mapview(PRW)

# get MFC and buffer by 300 m or 10 cells to ensure we 
  #   capture the eatershed divide 
MFC <- PRW %>% filter(Id == "79")

MFC_buf <- PRW %>% filter(Id == "79") %>% st_buffer(., 300)
  #   mapview(MFC)


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
# plot(mfc_rast)
raster::writeRaster(mfc_rast, paste0(dir_path,"/dem_mfc_mask.tif"))


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
    weights = hzthk.r/sum(hzthk.r)
  )
  
  data_sum <- (data %>% 
                 dplyr::select(all_of(cols))*data$weights) %>% 
    summarise_all( sum) %>% 
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
  
  layer_A_temp_2 = depth_weight(data =layer_A_temp, cols = hoz_comp,
                                layer = "A")
  rm(layer_A_temp)
  
  if(max_depth>30){
    layer_B_temp = temp_hoz %>% 
      filter(hzdepb.r > depth_a)%>% 
      filter(hzdept.r <= max_depth)%>% 
      arrange(hzdept.r)
    
    layer_B_temp$hzdept.r[1] <- depth_a
    layer_B_temp$hzdepb.r[nrow(layer_B_temp)] <- max_depth
    
    layer_B_temp_2 = depth_weight(layer_B_temp, cols = hoz_comp, layer = "B")
  }
  
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

NLCD_reclass <- reclassify(NLCD, NLCD_reclass_ma) %>% projectRaster(., NED_utm)
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

for(k in unique(iso_mask)){
    temp_path = mfc_path %>% gsub(., pattern = "MFC", replacement = paste0("sub_",k))
    if(!dir.exists(temp_path)){dir.create(temp_path)}

      temp <- tst %>% crop(., rast(iso_mask)) %>% 
        mask(., rast(iso_mask)==k, maskvalue = 1, inverse = T) %>% trim

      writeRaster(temp_mfc, filename = paste0(temp_path,"/",j,"_",i,".asc"),
        overwrite=T)

  }

# re imagine the building of iso build full set then load as a stack and mask all at once
#

print("The asc file package complete")

