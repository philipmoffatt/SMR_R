library(tidyverse)
library(mapview)
library(raster)
library(rasterVis)
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
# mfc mask, buffered outline
mfc_no_pullman <- st_intersection(mfc_no_pullman, MFC_buf)

# national elevation data for the mfc huc
mfc_ned <- raster('/Users/philipmoffatt/Dropbox/SMR_R/processed_data/NED/MFC_NED_UTM.tif')
plot(mfc_ned)
gplot(strm) + 
  geom_raster(aes(fill = value))

# elevation data, convert the raster to a df for ggplot
ned <- as.data.frame(as(mfc_ned, "SpatialPixelsDataFrame"))
colnames(ned) <- c("value", "x", "y")

# get stream cells
strm <- raster('/Users/philipmoffatt/Dropbox/SMR_R/processed_data/dem/MFC/dem_streams.tif')
plot(strm)
strmMsk <- mask(strm, MFC)
creek <- as.data.frame(as(strmMsk, "SpatialPixelsDataFrame"))
colnames(creek) <- c("value", "x", "y")

# set the breaks for contours
mybreaks <- seq(0,1500,50)

coordinate_cities <- data.frame(
  city = c("Pullman", "Moscow"),
  x = c(486700, 499000),
  y = c(5175600, 5175500)) 

# plot the contours with the mfc shape on top; change MFC to mfc_no_pullman for actual extent 
ggplot() +
  geom_contour_filled(data = ned,aes(x, y, z = value), breaks = mybreaks) +
  geom_sf(data = MFC, linewidth = 0.5, color = 'lightyellow', fill = NA) +
  geom_tile(data = creek, aes(x,y),  fill = "lightblue") +
  geom_point(data = coordinate_cities, aes(x = x, y = y), colour = "red") +
  geom_text(data = coordinate_cities, aes(x = x, y = y, label = city), 
            vjust = -1, hjust = 0.15,
            color = 'white') +
  labs(title = 'Missouri Flat Creek',
       x = "",
       y = "",
       fill = 'Elevation (m)')
# add the scale and north arrow

# add map of WA & ID

# add the map of US 


