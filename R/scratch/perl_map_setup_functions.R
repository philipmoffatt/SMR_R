library(tidyverse)
library(mapview)
library(raster)
library(FedData)
library(sf)
library(sp)
library(terra)
library(whitebox)
library(gridExtra)
library(ggplot2)

### This set of functions will use the maps that are pulled using R packages 
##   in the SMR_data_setup.R script and it will convert them to match the format
##   (units, naming, and extent) of the PERL SMR script.

## Variable Setups:
imitate_smr_path <- "/Users/duncanjurayj/Dropbox/SMR_R/processed_data/imitate_smr_setup"
shape_path <- "/Users/duncanjurayj/Dropbox/SMR_R/raw_data/template/mfc_no_pullman/layers/globalwatershed.shp"
basins_path <- "/Users/duncanjurayj/Dropbox/SMR_R/processed_data/dem/MFC/dem_basins.tif"
watershed_id <- 84
dem_breach_path <- "/Users/duncanjurayj/Dropbox/SMR_R/processed_data/dem/MFC/dem_breached.tif"
wshed_mask_path <- "/Users/duncanjurayj/Dropbox/SMR_R/processed_data/imitate_smr_setup/watershed.tif"
streams_raw_path <- "/Users/duncanjurayj/Dropbox/SMR_R/processed_data/dem/MFC/dem_streams.tif"
nlcd_simple_path <- "/Users/duncanjurayj/Dropbox/SMR_R/processed_data/ASC/MFC/NLCD_simple.asc"
mfc_folder_path <- "/Users/duncanjurayj/Dropbox/SMR_R/processed_data/ASC/MFC"
initial_depth_names <- c("depth_A.asc", "depth_B.asc")
output_depth_names <- c("soil_depth_A.tif", "soil_depth_B.tif")
landuse_path <- file.path(imitate_smr_path, "landuse.tif")
stream_path <- file.path(imitate_smr_path, "strms_30m.tif")

wshed_mask <- read_watershed_mask(shape_path, 
                                  basins_path, 
                                  watershed_id, 
                                  buffer=FALSE,
                                  buffer_amount=0)

write_watershed_mask(mask = wshed_mask, 
                     output_folder = imitate_smr_path, 
                     filename = "watershed")

crop_dem_to_watershed(dem_path = dem_breach_path, 
                      watershed_mask_path = wshed_mask_path,
                      output_folder = imitate_smr_path, 
                      file_name = "el.tif")

calculate_flow_direction(dem_path = file.path(imitate_smr_path, "el.tif"), 
                         output_folder = imitate_smr_path,
                         produce_plot = TRUE)

process_streams(raw_streams_path = streams_raw_path, 
                watershed_mask_path = wshed_mask_path,
                output_folder = imitate_smr_path, 
                output_filename = "strms_30m.tif",
                produce_plots = TRUE
                )

process_land_use(landuse_path = nlcd_simple_path, 
                 output_folder = imitate_smr_path,
                 output_filename = "landuse.tif",
                 produce_plots = TRUE,
                 watershed_mask_path = wshed_mask_path)

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


calculate_soil_depth(initial_depth_root_path = mfc_folder_path, 
                     initial_depth_names, 
                     output_depth_names, 
                     landuse_path, 
                     stream_path, 
                     output_dir = imitate_smr_path,
                     watershed_mask_path = wshed_mask_path)

combined_saturation_amount(
  root_folder_path = imitate_smr_path,
  sat_mc_A_name = "sat_mc_A.tif",
  sat_mc_B_name = "sat_mc_B.tif",
  soil_depth_A_name = "soil_depth_A.tif",
  soil_depth_B_name = "soil_depth_B.tif",
  sat_combined_name = "sat_amt.tif",
  watershed_mask_path = wshed_mask_path
)

combined_wiltpt_amount(
  root_folder_path = imitate_smr_path,
  wiltpt_mc_A_name = "wiltpt_mc_A.tif",
  wiltpt_mc_B_name = "wiltpt_mc_B.tif",
  soil_depth_A_name = "soil_depth_A.tif",
  soil_depth_B_name = "soil_depth_B.tif",
  wiltpt_combined_name = "wiltpt_amt.tif",
  watershed_mask_path = wshed_mask_path
)

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


calculate_Ksubsurface(
  root_folder_path = imitate_smr_path,
  Ksat_matrix_A_name = "Ksat_matrix_A.tif",
  output_name = "Ksubsurface.tif",
  watershed_mask_path = wshed_mask_path,
  produce_plots = TRUE
)

calculate_ETreduction_mc(
  root_folder_path = imitate_smr_path,
  fieldcap_amt_name = "fieldcap_amt.tif",
  soil_depth_name = "soil_depth.tif",
  output_name = "ETreduction_mc.tif",
  watershed_mask_path = wshed_mask_path,
  produce_plots = TRUE
)

calculate_Ksat_mpores(
  root_folder_path = imitate_smr_path,
  Ksat_matrix_A_name = "Ksat_matrix_A.tif",
  Ksat_matrix_B_name = "Ksat_matrix_B.tif",
  output_A_name = "Ksat_mpore_A.tif",
  output_B_name = "Ksat_mpore_B.tif",
  watershed_mask_path = wshed_mask_path,
  produce_plots = TRUE
)

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

calculate_Ksubsurface <- function(root_folder_path, Ksat_matrix_A_name, output_name,
                                  watershed_mask_path, produce_plots = TRUE) {
  
  watershed_mask <- raster(watershed_mask_path)
  
  Ksat_matrix_A <- raster(file.path(root_folder_path, Ksat_matrix_A_name))
  Ksat_matrix_A <- mask(Ksat_matrix_A, watershed_mask)
  crs(Ksat_matrix_A) <- crs(watershed_mask)


  Ksubsurface <- Ksat_matrix_A / 500
  
  output_path <- file.path(root_folder_path, output_name)
  raster::writeRaster(Ksubsurface, output_path, overwrite = TRUE)
  
  if (produce_plots) {
    plot(Ksubsurface, main = "Ksubsurface")
  }
}

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
  Ksat_mpore_B <- Ksat_matrix_B / 2
  
  output_path_A <- file.path(root_folder_path, output_A_name)
  output_path_B <- file.path(root_folder_path, output_B_name)
  
  raster::writeRaster(Ksat_mpore_A, output_path_A, overwrite = TRUE)
  raster::writeRaster(Ksat_mpore_B, output_path_B, overwrite = TRUE)
  
  if (produce_plots) {
    plot(Ksat_mpore_A, main = "Ksat_mpore_A")
    plot(Ksat_mpore_B, main = "Ksat_mpore_B")
  }
}

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





