# import libraries
library(tidyverse)
library(dplyr) 
library(lubridate)
library(AOI)
library(climateR)
library(rnoaa)
library(chillR)
library(humidity)
library(bigleaf)
library(mapview)
library(sf)
library(FAO56)
library(elevatr)

### File Description:
#     This file provides functions for downloading gridMET data, producing all the 
#     daily inputs required for running the SMR model, and formatting a CSV file
#     to be used by the SMR model in PERL. 
# --------------------------------------------------------------------------- #


# Purpose: Produce a data frame of gridMET weather data using the climateR package.
# Parameters: 
#   site_name: character -- The name of the site to pull data from (e.g., "Pullman, WA")
# Returns:
#   A data frame with:
#     prcp - cm
#     tmax - C
#     tmin - C
#     pet_grass - cm
#     srad - W/m^2
#     wind_vel - m/s
#     vpd - kPa
# Notes: 
#   The gridMET data is in mm and the function is hard coded to convert 
#   values to cm for PERL SMR. 
pull_gridMET <- function(site_name) {
  
  # establish AOI
  AOI::aoi_get(site_name) %>% 
    AOI::aoi_map(returnMap = TRUE)
  
  # pull lat and lon for site
  cords <- geocode(site_name)
  lat <- cords$lat
  lon <- cords$lon

  gridMETdata <- param_meta$gridmet
  
  # pull precipitation (mm), temp (C), pet (mm) from GRIDMET
  AOI_p = AOI::geocode(site_name, pt = TRUE)
  gridMET_data  = as.data.frame(getGridMET(AOI_p, param = c('prcp', 'tmax', 'tmin', 'pet_grass', 'srad', 'wind_vel', 'vpd'), 
                                           startDate = "2000-01-01", 
                                           endDate = "2022-12-31"))
  
  mapview(AOI_p)
  
  gridMET_data <- gridMET_data %>% 
    mutate(year = lubridate::year(date), 
           month = lubridate::month(date), 
           day = lubridate::day(date),
           tmax = tmax - 273.15,
           tmin = tmin - 273.15,
           tavg = (tmax+tmin)/2,
           doy = yday(date),
           l_turb = 2.5,
           prcp = prcp/10,
           pet_grass = pet_grass/10
           )
  
  return(gridMET_data)
  
}

# Purpose: Makes 24 hourly temperatures out of a daily maximum and minimum temperatures using the chillR package.
# Parameters: 
#   gridMET: The data frame returned by pull_gridMET().
#   latitude: The latitude of the measurement site (gridMET tile).
# Returns:
#   gridMET data frame with hourly temperature columns (24 new columns). 
# Notes:
get_hourly_temps <- function(gridMET, latitude) {
  
  gridMET_hourly <- gridMET[,c("tmin", "tmax", "year", "month", "day")] %>%
    rename(Year = year, Day = day, Month = month, Tmin = tmin, Tmax = tmax) %>%
    make_hourly_temps(latitude, year_file = ., keep_sunrise_sunset = FALSE) %>%
    rename(year = Year, day = Day, month = Month, tmin = Tmin, tmax = Tmax)
  
  # Merge the two data.frames
  gridMET_hourly <- merge(gridMET, gridMET_hourly)
  
  return(gridMET_hourly)
}

# Purpose: Produces a crop coefficient curve. 
# Parameters:
#   gridMET: The gridMET data frame.
#   time_params: A list containing every growth stage as the name and every corresponding growth length as the value.
#   Kc_params: A list containing every growth stage as the name and every corresponding crop coefficient as the value.
#   planting_offset: A numeric that causes the crop coefficient curve to be shifted to a different time of year. 
#                    generally this is set to 0. 
#   crop_name: A character for the name of the crop. This will determine the column name (e.g. "Wheat" --> Kc_Wheat)
# Returns:
#   The gridMET data frame with a crop coefficient curve added to it. The new columns 
#   will be named based on the crop_name.
# Notes:
#   The offset parameter should usually be set to 0. It has only been used so far 
#   for making the winter wheat curve where planting and growth occurs at a 
#   significantly different time of year. 
ann_Kc_curve <- function(gridMET, time_params, Kc_params, planting_offset, crop_name) {
  
  total_curve <- c(rep(Kc_params[["Kc_preplant"]], time_params[["preplant_len"]]),
                   rep(Kc_params[["Kc_init"]], time_params[["init_len"]]),
                   seq(Kc_params[["Kc_init"]], Kc_params[["Kc_mid"]], length.out = time_params[["dev_len"]]),
                   rep(Kc_params[["Kc_mid"]], time_params[["mid_len"]]),
                   seq(Kc_params[["Kc_mid"]], Kc_params[["Kc_end"]], length.out = time_params[["late_len"]]),
                   rep(Kc_params[["Kc_post_harv"]], time_params[["post_len"]]))
  
  # checking that all sections add up to a year
  stopifnot(length(total_curve) == sum(unlist(time_params)))
  
  doy <- as.numeric(1:365)
  Kc <- round(total_curve, 4)
  Kc_frame <- data.frame(doy = (doy + planting_offset) %% 365, Kc = Kc)
  
  # attach the Kc curve to gridMET
  col_name <- paste0("Kc_", crop_name)
  gridMET <- inner_join(gridMET, Kc_frame, by = c("doy" = "doy"))
  colnames(gridMET)[colnames(gridMET) == "Kc"] <- col_name
  
  return(gridMET)
}

# Purpose: Calculates the daily dew point temperature.
# Parameters:
#   gridMET: A data frame with daily weather data from gridMET.
# Returns:
#   The updated gridMET data frame with a column called "tdew".
# Notes:
#   This function assumes the following columns exist in the gridMET data frame:
#   tmin.
dewpoint_temperature <- function(gridMET) {
  
  # create a new column for dew point temperature
  gridMET$tdew <- case_when(gridMET$tmin <= 0 ~ gridMET$tmin * 0.7 - 2,
                            gridMET$tmin > 0 ~ gridMET$tmin * 0.97 - 0.53)
  
  return(gridMET)
}

# Purpose: Calculates the daily cloud fraction and adds it to the gridMET data frame.
# Parameters:
#   gridMET: A data frame with daily weather data from gridMET.
#   latitude: A numeric class with the latitude for the location of the gridMET data.
#   longitude: A numeric class with the longitude for the location of the gridMET data.
#   timezone: A numeric value representing the timezone offset in hours from Greenwich (I think but worth double checking)
# Returns:
#   the gridMET data frame with a new "cloud" column containing daily cloud fraction values.
# Notes:
calculate_cloud_fraction <- function(gridMET, latitude, longitude, timezone) {
  
  # setup sequence of days for potential radiation function
  julian_year <- seq(1, 365)
  
  # calculate the sunrise, sunset, and daytime hours for each day based on day of year
  # and latitude
  annual_day_lengths <- daylength(latitude = latitude, JDay = julian_year, notimes.as.na = FALSE)
  
  # setup solar radiation data frame for merging with gridMET frames
  srad_df <- 
    data.frame(
      doy = numeric(), 
      srad_potential = numeric()
    )
  
  for (day in julian_year) {
    
    # get the sunrise and sunset time for the corresponding day of year and build
    # a sequence from them --> the 0.1 step is arbitrary
    day_sequence <- 
      seq(
        from = annual_day_lengths$Sunrise[day], 
        to = annual_day_lengths$Sunset[day],
        by = 0.1
      )
    
    # calculate the potential radiation for the sequence
    potential_radiation <- 
      potential.radiation(
        doy = day,
        hour = day_sequence,
        latDeg = latitude,
        longDeg = longitude,
        timezone = timezone,
        useSolartime = TRUE
      )
    
    # take the average potential radiation to match a daily summary
    mean_daily_srad <- mean(potential_radiation)
    
    # add average to corresponding position in the potential_srad dataframe
    srad_df[day, "doy"] <- day
    srad_df[day, "srad_potential"] <- mean_daily_srad
    
  }
  
  # add the repeating potential solar radiation frame to the gridMET dataframe
  # based on doy
  gridMET_joined <- inner_join(gridMET, srad_df[, c("doy", "srad_potential")], by = "doy")
  
  # cloud fraction is the ratio of the observed solar radiation to 
  # potential solar radiation at that point
  gridMET_joined$cloud <- gridMET_joined$srad / gridMET_joined$srad_potential
  
  return(gridMET_joined)
}

# Purpose: Produces a daily atmospheric pressure value using the bigleaf package.
# Parameters:
#   gridMET: A data frame with daily weather data from gridMET. 
#   site_latitude: A numeric value specifying the latitude of the gridMET data.
#   site_longitude: A numeric values specifying the longitude of the gridMET data.
# Returns:
#   The gridMET data frame with a new column called "atmospheric pressure".
# Notes:
#   This didn't actually get used in the end but I'm leaving it in the code. 
pressure_from_el <- function(gridMET, site_latitude, site_longitude) {
  
  # formatting latitude and longitude inputs to match parameter requirements
  # for R package 'elevatr'
  latitude_longitude_frame <- data.frame(x = site_longitude, y = site_latitude)
  
  site_elevation_spatial <- 
    get_elev_point(
      locations = latitude_longitude_frame,
      prj = "EPSG:4326", # this could be made more flexible in the future
      src = 'aws', # not passing 'aws' or using 'epqs' produces 'lexical error'
    )
  
  site_elevation_value <- site_elevation_spatial$elevation # in meters
  
  # using R package 'bigleaf' function pressure.from.elevation()
  # Apply pressure.from.elevation() for each row in gridMET dataframe
  gridMET_updated <- gridMET %>%
    mutate(
      atmospheric_pressure = pressure.from.elevation(
        elev = site_elevation_value,
        Tair = tavg,
        VPD = vpd,
        constants = bigleaf.constants() # built in constants
      )
    )
  
  return(gridMET_updated)
}

# Purpose: Estimates daily wind speed at a specified height based on the measured wind speed at a different height.
# Parameters: 
#   actual_measurement_height: The actual measurement height of the gridMET wind speed data.
#   measured_wind_speed: Daily measured wind speed data.
#   measured_land_cover_height: The height of the land cover from Hansen 1993.
#   predicted_land_cover_height: The height of the new land cover.
# Returns:
#   The gridMET data frame with 
# Notes:
#   Didn't end up getting used but leaving it in for now. This is a helper 
#   function in land_cover_rh_wrapper() function in this file.
estimated_wind_speed <- function(
    actual_measurement_height = 2,
    measured_wind_speed,
    measured_land_cover_height = 0.9, # from Hansen 1993, not a big change
    predicted_land_cover_height # possibly a confusing parameter name (the wind 
    #measurement height not the land cover height is predicted though 
    # measurement height is defined by this land cover height)
    ) {
  
  measured_zero_plane_displacement <- 0.65 * measured_land_cover_height
  measured_momentum_roughness <- 0.1 * measured_land_cover_height
  
  predicted_wind_speed <- 
    measured_wind_speed * 
    log( (predicted_land_cover_height - measured_zero_plane_displacement) / measured_momentum_roughness) / 
    log( (measured_land_cover_height - measured_zero_plane_displacement) / measured_momentum_roughness)
  
  return(predicted_wind_speed)
}


# Purpose: Calculates daily heat transfer to resistance (rh) values.
# Parameters:
#   wind_measurement_height: A numeric value specifying the height of the wind measurement in meters.
#   temp_measurement_height: A numeric value specifying the height of the temperature measurement in meters.
#   von_karman: Von Karman's constant.
#   zero_plane_displacement_height: A numeric value specifying the zero plane displacemnet height of the land cover.
#   momentum_roughness: A numeric value specifying the momentum roughness parameter of the land cover.
#   heat_vapor_roughness: A numeric value specifying the heat vapor roughness of the land cover. 
#   wind_speed: A numeric value for the daily wind speed in meters per second. 
# Returns:
#    A single resitance to heat transfer value to be placed in a row in the new rh column.
# Notes:
#   This is a helper function used in the land_cover_rh() function and 
#   it is applied to an entire wind speed column in the land_cover_rh_wrapper() function.
heat_transfer_resistance <- 
  function(
    wind_measurement_height, 
    temp_measurement_height, 
    von_karman = 0.41, 
    zero_plane_displacement_height, 
    momentum_roughness, 
    heat_vapor_roughness, 
    wind_speed
    ) {
  
  rh <- (
    log((wind_measurement_height - zero_plane_displacement_height + momentum_roughness) / momentum_roughness) * 
    log((temp_measurement_height - zero_plane_displacement_height + heat_vapor_roughness) / heat_vapor_roughness)) / 
    (von_karman^2 * wind_speed) * (1/86400)
  
  return(rh)
  }

# Purpose: Calculates RH values for a specific land cover.
# Parameters:
#   land_cover_height: A numeric specifying the height in meters of the land cover.
#   wind_measurement_height: A numeric specifying the height in meters of the wind measurement.
#   temp_measurement_height: A numeric specifying the height in meters 
#   von_karman: A numeric specifying the Von Karman constant. 
#   wind_speed: A numeric specifying the wind speed in meters per second.
# Returns:
#   A single RH value.
# Notes:
land_cover_rh <- 
  function(
    land_cover_height, 
    wind_measurement_height, 
    temp_measurement_height,
    von_karman = 0.41,
    wind_speed
  ) {
    
    
    # the resistance to heat transfer function only works when:
    # the height of the wind speed measurement is greater than 75% of the land
    # cover height. 
    if ((land_cover_height*0.75) >= wind_measurement_height) {
      
      displacement_height <- 0.65 * land_cover_height
      momentum_roughness <- 0.1 * land_cover_height
      heat_vapor_roughness <- 0.2 * momentum_roughness

      estimated_wind_speed_val <- estimated_wind_speed(
        measured_wind_speed = wind_speed,
        predicted_land_cover_height = land_cover_height)
      
      new_wind_measurement_height <- land_cover_height
      new_temp_measurement_height <- land_cover_height
      
      rh_value <- heat_transfer_resistance(
        wind_measurement_height = new_wind_measurement_height,
        temp_measurement_height = new_temp_measurement_height,
        von_karman = von_karman,
        zero_plane_displacement_height = displacement_height,
        momentum_roughness = momentum_roughness,
        heat_vapor_roughness = heat_vapor_roughness,
        wind_speed = estimated_wind_speed_val
      )
      
    } else {
      
      displacement_height <- 0.65 * land_cover_height
      momentum_roughness <- 0.1 * land_cover_height
      heat_vapor_roughness <- 0.2 * momentum_roughness
      
      rh_value <- heat_transfer_resistance(
        wind_measurement_height = wind_measurement_height,
        temp_measurement_height = temp_measurement_height,
        von_karman = von_karman,
        zero_plane_displacement_height = displacement_height,
        momentum_roughness = momentum_roughness,
        heat_vapor_roughness = heat_vapor_roughness,
        wind_speed = wind_speed
      )
    }
    
    return(rh_value)
  }

# Purpose: Produces a new column in the gridMET data frame containing daily RH 
#   values for each land cover of interest.
# Parameters:
#   gridMET: A data frame containing daily weather data and SMR input data from gridMET.
#   land_cover_names: A character vector of land cover names.
#   land_cover_heights: A numeric vector of land cover heights in the same order as the land cover names.
#   wind_measurement_height: A numeric specifying the height of the wind measurement speed.
#   temp_measurement_height: A numeric specifying the height of the temperature measurement speed.
#   von_karman: A numeric specifying the Von Karman constant.
# Returns:
#   The gridMET data frame with 1-n new daily rh columns for each land cover.
# Notes:
land_cover_rh_wrapper <- 
  function(
    gridMET,
    land_cover_names,
    land_cover_heights,
    wind_measurement_height = 2,
    temp_measurement_height = 2,
    von_karman = 0.41
  ) {
    
    for (i in 1:length(land_cover_names)) {
      land_cover_name <- land_cover_names[i]
      land_cover_height <- land_cover_heights[i]
      
      # Use sapply instead of lapply, and fix the anonymous function syntax
      rh_values <- sapply(gridMET$wind_vel, function(wind_vel) {
        rh_val <- land_cover_rh(
          land_cover_height = land_cover_height,
          wind_measurement_height = wind_measurement_height,
          temp_measurement_height = temp_measurement_height,
          von_karman = von_karman,
          wind_speed = wind_vel
        )
        return(rh_val)
      })
      
      # Add a new column to gridMET with name 'rh_[land_coverage_name]'
      gridMET[paste0("rh_", land_cover_name)] <- rh_values
      
    }
    
    return(gridMET)
  }

# Purpose: Renames the columns of the gridMET data frame.
# Parameters:
#   gridMET: A data frame containing daily weather data and SMR input data from gridMET.
#   renamed_columns: A character vector containing the new names for the gridMET data frame. 
#                    Names should be in the same order as the current gridMET columns.
# Returns:
#   The gridMET data frame with new column names.
# Notes:
#   This is a helper function to the subset_order_columns() function. 
rename_columns <- function(gridMET, renamed_columns) {
  
  colnames(gridMET) <- renamed_columns
  
  return(gridMET)
  
}

# Purpose: Subset the gridMET data frame to only contain a specified set of columns and renames 
#          existing columns.
# Parameters:
#   gridMET: A data frame containing daily weather data and input SMR data from gridMET.
#   desired_columns: A character vector of the columns from the gridMET data frame to keep.
#   renamed_columns: A character vector of the new names for the subset of columns to keep.
# Returns:
#   The gridMET data fame with filtered, ordered, and renamed columns to match the PERL SMR format.
# Notes:
#   The sub-setting occurs BEFORE the renaming.  This means your renamed_columns
#   vector should only contain new names corresponding to the desired_columns vector.
#   The point of renaming is to ensure that the input data matches the expected names
#   in the PERL SMR script and that columns are in the correct order.
subset_order_columns <- function(gridMET, desired_columns, renamed_columns) {
  gridMET_subset <- gridMET %>%
    dplyr::select(all_of(desired_columns)) %>%
    rename_columns(renamed_columns = renamed_columns)
  return(gridMET_subset)
}


# Purpose: Writes out the gridMET for a specific period to the raw_data/weather folder.
# Parameters:
#   gridMET: A data frame containing daily weather data from gridMET and other daily input data to SMR.
#   start_date: A character class date in the yyyy-mm-dd format for the beginning of the modeling period.
#   end_date: A character class date in the yyyy-mm-dd format for the end of the modeling period.
#   location: A character class of the name of the location (e.g. "Pullman", "Moscow", etc.)
#             to make uniquely named output CSV files.
# Returns:
#   No returns. Writes out data.
# Notes:
#   I have been using the location parameter to make different sizes of input data be named 
#   based on their size. For example location = "large", location = "medium", location = "small"
write_weather_data <- function(gridMET, start_date, end_date, location) {
  
  gridMET <- gridMET[order(gridMET$date), ]
  
  # Create the directory if it doesn't exist
  dir_path <- "./raw_data/weather"
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  # Select only rows within the specified date range
  gridMET_clipped <- gridMET[as.Date(gridMET$date) >= as.Date(start_date) &
                                             as.Date(gridMET$date) <= as.Date(end_date), ]
  
  # Write the clipped data to the file
  output_file <- file.path(dir_path, paste0("gridMET_", location, ".csv"))
  write.table(gridMET_clipped,
              file = output_file,
              col.names = FALSE,
              row.names = FALSE,
              sep = " ")
  
}


## constant setup:
#initialize variables
lat <- 46.7312700
lon <- -117.1861


# setting up dictionaries for crop curves
water_time <- list('preplant_len'=60, 'init_len'=60, 'dev_len'=60, 'mid_len'=60, 'late_len'=60, 'post_len'=65)
water_Kc <- list('Kc_preplant'=0.30, 'Kc_init'=0.6525, 'Kc_mid'=0.6525, 'Kc_end'=0.6525, 'Kc_post_harv'=0.30)

urban_time <- list('preplant_len'=60, 'init_len'=60, 'dev_len'=60, 'mid_len'=60, 'late_len'=60, 'post_len'=65)
urban_Kc <- list('Kc_preplant'=0.30, 'Kc_init'=0.30, 'Kc_mid'=0.30, 'Kc_end'=0.30, 'Kc_post_harv'=0.30)

# using FAO conifer values for now meaning it's static but that can change
forest_time <- list('preplant_len'=60, 'init_len'=60, 'dev_len'=60, 'mid_len'=60, 'late_len'=60, 'post_len'=65)
forest_Kc <- list('Kc_preplant'=1, 'Kc_init'=1, 'Kc_mid'=1, 'Kc_end'=1, 'Kc_post_harv'=1)

# using deciduous orchard
shrub_time <- list('preplant_len'=60, 'init_len'=20, 'dev_len'=70, 'mid_len'=90, 'late_len'=30, 'post_len'=95)
shrub_Kc <- list('Kc_preplant'=0.30, 'Kc_init'=0.45, 'Kc_mid'=0.95, 'Kc_end'=0.70, 'Kc_post_harv'=0.30)

#relied on frost dates from AgWeatherNet: Pullman station: frost on --> May 1st, frost off --> September 19th --> splitting resiudal (125) in half for now
grass_time <- list('preplant_len'=113, 'init_len'=10, 'dev_len'=20, 'mid_len'=64, 'late_len'=62, 'post_len'=96)
grass_Kc <- list('Kc_preplant'=0.30, 'Kc_init'=0.30, 'Kc_mid'=0.75, 'Kc_end'=0.75, 'Kc_post_harv'=0.30)

# using spring wheat values for now
#row_crop_time <- list('preplant_len'=75, 'init_len'=20, 'dev_len'=25, 'mid_len'=60, 'late_len'=30, 'post_len'=155)
row_crop_time <- list('preplant_len'=0, 'init_len'=160, 'dev_len'=75, 'mid_len'=75, 'late_len'=25, 'post_len'=30)
row_crop_Kc <- list('Kc_preplant'=0.30, 'Kc_init'=0.7, 'Kc_mid'=1.15, 'Kc_end'=0.35, 'Kc_post_harv'=0.30)

#land_cover_names <- c('water', 'urban', 'forest', 'shrub', 'grass', 'row_crop')
#land_cover_heights <- c(0.01, 4, 11, 0.3, 0.9, 0.9)
land_cover_names <- c('snow')
land_cover_heights <- c(0.01)


# column names
desired_columns <- c("date", "year", "month", "day", "doy", "tmax", "tmin",
                     "tavg", "tdew", "prcp","pet_grass","Hour_1","Hour_6","Hour_12","Hour_18",
                     "l_turb","cloud","Kc_water","Kc_urban","Kc_forest","Kc_shrub",
                     "Kc_grass", "Kc_row_crop","rh_snow")


smr_columns <- c("date", "year", "month", "day", "doy", "tmax", "tmin",
                 "tavg", "tdew", "precip","pet","hour_1","hour_6","hour_12","hour_18",
                 "l_turb","cloud","cc_water","cc_urban","cc_forest","cc_shrub",
                 "cc_grass","cc_row_crop","rh_snow")


## calling functions
gm_m <- pull_gridMET('Moscow, ID') 
gm_p <- pull_gridMET('Pullman, WA') %>%
  pressure_from_el(site_latitude = lat, site_longitude = lon) %>%
  get_hourly_temps(lat) %>%
  ann_Kc_curve(water_time, water_Kc, 0, 'water') %>%
  ann_Kc_curve(urban_time, urban_Kc, 0, 'urban') %>%
  ann_Kc_curve(forest_time, forest_Kc, 0, 'forest') %>%
  ann_Kc_curve(shrub_time, shrub_Kc, 0, 'shrub') %>%
  ann_Kc_curve(grass_time, grass_Kc, 0, 'grass') %>%
  ann_Kc_curve(row_crop_time, row_crop_Kc, 274, 'row_crop') %>%
  dewpoint_temperature() %>%
  calculate_cloud_fraction(lat, lon, -7) %>% 
  land_cover_rh_wrapper(land_cover_names = land_cover_names,
                        land_cover_heights = land_cover_heights,
                        wind_measurement_height = 2,
                        temp_measurement_height = 2,
                        von_karman = 0.41) %>%
  subset_order_columns(desired_columns = desired_columns,
                       renamed_columns = smr_columns)

start_date <- "2018-05-23"
end_date <- "2022-06-10"

# write big historical
write_weather_data(
  gridMET = gm_p,
  start_date = start_date,
  end_date = end_date,
  location = "large"
)

# write mini historical
mini_start_date <- "2019-10-10"
mini_end_date <- "2019-10-20"

write_weather_data(
  gridMET = gm_p,
  start_date = mini_start_date,
  end_date = mini_end_date,
  location = "mini"
)

# write mini historical
medium_start_date <- "2019-10-01"
medium_end_date <- "2020-10-01"

write_weather_data(
  gridMET = gm_p,
  start_date = medium_start_date,
  end_date = medium_end_date,
  location = "medium"
)



