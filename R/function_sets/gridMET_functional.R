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


## constant setup:
# initialize variables
lat <- 46.7312700
lon <- -117.1861
high_el <- 1433
low_el <- 784

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

## functions for building weather data
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

get_hourly_temps <- function(gridMET, latitude) {
  
  gridMET_hourly <- gridMET[,c("tmin", "tmax", "year", "month", "day")] %>%
    rename(Year = year, Day = day, Month = month, Tmin = tmin, Tmax = tmax) %>%
    make_hourly_temps(latitude, year_file = ., keep_sunrise_sunset = FALSE) %>%
    rename(year = Year, day = Day, month = Month, tmin = Tmin, tmax = Tmax)
  
  # Merge the two data.frames
  gridMET_hourly <- merge(gridMET, gridMET_hourly)
  
  return(gridMET_hourly)
}


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


lapse_rates <- function(low_site, high_site, low_site_el, high_site_el) {
  
  # calculate elevation difference in meters
  elevation_diff <- low_site_el - high_site_el
  
  # calculate precipitation lapse rate
  pullman_mean_prcp <- mean(low_site$prcp, na.rm = TRUE)
  moscow_mean_prcp <- mean(high_site$prcp, na.rm = TRUE)
  prcp_lapse_rate <- ((pullman_mean_prcp - moscow_mean_prcp) / elevation_diff) * 1000
  
  # calculate temperature lapse rate
  pullman_mean_tavg <- mean(low_site$tavg, na.rm = TRUE)
  moscow_mean_tavg <- mean(high_site$tavg, na.rm = TRUE)
  tavg_lapse_rate <- ((pullman_mean_tavg - moscow_mean_tavg) / elevation_diff) * 1000
  
  # calculate PET (grass) lapse rate
  pullman_mean_pet <- mean(low_site$pet_grass, na.rm = TRUE)
  moscow_mean_pet <- mean(high_site$pet_grass, na.rm = TRUE)
  pet_grass_lapse_rate <- ((pullman_mean_pet - moscow_mean_pet) / elevation_diff) * 1000
  
  # produce a linear model for PET based on site elevations and PET
  # setup variables for elevation and PET
  moscow_mountain_el <- high_site_el
  moscow_mountain_pet <- moscow_mean_pet
  
  pullman_el <- low_site_el
  pullman_pet <- pullman_mean_pet
  
  # create data frame with xy values
  pet_lm_data <- data.frame(x = c(pullman_pet,pullman_el), y = c(moscow_mountain_pet,moscow_mountain_el))
  
  # fit linear model
  model_values <- summary(lm(y ~ x, data = pet_lm_data))
  
  # return data frame with lapse rates and PET linear model values
  return(data.frame(precip_lapse_rate = prcp_lapse_rate,
                    temp_lapse_rate = tavg_lapse_rate,
                    pet_grass_lapse_rate = pet_grass_lapse_rate,
                    pet_lm_intercept = model_values$coefficients[1],
                    pet_lm_slope = model_values$coefficients[2]))
}

# calculate the dew point temperature
dewpoint_temperature <- function(gridMET) {
  
  # create a new column for dew point temperature
  gridMET$tdew <- case_when(gridMET$tmin <= 0 ~ gridMET$tmin * 0.7 - 2,
                            gridMET$tmin > 0 ~ gridMET$tmin * 0.97 - 0.53)
  
  return(gridMET)
}

# calculate the cloud fraction
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

## pressure_from_el() function: --> not really needed anymore
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

## estimated_wind_speed() function: 
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

## heat transfer to resistance function
#  params: 
#   wind_measurement_height: height of wind speed measurements (m)
#   temp_measurement_height: height of the air temperature measurements (m)
#   zero_plane_displacement_height: height of zero-plane displacement (m)
#   momentum_roughness: the momentum roughness parameter
#   heat_vapor_roughness: the heat and vapor roughness parameter
#   von_karman: von Karman's constant (0.41)
#   wind_speed: measured wind speed (m/s)
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

rename_columns <- function(gridMET, renamed_columns) {
  
  colnames(gridMET) <- renamed_columns
  
  return(gridMET)
  
}

subset_order_columns <- function(gridMET, desired_columns, renamed_columns) {
  gridMET_subset <- gridMET %>%
    dplyr::select(all_of(desired_columns)) %>%
    rename_columns(renamed_columns = renamed_columns)
  return(gridMET_subset)
}


# write weather data out to correct location -- probably needs some work
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



# calling functions
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
  calculate_cloud_fraction(lat, lon, -8) %>%
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
mini_end_date <- "2019-11-20"

write_weather_data(
  gridMET = gm_p,
  start_date = mini_start_date,
  end_date = mini_end_date,
  location = "mini"
)

# write mini historical
medium_start_date <- "2000-10-01"
medium_end_date <- "2001-10-01"

write_weather_data(
  gridMET = gm_p,
  start_date = medium_start_date,
  end_date = medium_end_date,
  location = "medium"
)



