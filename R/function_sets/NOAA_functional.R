# set of functions for downloading and processing NOAA data. this will be site
# specific but eventually should incorporate the 'stations within' feature. 

library(tidyverse)
library(tidyr)
library(dplyr) 
library(lubridate)
library(AOI)
library(climateR)
library(rnoaa)
library(chillR)
library(humidity)
library(bigleaf)


# pull data with rnoaa package
pull_noaa <- function(site_id, date_min = NULL, date_max = NULL) {
  if (is.null(date_min) & is.null(date_max)) {
    wx <- meteo_tidy_ghcnd(site_id, var = 'all')
  } else {
    wx <- meteo_tidy_ghcnd(site_id, var = 'all', date_min = date_min, date_max = date_max)
  }
  
  wx <- wx %>%
    mutate(date = as.Date(date)) %>%
    complete(date = seq.Date(min(date), max(date), by="day")) %>%
    fill(c('prcp', 'snow', 'snwd', 'tmax', 'tmin', 'tobs'))
  
  return(wx)
}
 
# pulling and cleaning sunshine data
pull_sunshine_data <- function(sunshine_file) {
  sunshine <- read.csv(sunshine_file) %>%
    mutate(date = as.Date(Date),
           doy = yday(date),
           year = lubridate::year(date),
           month = lubridate::month(date),
           day = lubridate::day(date),
           tsun_hour = tsun_min / 60) %>%
    complete(date = seq.Date(min(date), max(date), by="day")) %>%
    fill(c('tsun_hour', 'tsun_min', 'psun_percent')) %>%

  return(sunshine)
}

# format units
convert_units <- function(noaa_data) {
  noaa_data <- 
    noaa_data[, c("date", "tmax", "tmin", "tobs", "prcp", "snow")] %>%
    mutate(doy = yday(date),
           year = lubridate::year(date),
           month = lubridate::month(date),
           day = lubridate::day(date),
           tmin = tmin / 10,
           tmax = tmax / 10,
           tobs = tobs / 10,
           tavg = (tmax + tmin) / 2,
           prcp = prcp / 100)
  
  return(noaa_data)
}

# add l_turb variable as a constant
l_turb <- function(noaa_data, l_turb_value = 2.5) {
  noaa_data <- noaa_data %>%
    mutate(l_turb = l_turb_value)
  return(noaa_data)
}

# build out hourly temperatures for the 6 hour temperature matrix in SMR
get_hourly_temps <- function(noaa_data, latitude, keep_sunrise_sunset=FALSE) {
  
  noaa_subset <- noaa_data[, c("tmin", "tmax", "year", "month", "day", "date")] %>%
    rename(Year = year,
           Day = day,
           Month = month,
           Tmin = tmin,
           Tmax = tmax)
  
  noaa_subset_hourly <- make_hourly_temps(latitude = latitude, year_file = noaa_subset, keep_sunrise_sunset = keep_sunrise_sunset) %>%
    rename(year = Year,
           day = Day,
           month = Month,
           tmin = Tmin,
           tmax = Tmax)
  
  noaa_output <- merge(noaa_data, noaa_subset_hourly, by="date") %>%
    dplyr::select(-ends_with(".y")) %>% 
    rename_at(vars(ends_with(".x")), ~str_remove(., ".x")) %>% 
    mutate(
      year = lubridate::year(date), 
      month = lubridate::month(date), 
      day = lubridate::day(date))
  
  return(noaa_output)
} 


## Hamon's method for pet and helpers

avg_temp_by_month_year <- function(noaa_data) {
  avg_temp_by_month_year <- noaa_data %>%
    group_by(year, month) %>%
    summarise(avg_temp = mean(tavg, na.rm = TRUE))
  
  # Keep all original columns
  avg_temp_by_month_year <- left_join(noaa_data %>% dplyr::select(year, month), avg_temp_by_month_year, by = c("year", "month"))
  
  return(avg_temp_by_month_year)
}

all_month_year_combinations <- function(noaa_data) {
  all_years <- seq(min(noaa_data$year), max(noaa_data$year), by = 1)
  all_months <- seq(1, 12, by = 1)
  all_month_year <- expand.grid(year = all_years, month = all_months)
  
  # Add other original columns back into data frame
  all_month_year_combinations <- left_join(all_month_year, noaa_data %>% 
                                             distinct(year, month), 
                                           by = c("year", "month"))
  
  return(all_month_year_combinations)
}

daily_avg_temps <- function(noaa_data) {
  
  avg_temp_by_month_year <- avg_temp_by_month_year(noaa_data)
  all_month_year         <- all_month_year_combinations(noaa_data)
  
  daily_avg_temps        <- left_join(all_month_year,
                                      avg_temp_by_month_year,
                                      by=c("year","month")) %>% 
    mutate(day = '01',
           date=ymd(paste0(year,"-",month,"-",day))
    ) %>%
    complete(date = seq(min(date), max(date), by = "day")) %>%
    fill(c('year', 'month', 'avg_temp')) %>%
    mutate(day = lubridate::day(date)) %>%
    filter(date >= as.Date("1960-01-01"))
  
  # Add other original columns back into data frame
  daily_avg_temps <- left_join(daily_avg_temps, noaa_data, by = c("year", "month", "day"))
  
  # Remove columns with suffix ".y"
  daily_avg_temps <- dplyr::select(daily_avg_temps, !contains(".y"))
  
  # Remove suffix ".x" from column names
  colnames(daily_avg_temps) <- sub("\\.x$", "", colnames(daily_avg_temps))
  
  return(daily_avg_temps)
}

# to apply saturation vapor pressure function 
# entire column with specific parameters
apply_svp <- function(x) {
  return(SVP(x, isK = FALSE, formula = "Murray"))
}

# too apply the day_length function in chillR to an entire column 
apply_day_length <- function(x) {
  return(
    daylength(latitude =  46.7312700, JDay = x, notimes.as.na = FALSE )[[3]]
  )
}

add_historical_columns <- function(noaa_data) {
  
  # Generate daily average temperature data frame using helper functions
  daily_avg_temps_df <- avg_temp_by_month_year(noaa_data) %>%
    all_month_year_combinations()
  
  # Apply saturation vapor pressure function using apply_svp helper function
  noaa_data$sat_vapor_pressure <- apply(noaa_data[, "tavg", drop = FALSE], 1, apply_svp)
  
  # Calculate day length using latitude and chillR package with apply_day_length helper function
  noaa_data$day_length <- apply(noaa_data[, "doy", drop = FALSE], 1, apply_day_length)
  
  # Merge daily average temperatures data frame with NOAA data frame on "date" column
  merged_data <- merge(daily_avg_temps_df, noaa_data, by=c("year", "month"))
  
  return(merged_data)
}

# join noaa data with sunshine data
join_noaa_sunshine <- function(noaa_data, sunshine) {
  
  # Find the latest common start and earliest common ending dates
  latest_start <- max(min(noaa_data$date), min(sunshine$date))
  earliest_end <- min(max(noaa_data$date), max(sunshine$date))
  
  # Clip noaa_data and sunshine data so their dates align
  noaa_data_clipped <- subset(noaa_data, date >= latest_start & date <= earliest_end)
  sunshine_clipped <- subset(sunshine, date >= latest_start & date <= earliest_end)
  
  # Merge the two data frames on "date" column
  joined_data <- merge(noaa_data_clipped, sunshine_clipped, by = c("date", "year", "month", "day", "doy"))
  
  return(joined_data)
}

# equation for hamon pet which can produce pet in cm or mm
hamon_pet_equation <- function(proportionality_coefficient=1, day_length, sat_vapor_pressure, tavg, units='cm') {
  pet_hamon <- ((proportionality_coefficient * 0.165 * 216.7) * (day_length / 12)) *
    (sat_vapor_pressure / (tavg + 273.3))
  
  if (units == 'cm') {
    pet_hamon <- pet_hamon / 10
  }
  
  return(pet_hamon)
}

# wrapper for equation to apply it to a passed data frame
hamon_pet <- function(noaa_data, proportionality_coefficient=1, units='cm') {
  noaa_data$pet_hamon <- with(noaa_data, hamon_pet_equation(proportionality_coefficient, day_length, sat_vapor_pressure, tavg, units))
  return(noaa_data)
}

# calculate the cloud fraction from psun_percent from sunshine data
noaa_cloud_fraction <- function(noaa_data) {
  noaa_data$cloud <- 1 - (noaa_data$psun_percent / 100)
  return(noaa_data)
}

# build the crop coefficient curves
noaa_Kc_curve <- function(noaa_data, time_params, Kc_params, planting_offset, crop_name) {
  
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
  
  # attach the Kc curve to gridmet
  col_name <- paste0("Kc_", crop_name)
  noaa_data <- inner_join(noaa_data, Kc_frame, by = c("doy" = "doy"))
  colnames(noaa_data)[colnames(noaa_data) == "Kc"] <- col_name
  
  return(noaa_data)
}

# dew point temperature based on threshold method
noaa_dew_point <- function(noaa_data) {
  noaa_data$tdew <- case_when(noaa_data$tmin <= 0 ~ noaa_data$tmin * 0.7 - 2,
                              noaa_data$tmin > 0 ~ noaa_data$tmin * 0.97 - 0.53)
  return(noaa_data)
}

# read and process wind_data
read_and_process_wind_data <- function(file_path) {
  # Read and add date information
  wind_data <- read.csv(file_path)
  wind_data$date <- as.Date(wind_data$DATE, "%m/%d/%Y")
  
  # Rename avgWspd_kmph column to avg_wind_speed
  names(wind_data)[names(wind_data) == "avgWspd_kmph"] <- "avg_wind_speed"
  
  # Calculate mean and standard deviation of non-zero wind speeds
  non_zero_wind <- wind_data$avg_wind_speed[wind_data$avg_wind_speed != 0]
  mean_wind <- mean(non_zero_wind, na.rm = TRUE)
  sd_wind <- sd(non_zero_wind, na.rm = TRUE)
  
  # Generate values for zero wind speed entries
  zero_indices <- which(wind_data$avg_wind_speed <= 0)
  
  generated_values <- rnorm(length(zero_indices), mean = mean_wind, sd = sd_wind)
  
  # Re-sample if any generated value is below 0
  while(any(generated_values <=0)) {
    indices_to_resample <- which(generated_values <=0)
    generated_values[indices_to_resample] <- rnorm(length(indices_to_resample), mean = mean_wind, sd = sd_wind)
    
  }
  
  wind_data$avg_wind_speed[zero_indices] <- generated_values
  
  return(wind_data)
}

# join noaa and wind_data
merge_noaa_wind <- function(noaa_data, wind_data) {
  # Format wind_data
  wind_data$avg_wind_speed <- wind_data$avg_wind_speed / 3.6
  wind_data$year <- year(wind_data$date)
  
  wind_data_formatted <- wind_data %>%
    complete(date = seq.Date(min(date), max(date), by="day")) %>%
    fill(c('avg_wind_speed')) %>%
    mutate(year = year(date),
           month = month(date),
           day = day(date),
           doy = yday(date)) %>%
    dplyr::filter(year >= min(noaa_data[, "year"]) & year <= max(noaa_data[, "year"]))

  # Format noaa_data
  noaa_data_formatted <- noaa_data %>%
    mutate(
      date = ymd(paste(year, month, day))
    ) %>%
    complete(date = seq.Date(min(date), max(date), by="day"))
  
  # Merge formatted dataframes
  merged_wind_noaa <- merge(noaa_data_formatted, wind_data_formatted)
  
  return(merged_wind_noaa)
}



## aerodynamic resistance to heat transfer calculation
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
      
      print(displacement_height)
      print(momentum_roughness)
      print(heat_vapor_roughness)
      
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
    noaa_data,
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
      rh_values <- sapply(noaa_data$avg_wind_speed, function(wind_vel) {
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
      noaa_data[paste0("rh_", land_cover_name)] <- rh_values
      
    }
    
    return(noaa_data)
  }

rename_columns <- function(noaa_data, renamed_columns) {
  
  colnames(noaa_data) <- renamed_columns
  
  return(noaa_data)
  
}

subset_order_columns <- function(noaa_data, desired_columns, renamed_columns) {
  noaa_subset <- noaa_data %>%
    dplyr::select(all_of(desired_columns)) %>%
    rename_columns(renamed_columns = renamed_columns)
  return(noaa_subset)
}

# write weather data out to correct location -- probably needs some work
write_weather_data <- function(noaa_joined, start_date, end_date, location) {
  
  # Create the directory if it doesn't exist
  dir_path <- "./raw_data/weather"
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  # Select only rows within the specified date range
  noaa_joined_clipped <- noaa_joined[as.Date(noaa_joined$date) >= as.Date(start_date) &
                                             as.Date(noaa_joined$date) <= as.Date(end_date), ]
  
  # Write the clipped data to the file
  output_file <- file.path(dir_path, paste0("noaa_", location, ".csv"))
  write.table(noaa_joined_clipped,
              file = output_file,
              col.names = FALSE,
              row.names = FALSE,
              sep = " ")
  
}


# variable setup
ghcnd <- 'USC00456789'
date_min <- '1959-10-01'
date_max <- '1979-10-01'
latitude <- 46.7312700
sunshine_path <- "./raw_data/weather/sunshine_min_wallaWalla.csv"
wind_path <- "./raw_data/weather/Estimated Pullman Historic Wind Speed.csv"
start_date <- "1965-10-01"
end_date <- "1968-10-01"
  
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

# rh calculation variables
land_cover_names <- c('water', 'urban', 'forest', 'shrub', 'grass', 'row_crop')
land_cover_heights <- c(0.01, 4, 11, 0.3, 0.9, 0.9)
land_cover_names <- c('snow')
land_cover_heights <- c(0.01)

# column names
desired_columns <- c("date", "year", "month", "day", "doy", "tm.x", "tmin",
                     "tavg", "tdew", "prcp","pet_hamon","Hour_1","Hour_6","Hour_12","Hour_18",
                     "l_turb","cloud","Kc_water","Kc_urban","Kc_forest","Kc_shrub",
                     "Kc_grass", "Kc_row_crop","rh_snow")

smr_columns <- c("date", "year", "month", "day", "doy", "tmax", "tmin",
                 "tavg", "tdew", "precip","pet","hour_1","hour_6","hour_12","hour_18",
                 "l_turb","cloud","cc_water","cc_urban","cc_forest","cc_shrub",
                 "cc_grass","cc_row_crop","rh_snow")


# function calls
sunshine_data <- pull_sunshine_data(sunshine_path)
wind_data <- read_and_process_wind_data(wind_path)

wx <- pull_noaa(site_id = ghcnd, date_min = date_min, date_max = date_max) %>%
  convert_units %>%
  l_turb() %>%
  get_hourly_temps(latitude = latitude, keep_sunrise_sunset = FALSE) %>%
  add_historical_columns() %>%
  join_noaa_sunshine(sunshine = sunshine_data) %>%
  hamon_pet() %>% 
  noaa_cloud_fraction() %>%
  noaa_Kc_curve(water_time, water_Kc, 0, 'water') %>%
  noaa_Kc_curve(urban_time, urban_Kc, 0, 'urban') %>%
  noaa_Kc_curve(forest_time, forest_Kc, 0, 'forest') %>%
  noaa_Kc_curve(shrub_time, shrub_Kc, 0, 'shrub') %>%
  noaa_Kc_curve(grass_time, grass_Kc, 0, 'grass') %>%
  noaa_Kc_curve(row_crop_time, row_crop_Kc, 274, 'row_crop') %>%
  noaa_dew_point() %>%
  merge_noaa_wind(wind_data = wind_data) %>%
  land_cover_rh_wrapper(land_cover_names = land_cover_names,
                        land_cover_heights = land_cover_heights,
                        wind_measurement_height = 2,
                        temp_measurement_height = 2,
                        von_karman = 0.41) %>% 
  subset_order_columns(desired_columns = desired_columns,
                       renamed_columns = smr_columns)

## ugly calibration fix -- will work on it tomorrow
wx$pet <- wx$pet * 1.69
d <- wx


d %>% 
  group_by(year) %>%
  summarise(sum_pet = sum(pet), tot_p = sum(precip)) %>%
  ggplot()+
  geom_point(aes(year,sum_pet))+
  geom_line(aes(year,tot_p))


# write big historical
write_weather_data(
  noaa_joined = wx,
  start_date = start_date,
  end_date = end_date,
  location = "pullman"
)
  
# write mini historical
mini_start_date <- "1965-10-01"
mini_end_date <- "1965-10-10"

write_weather_data(
  noaa_joined = wx,
  start_date = mini_start_date,
  end_date = mini_end_date,
  location = "pullman_mini"
)

# write mini historical
medium_start_date <- "1965-10-01"
medium_end_date <- "1966-10-01"

write_weather_data(
  noaa_joined = wx,
  start_date = medium_start_date,
  end_date = medium_end_date,
  location = "pullman_medium"
)

