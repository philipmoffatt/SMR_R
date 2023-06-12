# a set of functions for calculating lapse rates with gridMET data

library(tidyverse)
library(dplyr) 
library(lubridate)
library(AOI)
library(climateR)
library(mapview)
library(elevatr)

### File Description:
#     This function-set provides functions for downloading NOAA weather data,
#     finding the elevation of NOAA weather stations, and calculating lapse
#     rates between NOAA weather stations with two different methods.
# --------------------------------------------------------------------------- #


# Purpose: Downloads gridMET data using the climateR package based on the name
#          of a location.
# Parameters:
#   site_name: A character of the name of the site (e.g. "Pullman, WA" or "Moscow, ID")
# Returns:
#   A data frame containing the gridMET data: prcp, tmax, tmin, pet_grass, srad, wind_vel, vpd
# Notes:
#   There are no unit conversions except for converting temperature to Celsius
#   in this function so check param_meta$gridmet to see other units. 
pull_gridMET <- function(site_name) {
  
  # establish AOI
  AOI::aoi_get(site_name) %>% 
    AOI::aoi_map(returnMap = TRUE)
  
  # pull lat and lon for site
  cords <- geocode(site_name)
  lat <- cords$lat
  lon <- cords$lon
  print(lat)
  print(lon)
  
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
           l_turb = 2.5
    )
  
  return(gridMET_data)
  
}

# Purpose: Produces a linear model with a coefficient and intercept for the 
#          lapse rate of a variable with elevation.
# Parameters:
#   low_df: A data frame containing relevant weather data at the low elevation.
#   high_df: A data frame containing relevant weather data at the low elevation.
#   low_elev: A numeric elevation of the measurement site associated with low_df.
#   high_elev: A numeric elevation of the measurement site associated with high_df.
#   variable_name: A character class name of the variable for which the lapse rate 
#                  is being calculated. Used in plotting the lapse rate.
#   plot: A Boolean determining whether or not the function should produce a plot
#         of the lapse rate. Default is set to TRUE
# Returns:
#   A linear model object produced by the lm() function that contains a 
#   coefficient and intercept.
# Notes:
#   Double check units but they should be in the native gridMET units UNLESS
#   it is temperature (then they are in Celsius).
gridMET_model_lapse <- function(low_df, high_df, low_elev, high_elev, variable_name, plot = TRUE) {
  
  low_means <- mean(low_df[[variable_name]])
  high_means <- mean(high_df[[variable_name]])
  
  if (plot) {
    lm_data <- data.frame(x = c(low_elev, high_elev), y = c(low_means, high_means))
    
    plot(lm_data$x, lm_data$y,
         xlab = "Elevation", ylab = paste0("Mean ", variable_name),
         main = paste0("Mean ", variable_name , " vs. Elevation"))
    
    abline(lm(y ~ x, data = lm_data))
  }
  
  lm_data <- data.frame(x = c(low_elev, high_elev), y = c(low_means, high_means))
  
  lm_model <- lm(y ~ x, data = lm_data)
  
  return(lm_model)
}

# Purpose: Produces a lapse rate based on a simple 'difference' method.
# Parameters: 
#   low_df:  A data frame containing relevant weather data at the low elevation.
#   high_df: A data frame containing relevant weather data at the low elevation.
#   low_elev: A numeric elevation of the measurement site associated with low_df.
#   high_elev: A numeric elevation of the measurement site associated with high_df.
#   variable_name: A character class name of the variable for which the lapse rate 
#                  is being calculated. Used in plotting the lapse rate.
#   plot: A Boolean determining whether or not the function should produce a plot
#         of the lapse rate. Default is set to TRUE
# Returns:
#   A numeric class number representing the rate (in an assumed unit based on the 'units' parameter)
#   at which the variable being used changes with each meter of elevation. 
# Notes:
#   This is basically identical when using just two points. 
gridMET_difference_lapse <- function(low_df, high_df, low_elev, high_elev, variable_name, plot = TRUE) {
  
  low_means <- mean(low_df[[variable_name]])
  high_means <- mean(high_df[[variable_name]])
  
  if (plot) {
    lm_data <- data.frame(x = c(low_elev, high_elev), y = c(low_means, high_means))
    
    plot(lm_data$x, lm_data$y,
         xlab = "Elevation", ylab = paste0("Mean ", variable_name),
         main = paste0("Mean ", variable_name , " vs. Elevation"))
    
    abline(lm(y ~ x, data = lm_data))
    
  }
  
  lapse_rate <- (low_means - high_means) / (low_elev - high_elev)

  return(lapse_rate)
}



met_pullman <- pull_gridMET('Pullman, WA')
met_moscow <- pull_gridMET('Moscow Mountain, ID')

precip_lapse_model <- gridMET_model_lapse(met_pullman, met_moscow, variable_name = "prcp", 784, 1433)
pet_lapse_model <- gridMET_model_lapse(met_pullman, met_moscow, variable_name = "pet_grass", 784, 1433)
temp_lapse <- gridMET_difference_lapse(met_pullman, met_moscow, variable_name = "tavg", 784, 1433)



