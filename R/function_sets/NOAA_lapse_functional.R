library(tidyverse)
library(dplyr) 
library(lubridate)
library(AOI)
library(climateR)
library(mapview)
library(elevatr)
library(rnoaa)


# Purpose: Pull weather data from the NOAA database for a specific weather station
# Parameters:
#   site_id: Character string with the site ID of the weather station (e.g., "GHCND:USC00012345")
#   date_min: Optional parameter specifying the minimum date of the data to retrieve (default: NULL)
#   date_max: Optional parameter specifying the maximum date of the data to retrieve (default: NULL)
# Returns:
#   A data frame with the following columns:
#     date: Date of the observation
#     prcp: Precipitation, in tenths of mm
#     snow: Snowfall, in tenths of mm
#     snwd: Snow depth, in mm
#     tmax: Maximum temperature, in degrees Celsius
#     tmin: Minimum temperature, in degrees Celsius
#     tobs: Temperature at the time of observation, in degrees Celsius
# Notes:
#       If the optional date_min and date_max parameters are not provided, the function retrieves
#       all available data for the specified weather station. Otherwise, it retrieves data within
#       the specified date range. Missing values for specific variables are filled using the
#       preceding non-missing value.
pull_noaa <- function(site_id, date_min = NULL, date_max = NULL) {
  
  #site_id <- paste0("GHCND:", site_id)
  
  if (is.null(date_min) & is.null(date_max)) {
    wx <- meteo_tidy_ghcnd(site_id, var = 'all')
  } else {
    wx <- meteo_tidy_ghcnd(site_id, var = 'all', date_min = date_min, date_max = date_max)
  }

  to_fill <- c('prcp', 'snow', 'snwd', 'tmax', 'tmin', 'tobs')
  to_fill <- intersect(to_fill, colnames(wx))
  print(to_fill)
  
  wx <- wx %>%
    mutate(date = as.Date(date)) %>%
    complete(date = seq.Date(min(date), max(date), by="day")) %>%
    fill(to_fill)
  
  return(wx)
}

# Purpose: Finds the elevation of a NOAA station when given the station ID
# Parameters:
# Returns:
get_station_elevation <- function(site_id) {
  
  stations <- ghcnd_stations()
  stations_filtered <- stations[stations[, 'id'] == site_id, ]
  stations_elevations <- unique(stations_filtered$elevation)
  
  if (length(stations_elevations) == 1) {
    station_elevation <- stations_elevations[1]
    return(station_elevation)
  }
  else if (length(stations_elevations) > 1) {
    station_mean_elevation <- mean(stations_elevations, na.rm = TRUE)
    message("Multiple station elevations were found. Returning the mean of the elevations")
    return(station_mean_elevation)
  }
}
 
### NOAA linear model lapse rate
##  packages
#   parameters:
#     units: 'mm/10', 'mm', 'cm'
# returns
noaa_modeled_lapse <- function(low_df, high_df, low_elev, high_elev, variable_name, plot=TRUE, units="mm/10") {
  
  if (units == "mm/10") {
    low_means <- mean(low_df[[variable_name]], na.rm = TRUE)
    high_means <- mean(high_df[[variable_name]], na.rm = TRUE)
  }
  
  else if (units == "mm") {
    low_means <- mean(low_df[[variable_name]], na.rm = TRUE) / 10
    high_means <- mean(high_df[[variable_name]], na.rm = TRUE) / 10
  }
  
  else if (units == "cm") {
    low_means <- mean(low_df[[variable_name]], na.rm = TRUE) / 100
    high_means <- mean(high_df[[variable_name]], na.rm = TRUE) / 10
  }
  
  if (plot) {
    lm_data <- data.frame(x = c(low_elev, high_elev), y = c(low_means, high_means))
    
    plot(lm_data$x, lm_data$y,
         xlab = "Elevation", ylab = paste0("Mean ", variable_name),
         main = paste0("Mean ", variable_name, " (", units, ")", " vs. Elevation", " (meters)", " -- Modeled"))
    
    abline(lm(y ~ x, data = lm_data))
  }
  
  lm_data <- data.frame(x = c(low_elev, high_elev), y = c(low_means, high_means))
  
  lm_model <- lm(y ~ x, data = lm_data)
  
  return(lm_model)
}

## NOAA simple difference lapse rate
#  NOTE: this is identical to the modeled output so can probably be removed
noaa_difference_lapse <- function(low_df, high_df, low_elev, high_elev, variable_name, plot = TRUE, units="mm/10") {
  
  if (units == "mm/10") {
    low_means <- mean(low_df[[variable_name]], na.rm = TRUE)
    high_means <- mean(high_df[[variable_name]], na.rm = TRUE)
  }
  
  else if (units == "mm") {
    low_means <- mean(low_df[[variable_name]], na.rm = TRUE) / 10
    high_means <- mean(high_df[[variable_name]], na.rm = TRUE) / 10
  }
  
  else if (units == "cm") {
    low_means <- mean(low_df[[variable_name]], na.rm = TRUE) / 100
    high_means <- mean(high_df[[variable_name]], na.rm = TRUE) / 10
  }
  
  if (plot) {
    lm_data <- data.frame(x = c(low_elev, high_elev), y = c(low_means, high_means))
    
    plot(lm_data$x, lm_data$y,
         xlab = "Elevation", ylab = paste0("Mean ", variable_name),
         main = paste0("Mean ", variable_name, " (", units, ")", " vs. Elevation", " (meters)", " -- Differenced"))
    
    abline(lm(y ~ x, data = lm_data))
    
  }
  
  lapse_rate <- (low_means - high_means) / (low_elev - high_elev)
  
  return(lapse_rate)
}




## Example usage of these functions -- remove comments to test the code

### To find the station id for a NOAA site you are interested use this link:
#   https://www.ncdc.noaa.gov/cdo-web/search?datasetid=GHCND and take the value
#   that comes after 'GHCND:' in the box 'Network:ID'

#pullman_id <- 'USC00456789' # pullman station ID from the NOAA link above
#moscow_id <- 'USS0016C02S' # moscow moutain station ID from NOAA link above

#pullman_data <- pull_noaa(pullman_id) # pulling the pullman data using pullman_id
#moscow_data <- pull_noaa(moscow_id) # pulling the moscow mountain data using the moscow_id

#pullman_elevation <- get_station_elevation(pullman_id) # get elevation of pullman station using pullman_id
#moscow_elevation <- get_station_elevation(moscow_id) # get elevatino of moscow station using moscow_id

## calculating the lapse rate with the lm() function using the pulled data and elevations (units as 'cm')
#modeled_lapse_rate <- noaa_modeled_lapse( 
#  low_df = pullman_data,
#  high_df = moscow_data,
#  low_elev = pullman_elevation,
#  high_elev = moscow_elevation,
#  variable_name = "prcp",
#  plot = TRUE,
#  units = "cm"
#)

## calculating the lapse rate with a simple difference using the pulled data and elevations (units as 'cm')
#differenced_lapse_rate <- noaa_difference_lapse(
#  low_df = pullman_data,
#  high_df = moscow_data,
#  low_elev = pullman_elevation,
#  high_elev = moscow_elevation,
#  variable_name = "prcp",
#  plot = TRUE,
#  units = "cm"
#)




                         