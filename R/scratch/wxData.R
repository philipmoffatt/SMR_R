# Weather data for SMR model - pl script format

# wx data needed for mica creek script - doy	year	tmax	tmin	tavg	precip	pet_snotel	
# for ET coefficient calculation - cc_forest	cc_partial	cc_open	
# for controlling output maps - output	
# 6-hr temp data disaggregated using tempdisagg - t0	t1	t2	t3	t4

library(tidyverse)
library(dplyr) 
library(lubridate)
library(AOI)
library(climateR)
library(rnoaa)
library(chillR)
library(humidity)
library(bigleaf)

### ----- PULLING DATA FROM NOAA AND GRIDMET ----- ###

# establish AOI
site <- 'Pullman'

AOI::aoi_get(site) %>% 
  AOI::aoi_map(returnMap = TRUE)

# use noaa
# https://cran.r-project.org/web/packages/rnoaa/rnoaa.pdf

# pull lat and lon for site
cords <- geocode(site)
lat <- cords$lat
lon <- cords$lon

# pulling stations within 10 km of the site ####
# * not tested
# station_data <- ghcnd_stations() # Takes a while to run
# lat_lon_df <- data.frame(id = site,
#                          latitude = lat,
#                          longitude = lon)
# nearby_stations <- meteo_nearby_stations(lat_lon_df = lat_lon_df,
#                                          station_data = station_data, radius = 10)


# pull data from PCFS site ####
# PULLMAN 2 NW, WA US / 	GHCND:USC00456789 / 	46.76016°, -117.1861°
# Elevation	766.6 m / 1940-10-21 - 	2023-02-17
wx <- meteo_tidy_ghcnd('USC00456789', var = 'all', date_min = '1959-10-01', date_max = '1979-10-01') %>%
  mutate(date = as.Date(date)) %>%
  complete(date = seq.Date(min(date), max(date), by="day")) %>%
  fill(c('prcp', 'snow', 'snwd', 'tmax', 'tmin', 'tobs'))


# use gridMET ####
# https://github.com/mikejohnson51/climateR

gridMETdata <- param_meta$gridmet
# pull precipitation (mm), temp (C), pet (mm) from GRIDMET
AOI_p = AOI::geocode('Pullman', pt = TRUE)
gridMET_p  = as.data.frame(getGridMET(AOI_p, param = c('prcp', 'tmax', 'tmin', 'pet_grass', 'srad', 'wind_vel'), 
                                         startDate = "2000-01-01", 
                                         endDate = "2022-12-31"))

# select date, prcp, tmax, tmin and add tavg, year, month, day
pullman <- gridMET_p[,c(4:8)] %>% 
  mutate(year = lubridate::year(date), 
  month = lubridate::month(date), 
  day = lubridate::day(date),
  tmax = tmax - 273.15,
  tmin = tmin - 273.15,
  tavg = (tmax+tmin)/2)

# gridMET pull for moscow mountain 
AOI_m = AOI::geocode('Moscow Mountain', pt = TRUE)
gridMET_m  = as.data.frame(getGridMET(AOI_m, param = c('prcp', 'tmax', 'tmin', 'pet_grass', 'srad', 'wind_vel'), 
                                      startDate = "2000-01-01", 
                                      endDate = "2022-12-31"))

### --------------------------------------- ###


### ----- INITIAL FORMATTING TO FIT PERL ----- ###

# add in day-of-year and separate dates to get year on its own
gridMET_p <- gridMET_p %>%
  mutate(doy = yday(date),
         year = lubridate::year(date),
         month = lubridate::month(date),
         day = lubridate::day(date),
         tmin = tmin - 273.15,
         tmax = tmax - 273.15,
         tavg = (tmax+tmin)/2,
         output = ifelse(month == 2, 1, 0),
         l_turb = 2.5,
         prcp = prcp/10,
         pet_grass = pet_grass/10
         )

gridMET_m <- gridMET_m %>%
  mutate(doy = yday(date),
         year = lubridate::year(date),
         month = lubridate::month(date),
         day = lubridate::day(date),
         tmin = tmin - 273.15,
         tmax = tmax - 273.15,
         tavg = (tmax+tmin)/2,
         output = ifelse(month == 2, 1, 0),
         l_turb = 2.5,
         prcp = prcp/10,
         pet_grass = pet_grass/10
         )

### --------------------------------------- ###


### ----- DISAGGREGATING INTO HOURLY TEMPERATURES ----- ###

# leaving the original gridMET_p frame copying new one
data_p <- gridMET_p[,c("tmin", "tmax", "year", "month", "day")] %>% 
  rename(Year = year, 
         Day = day, 
         Month = month,
         Tmin = tmin, 
         Tmax = tmax
         )

# extracting hourly temperatures from the minimum and maximum temperatures
# pulling four evenly spaced hours (1, 6, 12, 18) to make the 'temperature matrix'
# joining the original gridMET_p data with the new temperature matrix
gridMET_p_hourly <- make_hourly_temps(latitude = 46.7312700, year_file = data_p, keep_sunrise_sunset = FALSE) %>%
  rename(year = Year, 
         day = Day, 
         month = Month,
         tmin = Tmin, 
         tmax = Tmax)

  gridMET_p_hourly <- gridMET_p_hourly[, c("year", "month", "day", "Hour_1", "Hour_6", "Hour_12", "Hour_18")]

# joining original gridMET_p with new hourly data to add in temperature matrix 
gridMET_joined_p <- inner_join(gridMET_p, gridMET_p_hourly)

# do the same hourly calculation for moscow
data_m <- gridMET_m[,c("tmin", "tmax", "year", "month", "day")] %>% 
  rename(Year = year, 
         Day = day, 
         Month = month,
         Tmin = tmin, 
         Tmax = tmax
  )

# extracting hourly temperatures from the minimum and maximum temperatures
# pulling four evenly spaced hours (1, 6, 12, 18) to make the 'temperature matrix'
# joining the original gridMET_p data with the new temperature matrix
gridMET_m_hourly <- make_hourly_temps(latitude = 46.7312700, year_file = data_m, keep_sunrise_sunset = FALSE) %>%
  rename(year = Year, 
         day = Day, 
         month = Month,
         tmin = Tmin, 
         tmax = Tmax)

# pull date information for merging and pulling 4 6-hour-spaced temperatures
gridMET_m_hourly <- gridMET_m_hourly[, c("year", "month", "day", "Hour_1", "Hour_6", "Hour_12", "Hour_18")] 

# join all gridMET data with the new hourly matrix based on date information
gridMET_joined_m <- inner_join(gridMET_m, gridMET_m_hourly)

### ----- ADDING IN THE CROP COEFFICIENTS----- ###
# https://www.fao.org/3/x0490e/x0490e0b.htm
  # 1. water: 0.6525 (always)
  # 2. urban/rock/barren: 0.30 (always) --> based on weather file and also within reason from FAO early bare soil coefficients 
  # 3. forest/woody wetlands: 1 (conifer, always), 1.05-1.10 (for short veg wetlands) --> going to start with conifer because its simple
  # 4. shrub: (time: deciduous orchard high latitudes, Kc: ) --> these should probably be improved upon
  # 5. grass/grassy wetlands/pasture: (FAO: grazing pasture) --> using "extensive" instead of "rotated" --> based on frost dates
  # 6. row crop: using wheat for now (FAO 35-45 latitude)

# function for building Kc curve given "growth stage length" and "Kc" dictionaries
ann_Kc_curve <- function(time_params, Kc_params) {
  
  preplant <- rep(Kc_params[["Kc_preplant"]], time_params[["preplant_len"]])
  print(length(preplant))
  
  init_section <- rep(Kc_params[["Kc_init"]], time_params[["init_len"]])
  print(length(init_section))
  
  dev_growth_rate <- (Kc_params[["Kc_mid"]] - Kc_params[["Kc_init"]]) / time_params[["dev_len"]]
  dev_arr <- seq(Kc_params[["Kc_init"]], Kc_params[["Kc_mid"]], length.out = time_params[["dev_len"]])
  dev_section <- round(dev_arr, 4)
  print(length(dev_section))
  
  mid_section <- rep(Kc_params[["Kc_mid"]], time_params[["mid_len"]])
  print(length(mid_section))
  
  late_growth_rate <- (Kc_params[["Kc_end"]] - Kc_params[["Kc_mid"]]) / time_params[["late_len"]]
  late_arr <- seq(Kc_params[["Kc_mid"]], Kc_params[["Kc_end"]], length.out = time_params[["late_len"]])
  late_section <- round(late_arr, 4)
  print(length(late_section))
  
  Kc_post_harv <- rep(Kc_params[["Kc_post_harv"]], time_params[["post_len"]])
  print(length(Kc_post_harv))
  
  total_curve = c(preplant, init_section, dev_section, mid_section, late_section, Kc_post_harv)
  
  # checking that all sections add up to a year
  print(length(total_curve))
  
  return(total_curve)
  
}

# setting up dictionaries for surface with static curves
water_time <- list('preplant_len'=60, 'init_len'=60, 'dev_len'=60, 'mid_len'=60, 'late_len'=60, 'post_len'=65)
water_Kc <- list('Kc_preplant'=0.6525, 'Kc_init'=0.6525, 'Kc_mid'=0.6525, 'Kc_end'=0.6525, 'Kc_post_harv'=0.6525)

urban_time <- list('preplant_len'=c(60), 'init_len'=c(60), 'dev_len'=c(60), 'mid_len'=c(60), 'late_len'=c(60), 'post_len'=c(65))
urban_Kc <- list('Kc_preplant'=0.30, 'Kc_init'=0.30, 'Kc_mid'=0.30, 'Kc_end'=0.30, 'Kc_post_harv'=0.30)

# using FAO conifer values for now meaning it's static but that can change
forest_time <- list('preplant_len'=60, 'init_len'=60, 'dev_len'=60, 'mid_len'=60, 'late_len'=60, 'post_len'=65)
forest_Kc <- list('Kc_preplant'=1, 'Kc_init'=1, 'Kc_mid'=1, 'Kc_end'=1, 'Kc_post_harv'=1)

# using deciduous orchard
shrub_time <- list('preplant_len'=60, 'init_len'=20, 'dev_len'=70, 'mid_len'=90, 'late_len'=30, 'post_len'=95)
shrub_Kc <- list('Kc_preplant'=0.45, 'Kc_init'=0.45, 'Kc_mid'=0.95, 'Kc_end'=0.70, 'Kc_post_harv'=0.45)

#relied on frost dates from AgWeatherNet: Pullman station: frost on --> May 1st, frost off --> September 19th --> splitting resiudal (125) in half for now
grass_time <- list('preplant_len'=113, 'init_len'=10, 'dev_len'=20, 'mid_len'=64, 'late_len'=62, 'post_len'=96)
grass_Kc <- list('Kc_preplant'=0.30, 'Kc_init'=0.30, 'Kc_mid'=0.75, 'Kc_end'=0.75, 'Kc_post_harv'=0.30)

# using spring wheat values for now
#row_crop_time <- list('preplant_len'=75, 'init_len'=20, 'dev_len'=25, 'mid_len'=60, 'late_len'=30, 'post_len'=155)
row_crop_time <- list('preplant_len'=0, 'init_len'=160, 'dev_len'=75, 'mid_len'=75, 'late_len'=25, 'post_len'=30)
row_crop_Kc <- list('Kc_preplant'=0.30, 'Kc_init'=0.7, 'Kc_mid'=1.15, 'Kc_end'=0.35, 'Kc_post_harv'=0.30)


# apply ann_Kc_curve to build annual curves for each land cover type
water_curve <- data.frame(cc_water = ann_Kc_curve(time_params = water_time, Kc_params = water_Kc), doy = as.numeric(1:365))
urban_curve <- data.frame(cc_urban = ann_Kc_curve(time_params = urban_time, Kc_params = urban_Kc), doy = as.numeric(1:365))
forest_curve <- data.frame(cc_forest = ann_Kc_curve(time_params = forest_time, Kc_params = forest_Kc), doy = as.numeric(1:365))
shrub_curve <- data.frame(cc_shrub = ann_Kc_curve(time_params = shrub_time, Kc_params = shrub_Kc), doy = as.numeric(1:365))
grass_curve <- data.frame(cc_grass = ann_Kc_curve(time_params = grass_time, Kc_params = grass_Kc), doy = as.numeric(1:365))
row_crop_curve <- data.frame(cc_row_crop = ann_Kc_curve(time_params = row_crop_time, Kc_params = row_crop_Kc), doy = as.numeric(1:365))
row_crop_curve$doy <- (row_crop_curve$doy + 274) %% 365 


# JOINING THE KC CURVES WITH THE HISTORICAL, PULLMAN, AND MOSCOW DATASETS

# add the Kc curves to Pullman
gridMET_joined_p = inner_join(gridMET_joined_p, water_curve)
gridMET_joined_p = inner_join(gridMET_joined_p, urban_curve)
gridMET_joined_p = inner_join(gridMET_joined_p, forest_curve)
gridMET_joined_p = inner_join(gridMET_joined_p, shrub_curve)
gridMET_joined_p = inner_join(gridMET_joined_p, grass_curve)
gridMET_joined_p = inner_join(gridMET_joined_p, row_crop_curve)

# add the Kc curves to Moscow
gridMET_joined_m = inner_join(gridMET_joined_m, water_curve)
gridMET_joined_m = inner_join(gridMET_joined_m, urban_curve)
gridMET_joined_m = inner_join(gridMET_joined_m, forest_curve)
gridMET_joined_m = inner_join(gridMET_joined_m, shrub_curve)
gridMET_joined_m = inner_join(gridMET_joined_m, grass_curve)
gridMET_joined_m = inner_join(gridMET_joined_m, row_crop_curve)


# FORMATTING COLUMN ORDER TO MATCH THE PERL SCRIPT

# make a vector of column names for sub-setting and reordering the data frames
gridMET_name_vector <- c(
  'year', 'month', 'day', 'doy', 'tmax', 'tmin', 'tavg', 'prcp', 'Hour_1', 'Hour_6', 'Hour_12', 'Hour_18', 'pet_grass', 'cc_water', 'cc_urban', 'cc_forest', 'cc_shrub', 'cc_grass', 'cc_row_crop', 'l_turb', 'output', 'srad', 'wind_vel'
)

# sub-setting and reording gridMET files
gridMET_joined_p <- select(gridMET_joined_p, all_of(gridMET_name_vector))
gridMET_joined_m <- select(gridMET_joined_m, all_of(gridMET_name_vector))

### ----- CALCULATING EVAPOTRANSPIRATION LAPSE RATES ----- ###
# calculate elevation difference in meters
# pullman station elevation: https://wsdot.com/travel/real-time/weather/4024
# Moscow mountain site elevation: comes with data

elevation_diff <- 1433 - 784

pullman_mean_pet <- mean(gridMET_joined_p$pet_grass, na.rm = TRUE)
moscow_mean_pet <- mean(gridMET_joined_m$pet_grass, na.rm = TRUE)

# pet / km
pet_grass_lapse_rate <- ((pullman_mean_pet - moscow_mean_pet) / elevation_diff) * 1000

# produce a linear model for PET based on site elevations and PET

# setup variables for elevation and PET
moscow_mountain_el <- 1433
moscow_mountain_pet <- mean(gridMET_joined_m[, "pet_grass"], na.rm=TRUE)
  
pullman_el <- 784
pullman_pet <- mean(gridMET_joined_p[, "pet_grass"], na.rm=TRUE)

# create dataframe with xy values
pet_lm_data <- data.frame(x = c(pullman_pet,pullman_el), y = c(moscow_mountain_pet,moscow_mountain_el))

# fit linear model
model_values <- summary(lm(y ~ x, data = pet_lm_data))

### ------------------------------------------------------- ###


### ----- DEW POINT TEMPERATURE ----- ###
# dew point temperature is calculated using (Gentilli, J. 1955) method based 
# and is based on minimum temperature

# threshold based calculation: gridMET_joined_p
gridMET_joined_p$tdew <- case_when(gridMET_joined_p$tmin <= 0 ~ gridMET_joined_p$tmin * 0.7 - 2,
                                    gridMET_joined_p$tmin > 0 ~ gridMET_joined_p$tmin * 0.97 - 0.53)

# threshold based calculation: gridMET_joined_m
gridMET_joined_m$tdew <- case_when(gridMET_joined_m$tmin <= 0 ~ gridMET_joined_m$tmin * 0.7 - 2,
                                    gridMET_joined_m$tmin > 0 ~ gridMET_joined_m$tmin * 0.97 - 0.53)

### ------------------------------------------------------- ###


### ----- CLOUD FRACTION ----- ###
# cloud fraction is calculated by the ratio of the observed solar radiation to 
# potential solar radiation at that point: using the bigleaf ecosystem 
# properties library

# setup sequence of days for potential radiation function
julian_year <- seq(1, 365)

# calculate the sunrise, sunset, and daytime hours for each day based on day of year
# and latitude
annual_day_lengths_p <- daylength(latitude = 46.76016, JDay = julian_year, notimes.as.na = FALSE)
annual_day_lengths_m <- daylength(latitude = 46.8035, JDay = julian_year, notimes.as.na = FALSE)

# setup solar radiation data frame for merging with gridMET frames
srad_df <- 
  data.frame(
    doy = numeric(), 
    srad_potential_p = numeric(), 
    srad_potential_m = numeric()
    )

for (day in julian_year) {
  
  # get the sunrise and sunset time for the corresponding day of year and build
  # a sequence from them --> the 0.1 step is arbitrary
  day_sequence_p <- 
    seq(
      from = annual_day_lengths_p$Sunrise[day], 
      to = annual_day_lengths_p$Sunset[day],
      by = 0.1
        )
  
  day_sequence_m <- 
    seq(
      from = annual_day_lengths_m$Sunrise[day], 
      to = annual_day_lengths_m$Sunset[day],
      by = 0.1
    )
  
  # calculate the potential radiation for the sequence
  potential_radiation_p <- 
    potential.radiation(
      doy = day,
      hour = day_sequence_p,
      latDeg = 46.76016,
      longDeg = -117.1861,
      timezone = -8,
      useSolartime = TRUE
    )
  
  potential_radiation_m <- 
    potential.radiation(
      doy = day,
      hour = day_sequence_m,
      latDeg = 46.8035,
      longDeg = -116.8688,
      timezone = -8,
      useSolartime = TRUE
    )
  
  # take the average potential radiation to match a daily summary
  mean_daily_srad_p <- mean(potential_radiation_p)
  mean_daily_srad_m <- mean(potential_radiation_m)
  
  
  # add average to corresponding position in the potential_srad dataframe
  srad_df[day, "doy"] <- day
  srad_df[day, "srad_potential_p"] <- mean_daily_srad_p
  srad_df[day, "srad_potential_m"] <- mean_daily_srad_m
  
  
}

# add the repeating potential solar radiation frame to the gridMET pullman frame
# based on doy
gridMET_joined_p <- inner_join(gridMET_joined_p, srad_df[, c("doy", "srad_potential_p")], by = "doy")
gridMET_joined_m <- inner_join(gridMET_joined_m, srad_df[, c("doy", "srad_potential_m")], by = "doy")

# cloud fraction for low site and high site is the:
# solar radiation / potential solar radiation
gridMET_joined_p$cloud <- 
  gridMET_joined_p$srad / gridMET_joined_p$srad_potential_p

gridMET_joined_m$cloud <- 
  gridMET_joined_m$srad / gridMET_joined_m$srad_potential_m

### ------------------------------------------------------- ###

### ----- RESISTANCE TO HEAT TRANSFER ----- ###
# NOTE!!: there is a problem with this equation in that with taller vegetation
#         you quickly get ln(-number) and it isn't useful depending on the 
#         height of the wind measurement
# resistance to heat transfer is dependent on land cover
# and is calculated with equations from Walter et al. (2004)
# more detailed info on equation: "Environmental Biophysics" Campbell (1977)

# variables needed for resistance to heat transfer (rh): 
#   1. zu: height of the wind speed measurement (constant, meters)
#   2. zt: height of the air temperature measurement (constant, meters)
#   3. d: height of the zero plane displacement 
#         (constant and a function of height)
#   4. zm: the momentum roughness parameter (calculated as a function 
#          of height)
#   5. zh: the heat and vapor roughness parameter (calculated as a function
#          of height)
#   6. k: Von Karman's constant (0.41)
#   7. u: the wind speed (dynamic part of the function, meters/second)

# setup constants that won't change between land covers (should confirm these)
zu_value <- 2
zt_value <- 1
k_value <- 0.41

# making a vector of the land covers heights (meters)
# order of covers: water, urban, forest, shrub, grass, row crop
land_cover_names <- c('water', 'urban', 'forest', 'shrub', 'grass', 'row_crop')
land_cover_heights <- c(0.01, 0.01, 1.3, 1.3, 1.3, 1.3)

heat_transfer_resistance <- function(zu, zt, k, d, zm, zh, u) {
  
  # maybe a better way to write this
  rh <- ( ( log((zu - d + zm) / (zm)) * log((zt - d + zh) / (zh)) ) / ( k^2 * u) ) #* (1 / 86400)
  
  return(rh)
  
}

# define zm outside of the for loop
zm_value <- 0.1 * land_cover_heights

for (landcover in 1:length(land_cover_heights)) {
  
  height <- land_cover_heights[landcover]
  
  # zero plane displacement height
  d_value <- 0.65 * height
  
  # use zm_value instead of zm in the function call
  rh <- heat_transfer_resistance(
    zu = zu_value,
    zt = zt_value,
    k = k_value,
    d = d_value,
    zm = zm_value[landcover],
    zh = 0.2 * zm_value[landcover],
    u = gridMET_joined_p$wind_vel
  )
  
  gridMET_joined_p[, paste0("rh_", land_cover_names[landcover])] <- rh
  
}

### ------------------------------------------------------- ###


### ----- INTERMEDIATE FORMATTING FOR WEATHER FILES ----- ###

gridMET_name_vector <- c(
  'date', 'year', 'month', 'day', 'doy', 'tmax', 'tmin', 'tavg', 'tdew', 'prcp', 'pet_grass', 'Hour_1', 'Hour_6', 'Hour_12', 'Hour_18', 'l_turb', 'cloud', 'cc_water', 'cc_urban', 'cc_forest', 'cc_shrub', 'cc_grass', 'cc_row_crop', 'rh_water', 'rh_urban', 'rh_forest', 'rh_shrub', 'rh_grass', 'rh_row_crop', 'output' 
)

gridMET_joined_p <- gridMET_joined_p %>%
  mutate(date = make_date(year, month, day))

gridMET_joined_p <- gridMET_joined_p[, gridMET_name_vector] %>%
  rename(pet = pet_grass, hour_1 = Hour_1, hour_6 = Hour_6, hour_12 = Hour_12, hour_18 = Hour_18)

### ------------------------------------------------------- ###


##### ---------- ALL HISTORICAL DATA PROCESSING --------- ######

### ----- SETUP FOR HISTORICAL WEATHER FILE ----- ###

# acting on new variable: historical_weather instead of wx from now on and adding constants:
# l_turb, output, doy, and splitting the date up
historical_weather <- wx[, c("date", "tmax", "tmin", "tobs", "prcp", "snow")] %>%
  mutate(doy = yday(date),
         year = lubridate::year(date),
         month = lubridate::month(date),
         day = lubridate::day(date),
         tmin = tmin/10,
         tmax = tmax/10,
         tobs = tobs/10,
         tavg = (tmax+tmin)/2,
         output = ifelse(month == 2, 1, 0),
         l_turb = 2.5)

historical_weather$prcp = historical_weather$prcp / 10

data_historical <- historical_weather[,c("tmin", "tmax", "year", "month", "day", "date")] %>% 
  rename(Year = year, 
         Day = day, 
         Month = month,
         Tmin = tmin, 
         Tmax = tmax
  )

historical_hourly <- make_hourly_temps(latitude = 46.7312700, year_file = data_historical, keep_sunrise_sunset = FALSE) %>%
  rename(year = Year, 
         day = Day, 
         month = Month,
         tmin = Tmin, 
         tmax = Tmax)

historical_hourly <- historical_hourly[, c("date", "year", "month", "day", "Hour_1", "Hour_6", "Hour_12", "Hour_18")]

historical_joined <- merge(historical_weather, historical_hourly, by="date") %>%
  #select(-day.y, -month.y, -year.y, -day.x, 
  #       -month.x, -year.x, -tmin.y, -tmax.y) %>%
  mutate(
    year = lubridate::year(date), 
    month = lubridate::month(date), 
    day = lubridate::day(date)) #%>%
#rename(
#  tmin = tmin.x,
#  tmax = tmax.x)

#historical_joined <- historical_joined[!apply(historical_joined, 1, function(row) any(is.na(row))), ]

### ----- HAMON METHOD FOR POTENTIAL EVAPOTRANSPIRATION ----- #####
# link to method: http://data.snap.uaf.edu/data/Base/AK_2km/PET/Hamon_PET_equations.pdf
# NOTE: this is all using Pullman data
#PET is potential evapotranspiration [mm day-1]
#k is proportionality coefficient = 1^1 [unitless]
#N is daytime length [x/12 hours]
#es is saturation vapor pressure [mb]
#T is average monthly temperature [°C]

# produce average monthly temperature columns and adding them as daily values to historical_joined
avg_temp_by_month_year <- historical_joined %>%
  group_by(year, month) %>%
  summarise(avg_temp = mean(tavg))

# create a new data frame with all possible combinations of years and months
all_years <- seq(min(historical_joined$year), max(historical_joined$year), by = 1)
all_months <- seq(1, 12, by = 1)
all_month_year <- expand.grid(year = all_years, month = all_months)

# left-join average temperature values from the original dataframe
daily_avg_temps <- left_join(all_month_year, avg_temp_by_month_year, by = c("year", "month")) %>%
  mutate(day = 1,
         date = ymd(paste(year, month, day, sep = "-"))
  ) %>%
  complete(date = seq(min(date), max(date), by = "day")) %>%
  fill(c('year', 'month', 'avg_temp')) %>%
  mutate(day = lubridate::day(date)) %>%
  filter(date >= as.Date("1960-01-01"))

# read in sunshine data
sunshine <- read.csv("./raw_data/weather/sunshine_min_wallaWalla.csv") %>%
  mutate(date = as.Date(Date),
         doy = yday(date),
         year = lubridate::year(date),
         month = lubridate::month(date),
         day = lubridate::day(date),
         tsun_hour = tsun_min / 60) %>%
  complete(date = seq.Date(min(date), max(date), by="day")) %>%
  fill(c('tsun_hour', 'tsun_min', 'psun_percent'))


# to apply SVP() function entire column with specific parameters
apply_svp <- function(x) {
  return(SVP(x, isK = FALSE, formula = "Clausius-Clapeyron"))
}

# too apply the day_length function in chillR to an entire column 
apply_day_length <- function(x) {
  return(
    daylength(latitude =  46.7312700, JDay = x, notimes.as.na = FALSE )[[3]]
  )
}

# HAMON METHOD VARIABLES
# add in monthly averages as a daily record
historical_joined <- merge(historical_joined, daily_avg_temps)
# calculate new saturation vapor pressure columnusing chillR package
historical_joined$sat_vapor_pressure <- apply(historical_joined[, "tavg", drop = FALSE], 1, apply_svp)
# calculate day_length using latitude and chillR package
historical_joined$day_length <- apply(historical_joined[, "doy", drop = FALSE], 1, apply_day_length)
# setup necessary variables for Hamon's method
proportionality_coefficient <- 1

# JOINING THE SUNSHINE AND HISTORICAL FILE SO DATES MATCH UP
# find the latest common start and the earliest common ending dates
latest_start <- max(min(historical_joined$date), min(sunshine$date))
earliest_end <- min(max(historical_joined$date), max(sunshine$date))

# clip historical and sunshine data so their dates align 
historical_joined <- subset(historical_joined, date >= latest_start & date <= earliest_end)
sunshine <- subset(sunshine, date >= latest_start & date <= earliest_end)

# merge the two frames
historical_joined <- merge(historical_joined, sunshine)

# HAMON PET AND CLOUD COVER FRACTION CALCULATION
# calculate pet_hamon (with Hamon formula) --> pet is cm not mm
historical_joined$pet_hamon <- proportionality_coefficient * 0.165 * 216.7 * (historical_joined$day_length / 12) * (historical_joined$sat_vapor_pressure / (historical_joined$avg_temp + 273.3)) / 10 # /10 to get cm

# calculate cloud fraction as: 1 - (percent_sunshine / 100) --> 1 = full cloud cover, 0 = no cloud cover --> for all dataframes
historical_joined$cloud <- 1 - (historical_joined$psun_percent / 100)

# add the Kc curves to historical_joined
historical_joined = inner_join(historical_joined, water_curve)
historical_joined = inner_join(historical_joined, urban_curve)
historical_joined = inner_join(historical_joined, forest_curve)
historical_joined = inner_join(historical_joined, shrub_curve)
historical_joined = inner_join(historical_joined, grass_curve)
historical_joined = inner_join(historical_joined, row_crop_curve)

# intermediate renaming for historical
historical_name_vector <- c(
  'year', 'month', 'day', 'doy', 'tmax', 'tmin', 'tavg', 'prcp',  'Hour_1', 'Hour_6', 'Hour_12', 'Hour_18', 'pet_hamon', 'cc_water', 'cc_urban', 'cc_forest', 'cc_shrub', 'cc_grass', 'cc_row_crop', 'cloud', 'l_turb', 'output' 
)

# sub-setting and reording historical file
historical_joined <- select(historical_joined, all_of(historical_name_vector))

# dew point historical calculation
# threshold based calculation: historical_joined
historical_joined$tdew <- case_when(historical_joined$tmin <= 0 ~ historical_joined$tmin * 0.7 - 2,
                                    historical_joined$tmin > 0 ~ historical_joined$tmin * 0.97 - 0.53)

# reading in the predicted wind file an adding date information
historical_predictd_wind <- read.csv("./raw_data/weather/Estimated Pullman Historic Wind Speed.csv")
historical_predictd_wind$date <- as.Date(historical_predictd_wind$DATE, "%m/%d/%Y")
historical_predictd_wind$avgWspd_kmph[historical_predictd_wind$avgWspd_kmph == 0] <- mean(historical_predictd_wind$avgWspd_kmph, na.rm = TRUE)

# converting units of wind speed
historical_predictd_wind$avgWspd_kmph <- historical_predictd_wind$avgWspd_kmph/3.6 # converts from km/hour to m/s
historical_predictd_wind$year <- year(historical_predictd_wind$date)

historical_predictd_wind <- historical_predictd_wind %>%
  complete(date = seq.Date(min(date), max(date), by="day")) %>%
  fill(c('avgWspd_kmph')) %>%
  mutate(year = year(date),
         month = month(date),
         day = day(date),
         doy = yday(date)
  ) %>%
  select(-DATE) %>%
  filter(year >= min(historical_joined[, "year"]) & year <= max(historical_joined[, "year"]))

historical_joined <- historical_joined %>%
  mutate(
    date = ymd(paste(year, month, day))
  ) %>%
  complete(date = seq.Date(min(date), max(date), by="day"))

# combine historical_joined and historical_wind
historical_joined <- merge(historical_joined, historical_predictd_wind) %>%
  mutate(
    prcp = prcp/10
  )

# looping through each land cover and calculating the heat transfer resistance
# based on historical wind file
for (landcover in 1:length(land_cover_heights)) {
  
  height <- land_cover_heights[landcover]
  
  # zero plane displacement height
  d_value <- 0.65 * height
  
  # use zm_value instead of zm in the function call
  rh <- heat_transfer_resistance(
    zu = zu_value,
    zt = zt_value,
    k = k_value,
    d = d_value,
    zm = zm_value[landcover],
    zh = 0.2 * zm_value[landcover],
    u = historical_joined$avgWspd_kmph
  )
  
  historical_joined[, paste0("rh_", land_cover_names[landcover])] <- rh
  
}

# final renaming and formatting
historical_name_vector <- c(
  'date', 'year', 'month', 'day', 'doy', 'tmax', 'tmin', 'tavg', 'tdew', 'prcp', 'pet_hamon','Hour_1', 'Hour_6', 'Hour_12', 'Hour_18', 'l_turb', 'cloud', 'cc_water', 'cc_urban', 'cc_forest', 'cc_shrub', 'cc_grass', 'cc_row_crop', 'rh_water', 'rh_urban', 'rh_forest', 'rh_shrub', 'rh_grass', 'rh_row_crop', 'output' 
)

historical_joined <- historical_joined[, historical_name_vector] %>%
  rename(precip = prcp, pet = pet_hamon, hour_1 = Hour_1, hour_6 = Hour_6, hour_12 = Hour_12, hour_18 = Hour_18)

### ------------------------------------------------------- ###


### ----- WRITING WEATHER FILES BACK OUT TO THE 'weather' FOLDER ----- ###

start_date <- as.Date("1965-10-01")
end_date <- as.Date("1970-10-01")


historical_joined <- historical_joined[historical_joined$date >= start_date & historical_joined$date <= end_date, ]

# write out the Pullman gridMET file
write.table(gridMET_joined_p[1:365,], 
            file="./raw_data/weather/pullman_gridMET.csv", 
            col.names=FALSE, 
            row.names=FALSE
)

# write out the historical Pullman file
write.table(historical_joined, 
            file="./raw_data/weather/pullman_historical.csv", 
            col.names=FALSE, 
            row.names=FALSE,
            append = FALSE,
)


# write out the mini historical Pullman file
write.table(historical_joined[387,], 
            file="./raw_data/weather/pullman_historical_mini.csv", 
            col.names=FALSE, 
            row.names=FALSE
)

### ------------------------------------------------------- ###


d <- historical_joined

ggplot(d, aes(date,pet))+
  geom_point()

ggplot(d, aes(date,rh_row_crop)) +
  geom_point()

ggplot(d, aes(date,rh_shrub)) +
  geom_point()

ggplot(d, aes(date,rh_forest)) +
  geom_point()

# below this d is grouped

d %>%
  group_by(year) %>%
  summarise(sum_pet = sum(pet), tot_p = sum(precip)) %>%
  ggplot()+
  geom_point(aes(year,sum_pet))+
  geom_line(aes(year,tot_p))

ggplot(historical_joined) +
  geom_point(aes(date, tdew)) +
  geom_point(aes(date, tavg, color="tavg"))
  
