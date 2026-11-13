library(yaml)
library(tidyverse)
library(dplyr)
library(lubridate)
library(climateR)
library(sf)
library(chillR)
library(humidity)
library(bigleaf)

## ---- load project config ----
config_path   <- Sys.getenv("SMR_CONFIG", unset = "config.yaml")
cfg           <- yaml::read_yaml(config_path)
root_dir      <- Sys.getenv("SMR_ROOT_DIR", unset = cfg$root_dir)
setwd(root_dir)
site          <- cfg$weather$site
gridmet_start <- cfg$weather$gridmet_start
gridmet_end   <- cfg$weather$gridmet_end
wshed_lat     <- cfg$watershed$lat
wshed_lon     <- cfg$watershed$lon
raw_path      <- cfg$data$raw_dir
## ---------------------------

# A bounded interval makes a live API smoke test cheap and reproducible.
gridmet_start <- Sys.getenv("SMR_GRIDMET_START", unset = gridmet_start)
gridmet_end   <- Sys.getenv("SMR_GRIDMET_END", unset = gridmet_end)
if (as.Date(gridmet_end) < as.Date(gridmet_start)) {
  stop("SMR_GRIDMET_END precedes SMR_GRIDMET_START.")
}

### GRIDMET pull (modern period only — no historical station data required)
AOI_site <- sf::st_as_sf(
  data.frame(site = site, lon = wshed_lon, lat = wshed_lat),
  coords = c("lon", "lat"), crs = 4326
)

gridmet_variables <- c(
  prcp = "pr", tmax = "tmmx", tmin = "tmmn",
  pet_grass = "pet", srad = "srad", wind_vel = "vs"
)
gridmet_parts <- lapply(names(gridmet_variables), function(output_name) {
  source_name <- gridmet_variables[[output_name]]
  pulled <- as.data.frame(climateR::getGridMET(
    AOI_site, varname = source_name,
    startDate = gridmet_start, endDate = gridmet_end
  ))
  pulled <- pulled[, c("date", source_name)]
  names(pulled)[[2]] <- output_name
  pulled
})
gridMET_raw <- Reduce(function(x, y) dplyr::left_join(x, y, by = "date"), gridmet_parts)

gridMET <- gridMET_raw %>%
  mutate(
    date   = as.Date(date),
    doy    = yday(date),
    year   = lubridate::year(date),
    month  = lubridate::month(date),
    day    = lubridate::day(date),
    tmax   = tmax - 273.15,
    tmin   = tmin - 273.15,
    tavg   = (tmax + tmin) / 2,
    prcp   = prcp / 10,        # mm -> cm
    pet_grass = pet_grass / 10, # mm -> cm
    output = ifelse(month == 2, 1, 0),
    l_turb = 2.5
  )

### 6-hourly temperature disaggregation
data_for_hourly <- gridMET[, c("tmin","tmax","year","month","day")] %>%
  rename(Year = year, Month = month, Day = day, Tmin = tmin, Tmax = tmax)

hourly <- make_hourly_temps(latitude = wshed_lat, year_file = data_for_hourly,
                            keep_sunrise_sunset = FALSE) %>%
  rename(year = Year, month = Month, day = Day, tmin = Tmin, tmax = Tmax)
hourly <- hourly[, c("year","month","day","Hour_1","Hour_6","Hour_12","Hour_18")]

gridMET <- inner_join(gridMET, hourly, by = c("year","month","day"))

### Dew point (Gentilli 1955)
gridMET$tdew <- case_when(
  gridMET$tmin <= 0 ~ gridMET$tmin * 0.7 - 2,
  TRUE              ~ gridMET$tmin * 0.97 - 0.53)

### Cloud fraction (observed / potential solar)
julian_year <- seq(1, 365)
day_lengths  <- daylength(latitude = wshed_lat, JDay = julian_year, notimes.as.na = FALSE)

srad_pot <- sapply(julian_year, function(day) {
  hours <- seq(day_lengths$Sunrise[day], day_lengths$Sunset[day], by = 0.1)
  mean(potential.radiation(doy = day, hour = hours,
                           latDeg = wshed_lat, longDeg = wshed_lon,
                           timezone = -8, useSolartime = TRUE))
})
srad_df <- data.frame(doy = julian_year, srad_potential = srad_pot)

gridMET <- inner_join(gridMET, srad_df, by = "doy")
gridMET$cloud <- gridMET$srad / gridMET$srad_potential
gridMET$cloud <- pmin(pmax(gridMET$cloud, 0), 1)

### Crop coefficients (FAO-based seasonal curves)
ann_Kc_curve <- function(time_params, Kc_params) {
  preplant  <- rep(Kc_params$Kc_preplant,  time_params$preplant_len)
  init_sec  <- rep(Kc_params$Kc_init,      time_params$init_len)
  dev_sec   <- round(seq(Kc_params$Kc_init, Kc_params$Kc_mid, length.out = time_params$dev_len), 4)
  mid_sec   <- rep(Kc_params$Kc_mid,       time_params$mid_len)
  late_sec  <- round(seq(Kc_params$Kc_mid, Kc_params$Kc_end, length.out = time_params$late_len), 4)
  post_sec  <- rep(Kc_params$Kc_post_harv, time_params$post_len)
  c(preplant, init_sec, dev_sec, mid_sec, late_sec, post_sec)
}

water_curve    <- data.frame(cc_water    = ann_Kc_curve(list(preplant_len=60,init_len=60,dev_len=60,mid_len=60,late_len=60,post_len=65), list(Kc_preplant=0.6525,Kc_init=0.6525,Kc_mid=0.6525,Kc_end=0.6525,Kc_post_harv=0.6525)), doy=1:365)
urban_curve    <- data.frame(cc_urban    = ann_Kc_curve(list(preplant_len=60,init_len=60,dev_len=60,mid_len=60,late_len=60,post_len=65), list(Kc_preplant=0.30,Kc_init=0.30,Kc_mid=0.30,Kc_end=0.30,Kc_post_harv=0.30)),         doy=1:365)
forest_curve   <- data.frame(cc_forest   = ann_Kc_curve(list(preplant_len=60,init_len=60,dev_len=60,mid_len=60,late_len=60,post_len=65), list(Kc_preplant=1.0,Kc_init=1.0,Kc_mid=1.0,Kc_end=1.0,Kc_post_harv=1.0)),               doy=1:365)
shrub_curve    <- data.frame(cc_shrub    = ann_Kc_curve(list(preplant_len=60,init_len=20,dev_len=70,mid_len=90,late_len=30,post_len=95),  list(Kc_preplant=0.45,Kc_init=0.45,Kc_mid=0.95,Kc_end=0.70,Kc_post_harv=0.45)),          doy=1:365)
grass_curve    <- data.frame(cc_grass    = ann_Kc_curve(list(preplant_len=113,init_len=10,dev_len=20,mid_len=64,late_len=62,post_len=96), list(Kc_preplant=0.30,Kc_init=0.30,Kc_mid=0.75,Kc_end=0.75,Kc_post_harv=0.30)),          doy=1:365)

# row crop: shift so that spring planting aligns to ~DOY 91 (April 1)
row_crop_raw   <- ann_Kc_curve(list(preplant_len=0,init_len=160,dev_len=75,mid_len=75,late_len=25,post_len=30), list(Kc_preplant=0.30,Kc_init=0.7,Kc_mid=1.15,Kc_end=0.35,Kc_post_harv=0.30))
shift_by       <- 91                                  # shift so init starts ~April 1
doy_shifted    <- ((seq_len(365) - 1 + shift_by) %% 365) + 1   # keeps values in 1-365
row_crop_curve <- data.frame(cc_row_crop = row_crop_raw, doy = doy_shifted)
row_crop_curve <- row_crop_curve[order(row_crop_curve$doy), ]   # sort so join works

gridMET <- gridMET %>%
  inner_join(water_curve,    by = "doy") %>%
  inner_join(urban_curve,    by = "doy") %>%
  inner_join(forest_curve,   by = "doy") %>%
  inner_join(shrub_curve,    by = "doy") %>%
  inner_join(grass_curve,    by = "doy") %>%
  inner_join(row_crop_curve, by = "doy")

### Heat transfer resistance (Campbell 1977)
zu_value <- 2; zt_value <- 1; k_value <- 0.41
land_cover_names   <- c("water","urban","forest","shrub","grass","row_crop")
land_cover_heights <- c(0.01, 0.01, 1.3, 1.3, 1.3, 1.3)
zm_value <- 0.1 * land_cover_heights

heat_transfer_resistance <- function(zu, zt, k, d, zm, zh, u) {
  (log((zu - d + zm) / zm) * log((zt - d + zh) / zh)) / (k^2 * u)
}

for (lc in seq_along(land_cover_heights)) {
  d_val <- 0.65 * land_cover_heights[lc]
  rh    <- heat_transfer_resistance(zu_value, zt_value, k_value,
                                    d_val, zm_value[lc], 0.2 * zm_value[lc],
                                    gridMET$wind_vel)
  gridMET[[paste0("rh_", land_cover_names[lc])]] <- rh
}

### Final column selection and rename
out_cols <- c(
  "date","year","month","day","doy",
  "tmax","tmin","tavg","tdew",
  "prcp","pet_grass",
  "Hour_1","Hour_6","Hour_12","Hour_18",
  "l_turb","cloud",
  "cc_water","cc_urban","cc_forest","cc_shrub","cc_grass","cc_row_crop",
  "rh_water","rh_urban","rh_forest","rh_shrub","rh_grass","rh_row_crop",
  "output"
)

gridMET_out <- gridMET[, out_cols] %>%
  rename(precip = prcp, pet = pet_grass,
         hour_1 = Hour_1, hour_6 = Hour_6, hour_12 = Hour_12, hour_18 = Hour_18)

### Write weather file
weather_dir  <- file.path(raw_path, "weather")
if (!dir.exists(weather_dir)) dir.create(weather_dir, recursive = TRUE)
out_path <- file.path(weather_dir, "gridmet_weather.csv")

write.table(gridMET_out,
            file      = out_path,
            col.names = FALSE,
            row.names = FALSE)

message("Weather file written to ", out_path)
message("Rows: ", nrow(gridMET_out), " | Period: ", min(gridMET_out$date), " to ", max(gridMET_out$date))
