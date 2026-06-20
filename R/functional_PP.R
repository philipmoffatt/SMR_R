library(yaml)
library(tidyverse)
library(dplyr)
library(lubridate)
library(hydroGOF)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(cowplot)

## ---- load project config ----
cfg             <- yaml::read_yaml("config.yaml")
setwd(cfg$root_dir)
discharge_file  <- cfg$validation$discharge_file
smr_output_file <- cfg$validation$smr_output_file
## ---------------------------

modeled_headers <- c(
  "wshed_id","date","year",
  "runoff_cm","precip_cm","rain_cm",
  "actualET_flow_cm","canopy_evap_cm","snowmelt_cm",
  "storage_amt_cm","throughfall_cm","canopy_storage_amt_cm",
  "perc_cm","Q",
  "swe_cm","condens_cm","snow_cm","baseflow",
  "srad","latent","sensible","lw","q_rain_ground","q_total"
)

preprocessing <- function(validation_path, modeled_path, start_date, end_date, headers) {
  validation_data <- read.csv(validation_path) %>%
    mutate(
      mean_discharge_cms = mean_discharge_cfs * 0.028316832,
      date = as.Date(date, "%m/%d/%Y"))

  modeled_data <- read.csv(modeled_path, sep = " ", header = FALSE,
                           col.names = headers) %>%
    mutate(date = as.Date(date))

  common_start <- max(min(validation_data$date), min(modeled_data$date))
  common_end   <- min(max(validation_data$date), max(modeled_data$date))

  validation_data <- validation_data %>%
    filter(date >= common_start, date <= common_end) %>%
    select(date, validation_Q = mean_discharge_cms)

  modeled_data <- modeled_data %>%
    filter(date >= common_start, date <= common_end) %>%
    merge(validation_data)

  modeled_data
}

Q_comparison <- function(d) {
  p1 <- ggplot(d) +
    geom_line(aes(date, validation_Q, color = "Observed")) +
    geom_line(aes(date, Q, color = "Modeled")) +
    scale_y_continuous(trans = "log10") +
    ylab("Streamflow (m3/s)") +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
    guides(color = guide_legend(title = "Q type"))
  p2 <- ggplot(d) +
    geom_line(aes(date, precip_cm, color = "Precipitation")) +
    geom_line(aes(date, rain_cm,   color = "Rain")) +
    geom_line(aes(date, snow_cm,   color = "Snow")) +
    scale_y_continuous(trans = "log10") +
    ylab("Depth (cm)") + xlab("Date") +
    guides(color = guide_legend(title = "Precip type"))
  cowplot::plot_grid(p1, p2, align = "v", ncol = 1, rel_heights = c(0.6, 0.4))
}

SAM_check <- function(d) {
  ggplot(d) +
    scale_y_continuous(trans = "log10") +
    geom_line(aes(date, snow_cm,    color = "Snow")) +
    geom_line(aes(date, swe_cm,     color = "SWE")) +
    geom_line(aes(date, snowmelt_cm,color = "Snowmelt")) +
    geom_line(aes(date, condens_cm, color = "Condensation")) +
    ggtitle("Snow Accumulation and Melt Variables") +
    theme(plot.title = element_text(hjust = 0.5)) +
    guides(color = guide_legend(title = "Variable"))
}

nse_Q <- function(d) NSE(d$Q, d$validation_Q)
kge_Q <- function(d) KGE(d$Q, d$validation_Q)

flux_ts <- function(d, flux) {
  n_cells   <- nrow(d) # placeholder — divide by watershed cell count if known
  vals      <- d[[flux]]
  ggplot() +
    geom_line(aes(x = d$date, y = vals)) +
    ggtitle(paste0("Basin-Averaged ", flux)) +
    ylab(paste0(flux, " (cm)")) + xlab("Date") +
    theme(plot.title = element_text(hjust = 0.5))
}

flux_ts_loop <- function(d, fluxes) {
  for (flux in fluxes) print(flux_ts(d, flux))
}

annual_fluxes <- function(d, fluxes) {
  d %>%
    select(all_of(c(fluxes, "date", "year"))) %>%
    pivot_longer(cols = -c(date, year), names_to = "variable", values_to = "value") %>%
    group_by(year, variable) %>%
    summarise(sum_value = sum(value, na.rm = TRUE), .groups = "drop") %>%
    ggplot(aes(year, sum_value, color = variable, fill = variable)) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_y_continuous(trans = "log10") +
    xlab("Year") + ylab("Total (cm)") +
    ggtitle("Annual Accumulated Fluxes") +
    theme(plot.title = element_text(hjust = 0.5))
}

radiation_ts <- function(d) {
  # scale from summed cell values to basin average
  n_cells <- max(d$wshed_id)  # rough proxy; replace with actual cell count if known
  d2 <- d %>% mutate(
    srad     = srad     / n_cells,
    latent   = latent   / n_cells,
    sensible = sensible / n_cells,
    lw       = lw       / n_cells,
    q_rain_ground = q_rain_ground / n_cells)
  ggplot(d2) +
    geom_line(aes(date, srad,          color = "Shortwave")) +
    geom_line(aes(date, latent,        color = "Latent")) +
    geom_line(aes(date, sensible,      color = "Sensible")) +
    geom_line(aes(date, lw,            color = "Longwave")) +
    geom_line(aes(date, q_rain_ground, color = "Rain-ground")) +
    ggtitle("Components of Radiation Budget") +
    theme(plot.title = element_text(hjust = 0.5)) +
    guides(color = guide_legend(title = "Component"))
}

# Annual mass balance: P = Q + ET + ΔS  (all terms in cm, basin-summed)
annual_mass_balance <- function(d) {
  d %>%
    group_by(year) %>%
    select(-c(date, wshed_id)) %>%
    summarise_all(sum) %>%
    mutate(
      total_ET     = actualET_flow_cm + canopy_evap_cm,
      delta_S      = storage_amt_cm + canopy_storage_amt_cm + swe_cm,
      mass_balance = precip_cm - runoff_cm - total_ET - perc_cm - delta_S
    ) %>%
    ggplot(aes(year, mass_balance)) +
    geom_line() + geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ggtitle("Annual Mass Balance (should be near 0)") +
    ylab("Residual (cm)") + xlab("Year") +
    theme(plot.title = element_text(hjust = 0.5))
}

# ---- Run ----
flux_vector <- c("runoff_cm","precip_cm","rain_cm",
                 "actualET_flow_cm","canopy_evap_cm",
                 "snowmelt_cm","storage_amt_cm",
                 "throughfall_cm","canopy_storage_amt_cm",
                 "perc_cm","swe_cm","condens_cm",
                 "snow_cm","baseflow")

VaM_data <- preprocessing(
  validation_path = discharge_file,
  modeled_path    = smr_output_file,
  start_date      = NULL,
  end_date        = NULL,
  headers         = modeled_headers)

cat("NSE:", nse_Q(VaM_data), "\n")
cat("KGE:", kge_Q(VaM_data), "\n")

annual_mass_balance(VaM_data)
flux_ts_loop(VaM_data, flux_vector)
annual_fluxes(VaM_data, flux_vector)
radiation_ts(VaM_data)
SAM_check(VaM_data)
Q_comparison(VaM_data)
