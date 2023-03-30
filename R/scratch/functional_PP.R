# import libraries
library(tidyverse)
library(dplyr)
library(lubridate)
library(hydroGOF)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(cowplot)


# data pre-processing and alignment function
preprocessing <- function(
    validation_path, 
    modeled_path, 
    start_date, 
    end_date, 
    modeled_headers
    ) 
  {
  
  validation_data <- read.csv(validation_path) %>% 
    mutate(
      mean_discharge_cms = mean_discharge_cfs*0.028316832,
      date = as.Date(date, "%m/%d/%Y")
      )

  modeled_data <- read.csv(modeled_path, sep=' ', header = FALSE, col.names = modeled_headers) %>% 
    mutate(date = as.Date(date))
  
  common_start <- max(min(validation_data$date), min(modeled_data$date))
  common_end <- min(max(validation_data$date), max(modeled_data$date))

  validation_data <- validation_data[validation_data$date >= common_start & validation_data$date <= common_end, ]
  
  modeled_data <- modeled_data[modeled_data$date >= common_start & modeled_data$date <= common_end, ]
  
  validation_data <- validation_data[, c('date', 'mean_discharge_cms')] %>%
    rename(validation_Q = mean_discharge_cms)
  
  modeled_data <- merge(modeled_data, validation_data)

  return(modeled_data)
  
  }


# stream flow comparison function
Q_comparison <- function(combined_data) {
  
  comparison_plot <- ggplot(data=combined_data) +
    geom_line(aes(x=date, y=validation_Q, color='Validation')) +
    geom_line(aes(x=date, y=Q, color='Modeled')) +
    ggtitle('Observed Versus Simulated Streamflow (with precip, rain, and snow)') +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.title.x=element_blank()) +
    theme(axis.text.x=element_blank()) +
    ylab('Streamflow (m3/s)') +
    scale_y_continuous(trans='log10') +
    guides(color=guide_legend(title='Q Type'))

  precip_plot <- ggplot(data=combined_data) +
    geom_line(aes(x=date, y=precip_cm, color='Precipitation')) +
    geom_line(aes(x=date, y=rain_cm, color='Rain')) +
    geom_line(aes(x=date, y=snow_cm, color='Snow')) +
    xlab('Date') +
    ylab('Precipitation, Rain, and Snow') +
    guides(color=guide_legend(title='Precip Type')) +
    scale_y_continuous(trans='log10')
  
  cowplot::plot_grid(comparison_plot, precip_plot, align = "v", ncol = 1, rel_heights = c(0.60, 0.40))
}


# function for plotting important components of SAM
SAM_check <- function(combined_data) {
  
  ggplot(data=combined_data) +
    scale_y_continuous(trans='log10') +
    geom_line(aes(x=date, y=snow_cm, color='Snow')) +
    geom_line(aes(x=date, y=swe_cm, color='SWE')) +
    geom_line(aes(x=date, y=snowmelt_cm, color='Snowmelt')) + 
    geom_line(aes(x=date, y=condens_cm, color='Condens')) +
    ggtitle('Relevant Snowmelt and Accumulation Variables') +
    theme(plot.title = element_text(hjust = 0.5)) +
    guides(color=guide_legend(title='Precip Type'))
  
}

# NSE function for entire period
nse_Q <- function(combined_data) {
  return(NSE(combined_data$Q, combined_data$validation_Q))
}

# KGE function for entire period
kge_Q <- function(combined_data) {
  return(KGE(combined_data$Q, combined_data$validation_Q))
}

# flux time-series function (looped through in flux_ts_loop())
flux_ts <- function(combined_data, flux) {
  
  flux_series <- combined_data[, flux] / 79384
  date_series <- combined_data[, 'date']
  max_date <- max(date_series)
  min_date <- min(date_series)
  
  ggplot() +
    geom_line(aes(x=date_series, y=flux_series)) +
    ggtitle(
      paste0('Basin-Averaged ',paste0(
        paste0(
          paste0(
            paste0(
              flux, ' from '), min_date), ' to '), max_date))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab(paste0(flux, ' (cm)')) +
    xlab(paste0('Date'))
  
}

# produces individual time series of fluxes by looping the flux_ts() function
flux_ts_loop <- function(combined_data, fluxes) {
  
  combined_date <- combined_data
  
  for (flux in fluxes) {
    print(flux_ts(combined_data, flux))
  }
  
}

# produces bar plot of the accumulated annual fluxes
annual_fluxes <- function(combined_data, fluxes) {
  
  fluxes <- c(fluxes, 'date', 'year')
  combined_data <- combined_data[, fluxes] %>% 
    pivot_longer(cols=-c(date, year), names_to = "variable", values_to = "value") %>%
    mutate(value = value) %>%
    group_by(year, variable) %>%
    summarize(sum_value = sum(value, na.rm=TRUE))
  
  ggplot(data = combined_data, aes(x = year, y = sum_value, group = variable, color = variable, fill=variable)) +
    geom_bar(position='dodge', stat='identity') + 
    scale_y_continuous(trans='log10') +
    xlab("Year") +
    ylab("Total Value") +
    #ylim(min(combined_data$sum_value), max(combined_data$sum_value)) +
    ggtitle("Annual Accumulateed Fluxes") +
    theme(plot.title = element_text(hjust = 0.5))

}

# produces a time series of the different radiation components
radiation_ts <- function(combined_data) {
  
  combined_data <- combined_data %>%
    mutate(
      srad = srad/79384,
      latent = latent/79384,
      sensible = sensible/79384,
      lw = lw/79384,
      q_rain_ground
    )
  
  ggplot(data=combined_data) +
    #scale_y_continuous(trans='log10') +
    geom_line(aes(x=date, y=srad, color='Srad')) +
    geom_line(aes(x=date, y=latent, color='Latent')) +
    geom_line(aes(x=date, y=sensible, color='Sensible')) + 
    geom_line(aes(x=date, y=lw, color='Longwave')) +
    geom_line(aes(x=date, y=q_rain_ground, color='Rain Ground')) +
    ggtitle('Components of Radiation Budget') +
    theme(plot.title = element_text(hjust = 0.5)) +
    guides(color=guide_legend(title='Precip Type'))
  
}

# annual mass balance function
annual_mass_balance <- function(combined_data) {
  
  annual_combined <- combined_data %>%
    group_by(year) %>%
    select(-c(date, wshed_id)) %>%
    summarize_all(.funs=sum)

  annual_combined$mass_balance <-
    annual_combined$baseflow +
    annual_combined$precip_cm - 
    annual_combined$storage_amt_cm - 
    annual_combined$canopy_storage_amt_cm - 
    annual_combined$swe_cm - 
    annual_combined$canopy_evap_cm
  
  ggplot(data=annual_combined) +
    #scale_y_continuous(trans='log10') +
    geom_line(aes(x=year, y=mass_balance)) +
    ggtitle('Annual Mass Balance') +
    theme(plot.title = element_text(hjust = 0.5))
}

# calling functions all above functions

flux_vector <- c('runoff_cm', 'precip_cm', 'rain_cm', 
                 'actualET_flow_cm', 'canopy_evap_cm', 
                 'snowmelt_cm', 'storage_amt_cm', 
                 'throughfall_cm', 'canopy_storage_amt_cm', 
                 'perc_cm', 'swe_cm', 'condens_cm', 
                 'snow_cm', 'baseflow')

VaM_data <- preprocessing(
  "./raw_data/discharge_historical/USGSdischarge.csv",
  "./raw_data/smr_output/MFC_mass_balance_79_1.csv",
  "1965-10-01",
  "1970-10-01",
  c('wshed_id', 'date', 'year', 'runoff_cm', 'precip_cm', 'rain_cm', 'actualET_flow_cm', 'canopy_evap_cm', 'snowmelt_cm', 'storage_amt_cm', 'throughfall_cm', 'canopy_storage_amt_cm', 'perc_cm', 'Q', 'swe_cm', 'condens_cm', 'snow_cm', 'baseflow', 'srad', 'latent', 'sensible', 'lw', 'q_rain_ground', 'q_total')
)

NSE_Q <- nse_Q(VaM_data)
KGE_Q <- kge_Q(VaM_data)

annual_mass_balance(VaM_data)
flux_ts_loop(VaM_data, flux_vector)

annual_fluxes(VaM_data, flux_vector)
radiation_ts(VaM_data)
SAM_check(VaM_data)
Q_comparison(VaM_data)

