# sourcing the functional post processing script
source("./R/scratch/functional_PP.R")

# fluxes to be included in the single flux time series plots
flux_vector <- c('runoff_cm', 'precip_cm', 'rain_cm', 
                 'actualET_flow_cm', 'canopy_evap_cm', 
                 'snowmelt_cm', 'storage_amt_cm', 
                 'throughfall_cm', 'canopy_storage_amt_cm', 
                 'perc_cm', 'swe_cm', 'condens_cm', 
                 'snow_cm', 'baseflow')


# date of the model run
model_run_date <- "2023-04-12"

VaM_data <- preprocessing(
  "./raw_data/discharge_historical/USGSdischarge.csv",
  "./raw_data/smr_output/mfc_mb_2023-04-12.csv",
  "1965-10-01",
  "1970-10-01",
  TRUE,
  model_run_date,
  c('wshed_id', 'date', 'year', 'runoff_cm', 'precip_cm', 'rain_cm', 'actualET_flow_cm', 'canopy_evap_cm', 'snowmelt_cm', 'storage_amt_cm', 'throughfall_cm', 'canopy_storage_amt_cm', 'perc_cm', 'Q', 'swe_cm', 'condens_cm', 'snow_cm', 'baseflow', 'srad', 'latent', 'sensible', 'lw', 'q_rain_ground', 'q_total', 'ice_content', 'liquid_water', 'refreeze', 'vap_d_air', 'vap_d_snow', 'u_surface', 'empty')
)

flux_vector <- c( 'runoff_cm', 'precip_cm', 'rain_cm', 'actualET_flow_cm', 'canopy_evap_cm', 'snowmelt_cm', 'storage_amt_cm', 'throughfall_cm', 'canopy_storage_amt_cm', 'perc_cm', 'Q', 'swe_cm', 'condens_cm', 'snow_cm', 'baseflow', 'srad', 'latent', 'sensible', 'lw', 'q_rain_ground', 'q_total', 'ice_content', 'liquid_water', 'refreeze', 'vap_d_air', 'vap_d_snow', 'u_surface')

# all post processing functions
NSE_Q <- nse_Q(VaM_data)
KGE_Q <- kge_Q(VaM_data)

annual_mass_balance(VaM_data)
flux_ts_loop(VaM_data, flux_vector, log_transform = F)

annual_fluxes(VaM_data, flux_vector, log_transform = F)
radiation_ts(VaM_data, log_transform = F)
SAM_check(VaM_data, log_transform = T)
swe_debug(VaM_data, log_transform = F)
q.latent_debug(VaM_data, log_transform = F)
Q_comparison(VaM_data, log_transform = T)

get_map_outputs(
  '/Users/duncanjurayj/Documents/SMR_R/raw_data/smr_output/map_outputs/feb_outputs',
  model_run_date
)

storage_amt_feb(`avg_A_amt_feb_1969_2023-04-12`, 'A')





