# calling functions all above functions
# source("./R/load_libs.R") # do this. 
source("./R/scratch/functional_PP.R")


flux_vector <- c('runoff_cm', 'precip_cm', 'rain_cm', 
                 'actualET_flow_cm', 'canopy_evap_cm', 
                 'snowmelt_cm', 'storage_amt_cm', 
                 'throughfall_cm', 'canopy_storage_amt_cm', 
                 'perc_cm', 'swe_cm', 'condens_cm', 
                 'snow_cm', 'baseflow')

getwd() 


VaM_data <- preprocessing(
  "./raw_data/discharge_historical/USGSdischarge.csv",
  "./raw_data/smr_output/mfc_mb_2023-03-31.csv",
  "1965-10-01",
  "1970-10-01",
  c('wshed_id', 'date', 'year', 'runoff_cm', 'precip_cm', 'rain_cm', 'actualET_flow_cm', 'canopy_evap_cm', 'snowmelt_cm', 'storage_amt_cm', 'throughfall_cm', 'canopy_storage_amt_cm', 'perc_cm', 'Q', 'swe_cm', 'condens_cm', 'snow_cm', 'baseflow', 'srad', 'latent', 'sensible', 'lw', 'q_rain_ground', 'q_total', 'ice_content', 'liquid_water', 'refreeze', 'vap_d_air', 'vap_d_snow')
)

NSE_Q <- nse_Q(VaM_data)
KGE_Q <- kge_Q(VaM_data)

annual_mass_balance(VaM_data)
flux_ts_loop(VaM_data, flux_vector)

annual_fluxes(VaM_data, flux_vector)
radiation_ts(VaM_data)
SAM_check(VaM_data)
Q_comparison(VaM_data)

