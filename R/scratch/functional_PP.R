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
library(raster)
library(mapview)

# data pre-processing and alignment function
version_tracker <- function(modeled_path, simulation_note = "") {
  outputDir <- dirname(modeled_path) %>% gsub(., pattern = "raw_data", replacement = "processed_data")
  
  if(!dir.exists(outputDir)) {
    dir.create(outputDir)
  }
  
  SMR_log_file <- file.path(outputDir, "SMR_run_log.txt")
  
  if(simulation_note == "") {
    simulation_note <- readline(prompt="Simulation description: ")
  }
  
  if(!file.exists(SMR_log_file)){
    runID <- "1"
  } else {
    table_temp <- (read.table(file= SMR_log_file) %>% row.names())
    runID <- table_temp[length(table_temp)] %>% as.integer() + 1
  }  
  
  runID <- paste0("MFC_", rep(0,(3-nchar(runID))) %>% paste0(., collapse = ""), runID)
  
  log_entry <- list(runID, Sys.Date(), simulation_note %>% as.character())
  
  write.table(x = log_entry, file = SMR_log_file, append = TRUE, quote = TRUE, col.names = FALSE,
              row.names = FALSE)
  
  
  out_path <- file.path(outputDir, paste0(runID, "_", Sys.Date()))
  
  if(!dir.exists(out_path)) {
    dir.create(out_path)
  }
  
  map_path <- file.path(out_path, "maps") 
  
  if(!dir.exists(map_path)) {
    dir.create(map_path)
  }
  
  return(out_path)
  
}

# allows for dynamic text based on model outputs in the markdown
get_run_dates <- function(weather_data_path) {
  
  weather_data <- read.csv(data_path)
  min_date <- min(weather_data$date)
  max_date <- max(weather_data$date)
  date_range <- c(as.character(min_date), as.character(max_date))
    
}

get_print_out_line <- function(perl_script_path) {
  # Read entire Perl script into a character vector of lines
  perl_script_lines <- readLines(perl_script_path)
  
  # find index of line containing print OUT statement
  print_out_line_index <- grep('print\\s+OUT', perl_script_lines)
  
  if (length(print_out_line_index) > 0) {
    # get full line containing print OUT statement
    full_line <- perl_script_lines[print_out_line_index]
    
    # extract part of line starting after "print OUT" and remove ' \n";'
    partial_line <- gsub('^.*?print\\s+OUT\\s+"', '', full_line)
    cleaned_line_1 <- gsub('\\s*\\\\n\";$', '', partial_line)
    
    # remove "\\" before "$wshed_id" and remove "_{$wshed_id}" from each word in the line
    cleaned_line_2 <- gsub('\\\\\\$', '', cleaned_line_1)
    final_cleaned_string <- gsub('_\\{\\$wshed_id\\}', '', cleaned_line_2)
    
    # split final_cleaned_string on space to obtain a list of variables, removing the "$" before each variable name
    variable_list <- strsplit(final_cleaned_string, ' ')[[1]] %>% 
      gsub('\\$', '', .)
    
    return(variable_list)
    
  } else {
    stop("No matching print OUT line was found.")
  }
}

preprocessing <- function(
    validation_path, 
    modeled_path, 
    modeled_headers,
    version_tracked_outpath
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
  
  out_path = file.path(version_tracked_outpath, 'mfc_mb.csv')
  write.csv(x = modeled_data, file = out_path)

  file.remove(modeled_path)
  
  return(modeled_data)
  
  }


# stream flow comparison function --> baseflow is subtracted out right now 
# stream flow comparison function 
Q_comparison <- function(combined_data, log_transform=FALSE) {
  
  if (log_transform) {
    y_trans <- 'log10'
  } else {
    y_trans <- 'identity'
  }
  
  combined_data <- combined_data %>%
    mutate(date = as.Date(date))
  
  comparison_plot <- ggplot(data=combined_data) +
    geom_line(aes(x=date, y=validation_Q, color='Validation')) +
    geom_line(aes(x=date, y=Q, color='Modeled')) +
    ggtitle('Observed Versus Simulated Streamflow (with precip, rain, and snow)') +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.title.x=element_blank()) +
    theme(axis.text.x=element_blank()) +
    ylab(expression(paste(
      "Streamflow (",m^3, "/", s,
      ")", sep=""))) +
    scale_y_continuous(trans=y_trans,labels = scales::comma) +
    guides(color=guide_legend(title='Streamflow'))
  
  precip_plot <- ggplot(data=combined_data) +
    geom_line(aes(x=date, y=precip_cm, color='Precipitation')) +
    geom_line(aes(x=date, y=rain_cm, color='Rain')) +
    geom_line(aes(x=date, y=snow_cm, color='Snow')) +
    xlab('Date') +
    ylab('Precipitation, Rain, and Snow') + 
    scale_x_date(date_breaks = "1 year") + 
    guides(color=guide_legend(title='Precip Type'))+
    scale_y_continuous(labels = scales::comma)
  
  cowplot::plot_grid(comparison_plot, precip_plot, align = "v", ncol = 1, rel_heights = c(0.60, 0.40))
}


# function for plotting important components of SAM
SAM_check <- function(combined_data, log_transform=FALSE) {
  
  if (log_transform) {
    y_trans <- 'log10'
  } else {
    y_trans <- 'identity'
  }
  
  ggplot(data=combined_data) +
    geom_line(aes(x=date, y=snow_cm, color='Snow')) +
    geom_line(aes(x=date, y=swe_cm, color='SWE')) +
    geom_line(aes(x=date, y=snowmelt_cm, color='Snowmelt')) + 
    geom_line(aes(x=date, y=condens_cm, color='Condens')) +
    ggtitle('Relevant Snowmelt and Accumulation Variables') +
    theme(plot.title = element_text(hjust = 0.5)) +
    guides(color=guide_legend(title='Precip Type')) + 
    scale_y_continuous(trans=y_trans)
  
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
flux_ts <- function(combined_data, flux, log_transform = FALSE) {
  
  if (log_transform) {
    y_trans <- 'log10'
  } else {
    y_trans <- 'identity'
  }
  
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
    xlab(paste0('Date')) + 
    scale_y_continuous(trans=y_trans)
  
}

# produces individual time series of fluxes by looping the flux_ts() function
flux_ts_loop <- function(combined_data, fluxes, log_transform=FALSE) {
  
  for (flux in fluxes) {
    print(flux_ts(combined_data, flux, log_transform = log_transform))
  }
  
}

# produces bar plot of the accumulated annual fluxes
annual_fluxes <- function(combined_data, fluxes, log_transform = FALSE) {
  
  if (log_transform) {
    y_trans <- 'log10'
  } else {
    y_trans <- 'identity'
  }
  
  fluxes <- c(fluxes, 'date', 'year')
  combined_data <- combined_data[, fluxes] %>% 
    pivot_longer(cols=-c(date, year), names_to = "variable", values_to = "value") %>%
    mutate(value = value) %>%
    group_by(year, variable) %>%
    summarize(sum_value = sum(value, na.rm=TRUE))
  
  ggplot(data = combined_data, aes(x = year, y = sum_value, group = variable, color = variable, fill=variable)) +
    geom_bar(position='dodge', stat='identity') + 
    scale_y_continuous(trans=y_trans) +
    xlab("Year") +
    ylab("Total Value") +
    ggtitle("Annual Accumulateed Fluxes") +
    theme(plot.title = element_text(hjust = 0.5))

}

# produces a time series of the different radiation components
radiation_ts <- function(combined_data, log_transform = FALSE) {
  
  # Set log transform for y-axis scale
  if (log_transform) {
    y_trans <- 'log10'
  } else {
    y_trans <- 'identity'
  }
  
  # Divide all values in selected columns by 79384
  combined_data <- combined_data %>%
    mutate(
      srad = srad/79384,
      latent = latent/79384,
      sensible = sensible/79384,
      lw = lw/79384,
      q_rain_ground_cm
    )
  
  ggplot(data=combined_data) +
    geom_line(aes(x=date, y=srad, color='Srad')) +
    geom_line(aes(x=date, y=latent, color='Latent')) +
    geom_line(aes(x=date, y=sensible, color='Sensible')) + 
    geom_line(aes(x=date, y=lw, color='Longwave')) +
    geom_line(aes(x=date, y=q_rain_ground_cm, color='Rain Ground')) +
    scale_y_continuous(trans=y_trans) + 
    ggtitle('Components of Radiation Budget') +
    theme(plot.title = element_text(hjust = 0.5)) +
    guides(color=guide_legend(title='Precip Type'))
  
}
# annual mass balance function
annual_mass_balance <- function(combined_data, log_transform = FALSE) {
  
  # Set log transform for y-axis scale
  if (log_transform) {
    y_trans <- 'log10'
  } else {
    y_trans <- 'identity'
  }
  
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
    scale_y_continuous(trans=y_trans) +
    geom_line(aes(x=year, y=mass_balance)) +
    ggtitle('Annual Mass Balance') +
    theme(plot.title = element_text(hjust = 0.5))
}

# function just for determing if ice.content or liquid.water is contributing to
# swe inflation
swe_debug <- function(combined_data, log_transform = FALSE) {
  
  # Set log transform for y-axis scale
  if (log_transform) {
    y_trans <- 'log10'
  } else {
    y_trans <- 'identity'
  }
  
  ggplot(data=combined_data) +
    geom_line(aes(x=date, y=ice_content, color='Ice Content')) +
    geom_line(aes(x=date, y=liquid_water, color='Liquid Water')) +
    geom_line(aes(x=date, y=refreeze, color='Refreeze')) +
    geom_line(aes(x=date, y=snow_cm, color='Snow')) +
    geom_line(aes(x=date, y=swe_cm, color='SWE')) + 
    scale_y_continuous(trans=y_trans) + 
    ggtitle('Snow Water Equivalent and its Components') +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab('cm') + 
    guides(color=guide_legend(title='SWE Variables'))
  
}


# function to confirm that value inflation comes from tiny rh rather than large
# vap.d.air which is in the numerator
q.latent_debug <- function(combined_data, log_transform = FALSE) {
  
  # Set log transform for y-axis scale
  if (log_transform) {
    y_trans <- 'log10'
  } else {
    y_trans <- 'identity'
  }
  
  ggplot(data=combined_data) +
    geom_line(aes(x=date, y=vap_d_air, color='Vap_d_air')) +
    geom_line(aes(x=date, y=vap_d_snow, color='Vap_d_snow')) +
    scale_y_continuous(trans=y_trans) + 
    ggtitle('Size of vap.d.air compared to vap.d.snow') +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab('Vapor Density (kg/m^3)') + 
    guides(color=guide_legend(title='Vapor Density Variable'))
  
}

# this should take in map outputs a date of the model run and 
# output those 
# TO ADD: needs a feature that only adds a date if a date already exists
#get_map_outputs <- function(map_dir, model_run_date) {
  
 # tif_files <- list.files(map_dir, pattern = ".tif$")
  
#  for (file in tif_files) {
    
 #   var_name <- gsub(".tif$", "", file)
  #  print(var_name)
    
  #  if (exists(var_name)) {
  #    next
  #  }
    
  #  if (grepl(model_run_date, file)) {
  #    new_file_name <- file
  #  } else {
  #    new_file_name <- paste0(var_name, "_", model_run_date, ".tif")
  #    file.rename(paste0(map_dir, "/", file), paste0(map_dir, "/", new_file_name))
  #  }
  #  
  #  if (grepl(model_run_date, new_file_name)) {
  #    assign(var_name, raster(paste0(map_dir, "/", new_file_name)), envir=globalenv())
  #  }
    
#  }
#}

get_map_outputs <- function(raw_map_dir, version_tracked_map_dir){
  
  file.copy(from = raw_map_dir, to = version_tracked_map_dir)
  
  file.remove(list.files(path = raw_map_dir, pattern = ".tif", full.names = TRUE))
  
  list_of_files <- list.files(version_tracked_map_dir, pattern=".tif", full.names=TRUE)
  
  for (file in list_of_files) {
    var_name <- sub("\\.tif$", "", basename(file))
    assign(var_name, raster::raster(file))
  }
  
}
# visualize soil moisture in a soil horizon
storage_amt_feb <- function(soil_storage_rast, horizon_letter) {
  
  par(mar=c(1.8, 1.8, 1.8, 1.8))
  plot(
    soil_storage_rast, 
    main=paste0('Average Soil Moisture in the ', horizon_letter, ' Horizon for February (cm)'),
    axes=FALSE
    )
  
}

# visualize average February runoff
runoff_feb <- function(avg_runoff_rast) {
  
  par(mar=c(1.8, 1.8, 1.8, 1.8))
  plot(
    avg_runoff_rast, 
    main=paste0('Average Daily Runoff in February (cm)'),
    axes=FALSE
    )
  
}

# visualize annual precipitation
annual_precip_rast <- function(ann_precip) {
  
  par(mar=c(1.8, 1.8, 1.8, 1.8))
  plot(
    ann_precip, 
    main=paste0('Annual Precipitation (cm)'),
    axes=FALSE
  )
  
}

# visualize annual potential evapotranspiration
annual_pet_rast <- function(ann_pet) {
  
  par(mar=c(1.8, 1.8, 1.8, 1.8))
  plot(
    ann_precip, 
    main=paste0('Annual PET (cm)'),
    axes=FALSE
  )
  
}

