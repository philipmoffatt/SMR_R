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

### File Description:
#     This file provides functions that organize and manage multiple model
#     over time and functions for plotting model outputs (both tabulated and spatial),
#     debugging certain variables in the SMR model, and conducting simple evaluations
#     of the SMR model's performance against validation data.
# --------------------------------------------------------------------------- #


# Purpose: This function handles the version tracking portion of multiple SMR
#          runs. This involves creating a new folder based on the current date
#          and the ID number of the run. Additionally, this function takes in
#          in a simulation note that allows the user to record any changes to the
#          model for this run. 
# Parameters:
#   modeled_path: A character class file path to the raw modeled output data. 
#                 This will go to the raw_data/smr_outputs folder. 
#   simulation_note: A character class written note describing the changes that
#                    were made to the model. 
# Returns:
#   The file path to the new version-tracked folder that will contain modeled
#   and validation data, map outputs, and the markdown report associated with 
#   that model run.
# Notes:
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
  #map_path <- file.path(out_path, paste0("maps_", Sys.info()[["user"]])) # trying to make things user dependent so multiple runs can happen at once
  print(paste0("new user dependent map path name: ", map_path))
  
  if(!dir.exists(map_path)) {
    dir.create(map_path)
  }
  
  return(out_path)
  
}

# Purpose: Find the first and last date in a csv file of weather data. 
# Parameters: 
#   weather_data_path: The file path to the input weather data being used in the
#                      model run. 
# Returns: 
#   A vector of two date objects: the earliest date from the input weather file
#   and the latest date from the input weather file. 
# Notes: 
#   THIS FUNCTION IS NOT BEING USED!! -- In the end it didn't seem worth adding
#   this automation. 
get_run_dates <- function(weather_data_path) {
  
  weather_data <- read.csv(data_path)
  min_date <- min(weather_data$date)
  max_date <- max(weather_data$date)
  date_range <- c(as.character(min_date), as.character(max_date))
    
  return(date_range)
  
}

# Purpose: Uses the SMR model PERL file to determine which variables the model 
#          writing out to the output CSV. 
# Parameters: 
#   perl_script_path: The file path to the PERL script being used in the model run
# Returns:
#   a vector of character class objects that are the names of the variables that 
#   the SMR model is writing out. This output is then used in the preprocessing()
#   function to assign names to the columns in the modeled output CSV.
# Notes:
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

# Purpose: Reads in a CSV file and uses the first line to determine if the file
#          should be read in with sep = " " or sep = ","
# Parameters:
#   file_path: A file path to a CSV file to be read in. 
# Returns:
#   The CSV file associated with the file path. Read in with spaces or commas
#   as the separator.
# Notes:
#   This was a recent addition to this function-set but in the testing different
#   validation data sets and doing analysis in R AND Excel it gives some useful
#   flexibility to the input format of the validation data. 
read_dynamic_csv <- function(file_path) {
  # Read the first line of the file
  first_line <- readLines(file_path, n = 1)
  
  # Check if the first line contains commas
  if (grepl(",", first_line)) {
    # If commas are found, use comma as the separator
    df <- read.csv(file_path, sep = ",")
  } else {
    # If no commas are found, use space as the separator
    df <- read.csv(file_path, sep = " ")
  }
  
  return(df)
}

# Purpose: Combines the modeled output data with the validation data, adds column
#          names to make analysis possible, and moves the raw_data/smr_output
#          generic modeled output data to the processed_data/smr_output folder
#          into its correct version-tracked output folder. 
# Parameters: 
#   validation_path: A file path to the validation data. Likely in raw_data/validation_data
#   modeled_path: A file path to the modeled data. Likely in raw_data/smr_output
#   modeled_headers: A character vector of variable names in the order that they appear 
#                    in the output line of the PERL SMR model. This is produced
#                    get_print_out_line() function in this function-set. 
#   version_tracked_outpath: A file path to the new version-tracked folder in
#                            processed_data/smr_output. This will be produced with 
#                            the version_tracker() function in this function-set.
# Returns:
#   A data frame object with the modeled data and validation data combined and with 
#   column names added to the data.
# Notes: 
preprocessing <- function(
    validation_path, 
    modeled_path, 
    modeled_headers,
    version_tracked_outpath
) {
  
  #validation_data <- read.csv(validation_path) %>% 
  #  mutate(
  #    mean_discharge_cms = as.numeric(mean_discharge_cms), # Convert mean_discharge_cms to numeric
  #    date = as.Date(date, "%m/%d/%Y")
  #  )
  
  validation_data <- read_dynamic_csv(validation_path) %>%
    mutate(
      mean_discharge_cms = as.numeric(mean_discharge_cms),
      date = as.Date(date)#, "%m/%d/%Y") --> removed because the new validation data is in the yyyy-mm-dd format
    ) %>%
    fill(date, mean_discharge_cms)
  print(lapply(validation_data, class))
  validation_data <- validation_data[order(validation_data$date), ]
  

  modeled_data <- read.csv(modeled_path, sep=' ', header = FALSE, col.names = modeled_headers) %>% 
    mutate(date = as.Date(date))
  

  common_start <- max(min(validation_data$date), min(modeled_data$date))
  print(common_start)
  common_end <- min(max(validation_data$date), max(modeled_data$date))
  print(common_end)

  validation_data <- validation_data[validation_data$date >= common_start & validation_data$date <= common_end, ]
  modeled_data <- modeled_data[modeled_data$date >= common_start & modeled_data$date <= common_end, ]
  
  validation_data <- validation_data[, c('date', 'mean_discharge_cms')] %>%
    rename(validation_Q = mean_discharge_cms)
  
  modeled_data <- merge(modeled_data, validation_data)
  
  out_path <- file.path(version_tracked_outpath, 'mfc_mb.csv')
  write.csv(x = modeled_data, file = out_path)
  
  file.remove(modeled_path)
  
  return(modeled_data)
}


# Purpose: Produces a two pane plot. The upper pane compares modeled vs validation 
#          stream flow. The lower pane plots the associated precipitation, snow, and rain. 
# Parameters:
#   combined_data: The data frame returned by the preprocessing() function in this
#                  function-set. 
#   log_transform: A Boolean determining whether or not to transform the y-axis
#                  of the graph with a logarithmic function. Defaults to FALSE.
# Returns:
#   No returns. Produces a plot.
# Notes:
#   There variable names like "Q", "date", "precip_cm", "rain_cm", and "snow_cm"
#   hard coded into this function. If those ever need to be changed it is a 
#   simple fix but be aware of that. 
Q_comparison <- function(combined_data, log_transform=FALSE) {
  
  if (log_transform) {
    y_trans <- 'log10'
  } else {
    y_trans <- 'identity'
  }
  
  combined_data <- combined_data %>%
    mutate(date = as.Date(date),
           validation_Q = as.numeric(validation_Q))


  comparison_plot <- ggplot(data=combined_data) +
    geom_line(aes(x=date, y=validation_Q, color='Validation')) +
    geom_line(aes(x=date, y=Q, color='Modeled')) +
    ggtitle('Observed Versus Simulated Streamflow') +
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
    guides(color=guide_legend(title='Precip Type')) +
    scale_y_continuous(labels = scales::comma)
  
  cowplot::plot_grid(comparison_plot, precip_plot, align = "v", ncol = 1, rel_heights = c(0.60, 0.40))
}

# Purpose: Plots the main variables in the SAM component of the SMR model.
# Parameters:
#   combined_data: The data frame returned by the preprocessing() function in this
#                  function-set. 
#   log_transform: A Boolean determining whether or not to transform the y-axis
#                  of the graph with a logarithmic function. Defaults to FALSE.
# Returns:
#   No returns. Produces a plot.
# Notes:
#   There are hard coded variables in this function. If you are going to change 
#   the names of variable outputs in the PERL script know that the new names will
#   change the get_print_out_line() function names, which will change the combined_data
#   names, which will break the function. 
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

# Purpose: Calculates the Nash-Sutcliffe Efficiency (NSE) based on modeled and validation
#          stream flow using the hydroGOF package.
# Parameters:
#   combined_data: The data frame returned by the preprocessing() function in this
#                  function-set. 
# Returns:
#   A numeric class number representing the NSE value for the model run over the 
#   entire period.
# Notes:
#   Right now this doesn't remove the first year from consideration. It may make 
#   more sense to add a filter that removes the first year as those modeled values
#   are generally not as good. 
nse_Q <- function(combined_data) {
  return(NSE(combined_data$Q, combined_data$validation_Q))
}

# Purpose: Calculates the Kling-Gupta Efficiency (KGE) based on modeled and validation
#          stream flow using the hydroGOF package.
# Parameters:
#   combined_data: The data frame returned by the preprocessing() function in this
#                  function-set. 
# Returns: 
#   A numeric class number representing the KGE value for the model run over the 
#   entire period.
# Notes:
#   Right now this doesn't remove the first year from consideration. It may make 
#   more sense to add a filter that removes the first year as those modeled values
#   are generally not as good. 
kge_Q <- function(combined_data) {
  return(KGE(combined_data$Q, combined_data$validation_Q))
}

# Purpose: Produces a single plot of a single flux over the entire modeled period
#          in the combined_data data frame. 
# Parameters:
#   combined_data: The data frame returned by the preprocessing() function in this
#                  function-set. 
#   flux: A character class name of the variable to plot.
#   log_transform: A Boolean determining whether or not to transform the y-axis
#                  of the graph with a logarithmic function. Defaults to FALSE.
# Returns:
#   No returns. Produces a plot.
# Notes:
#   You can use this function in isolation but it was designed to be used as a helper 
#   function in the flux_ts_loop() function. This allows multiple variables
#   to be plotted individually  without needing to call flux_ts() many times.
flux_ts <- function(combined_data, flux, log_transform = FALSE) {
  
  if (log_transform) {
    y_trans <- 'log10'
  } else {
    y_trans <- 'identity'
  }
  
  flux_series <- combined_data[, flux] / 72678
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

# Purpose: Loops over the flux_ts() function for a vector of variables to be plotted.
# Parameters:
#   log_transform: A Boolean determining whether or not to transform the y-axis
#                  of the graph with a logarithmic function. Defaults to FALSE.
# Returns: 
#   combined_data: The data frame returned by the preprocessing() function in this
#                  function-set.
#   fluxes: A character vector of variable names for plotting
#   log_transform: A Boolean determining whether or not to transform the y-axis
#                  of the graph with a logarithmic function. Defaults to FALSE.
# Notes:
flux_ts_loop <- function(combined_data, fluxes, log_transform=FALSE) {
  
  for (flux in fluxes) {
    print(flux_ts(combined_data, flux, log_transform = log_transform))
  }
  
}

# Purpose: Produces a bar plot of annual accumulated fluxes for a set of variables
#          recorded by the model. 
# Parameters
#   combined_data: The data frame returned by the preprocessing() function in this
#                  function-set.
#   fluxes: A character vector of the variables to include in the plot.
#   log_transform: A Boolean determining whether or not to transform the y-axis
#                  of the graph with a logarithmic function. Defaults to FALSE.
# Returns:
#   No returns. Produces a plot.
# Notes: 
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
    ggtitle("Annual Accumulated Fluxes") +
    theme(plot.title = element_text(hjust = 0.5))

}

# Purpose: Produces a plot of the variables relevant to the radiation balance of the 
#          surface: srad, latent, sensible, long wave, and q_rain_ground
# Parameters:
#   combined_data: The data frame returned by the preprocessing() function in this
#                  function-set. 
#   log_transform: A Boolean determining whether or not to transform the y-axis
#                  of the graph with a logarithmic function. Defaults to FALSE.
# Returns:
#   No returns. Produces a plot.
# Notes:
#   There are hard coded variables in this function. If you are going to change 
#   the names of variable outputs in the PERL script know that the new names will
#   change the get_print_out_line() function names, which will change the combined_data
#   names, which will break the function. 
radiation_ts <- function(combined_data, log_transform = FALSE) {
  
  # Set log transform for y-axis scale
  if (log_transform) {
    y_trans <- 'log10'
  } else {
    y_trans <- 'identity'
  }
  
  # Divide all values in selected columns by 72678
  combined_data <- combined_data %>%
    mutate(
      srad = srad/72678,
      latent = latent/72678,
      sensible = sensible/72678,
      lw = lw/72678,
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

# Purpose: Produces a plot of the mass balance of the model. 
# Parameters:
#   combined_data: The data frame returned by the preprocessing() function in this
#                  function-set. 
#   log_transform: A Boolean determining whether or not to transform the y-axis
#                  of the graph with a logarithmic function. Defaults to FALSE.
# Returns:
#   No returns. Produces a plot.
# Notes:
mass_balance <- function(combined_data, log_transform = FALSE) {
  
  # Set log transform for y-axis scale
  if (log_transform) {
    y_trans <- 'log10'
  } else {
    y_trans <- 'identity'
  }
  
  annual_combined <- combined_data #%>%
    #group_by(year) #%>%
    #select(-c(date, wshed_id)) %>% # seeing if I need this
    #summarize_all(.funs=sum)
  
  view(annual_combined)
  annual_combined$mass_balance <-
    annual_combined$precip_cm +
    annual_combined$actual_ET_daily -
    annual_combined$runoff_cm -
    annual_combined$perc_cm
  
  ggplot(data=annual_combined) +
    scale_y_continuous(trans=y_trans) +
    geom_line(aes(x=date, y=mass_balance)) +
    ggtitle('Mass Balance') +
    theme(plot.title = element_text(hjust = 0.5))
}

# Purpose: Produces a plot of the variables that make up the snow water equivalent 
#          variable. This was an especially sensitive variable so being to study its
#          components each time can be useful. 
# Parameters:
#   combined_data: The data frame returned by the preprocessing() function in this
#                  function-set. 
#   log_transform: A Boolean determining whether or not to transform the y-axis
#                  of the graph with a logarithmic function. Defaults to FALSE.
# Returns:
#   No returns. Produces a plot. 
# Notes:
#   There are hard coded variables in this function. If you are going to change 
#   the names of variable outputs in the PERL script know that the new names will
#   change the get_print_out_line() function names, which will change the combined_data
#   names, which will break the function. 
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

# Purpose: Produces a plot of the two variables that directly make up the latent 
#          heat component of the SMR-SAM surface radiation budget. 
# Parameters:
#   combined_data: The data frame returned by the preprocessing() function in this
#                  function-set. 
#   log_transform: A Boolean determining whether or not to transform the y-axis
#                  of the graph with a logarithmic function. Defaults to FALSE.
# Returns:
#   No returns. Produces a plot.
# Notes:
#   There are hard coded variables in this function. If you are going to change 
#   the names of variable outputs in the PERL script know that the new names will
#   change the get_print_out_line() function names, which will change the combined_data
#   names, which will break the function. 
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

# Purpose: Produces a simple plot of a single raster file.
# Parameters:
#   map: A raster() object or really any spatial object.
#   title: A character for the title of the plot.
# Returns:
#   No returns. Produces a plot.
# Notes:
#   You can use this function individually, it is used in the modeling workflow as 
#   helper function to the get_map_outputs() function which loops through a folder
#   of maps and calls this function on each map.
visualize_map_outputs <- function(map, title) {

  par(mar=c(1.8, 1.8, 1.8, 1.8))
  plot(
    map, 
    main=paste0(title),
    axes=FALSE
  )
  
}

# Purpose: Produces plots for all the raster files in a folder.
# Parameters:
#   raw_map_dir: A file path to the folder containing all the maps from the most
#                recent model run. Likely this will be something like: raw_data/smr_output/maps_[user]
# Returns:
#   No returns. Produces multiple plots.
# Notes:
#   This function names the map plots just based on the name of the file it finds. 
#   this means that if you want different plot titles you need to change the output
#   map names in the PERL SMR script itself.
get_map_outputs <- function(raw_map_dir, version_tracked_map_dir) {
  
  print(paste0("version tracked map path: ", version_tracked_map_dir))
  
  tif_files <- list.files(path = raw_map_dir, pattern = "\\.tif$", full.names = TRUE)
  
  if (length(tif_files) == 0) {
    cat("no .tif files found in raw_map_dir.\n")
    return()
  }
  
  for (tif_file in tif_files) {
    success <- file.copy(from = tif_file, to = file.path(version_tracked_map_dir, basename(tif_file)))
    if (!success) {
      cat("failed to copy", tif_file, "to", version_tracked_map_dir, "\n")
    } else {
      var_name <- sub("\\.tif$", "", basename(tif_file))
      assign(var_name, raster::raster(file.path(version_tracked_map_dir, basename(tif_file))))
    }
  }
  
  #generic_version_tracked_map_dir <- gsub("/maps_duncanjurayj$", "/maps", version_tracked_map_dir)
  
  #file.rename(version_tracked_map_dir, generic_version_tracked_map_dir)
  
  for (var_name in ls()) {
    if (class(get(var_name)) == "RasterLayer") {
      title <- paste0(var_name, " (cm)")
      visualize_map_outputs(get(var_name), title)
    }
  }
}






