library(rmarkdown)
library(dplyr)

### File Description:
#     This file contains functions that are designed to facilitate the 
#     knitting of a markdown report which in turn triggers the running of the 
#     SMR model using the source_SMR_function.R function-set file.
# --------------------------------------------------------------------------- #


# Purpose: 
# Parameters:
# Returns:
# Notes:
get_next_output_path <- function(processed_path="processed_data/smr_output") {
  
  outputDir <- processed_path
  
  SMR_log_file <- file.path(outputDir, "SMR_run_log.txt")
  print(paste0("SMR Log File Path: ", SMR_log_file))
  
  if(!file.exists(SMR_log_file)){
    runID <- "1"
  } else {
    table_temp <- (read.table(file=SMR_log_file) %>% row.names())
    runID <- table_temp[length(table_temp)] %>% as.integer() + 1
  }  
  
  runID <- paste0("MFC_", rep(0,(3-nchar(runID))) %>% paste0(., collapse = ""), runID)
  
  out_path <- file.path(outputDir, paste0(runID, "_", Sys.Date()))
  
  print(paste0("This is the final outpath: ", out_path))
  
  return(out_path)
  
}

# Purpose: Find the path to the grassdata folder that contains the locations
#          and the mapsets for GRASS GIS.
# Parameters:
#   no parameters
# Returns:
#   The full file path to the "grassdata" folder containing the maps for 
#   GRASS GIS.
# Notes:
#   This function is designed to work for windows, max, and linux, but for
#   the windows case it is hard coded to Philip Moffat's windows path. The 
#   primary assumption of this function is that maps for GRASS GIS are held in 
#   a folder called "grassdata"
get_grassdata_dir <- function() {
  os <- Sys.info()["sysname"]
  grass_data_glob <- "grassdata"
  
  if (os == "Linux") {
    grass_data_dir <- Sys.glob(file.path("/home", "*", "Documents", grass_data_glob))
  } else if (os == "Darwin") {
    grass_data_dir <- Sys.glob(file.path("/Users", "*", grass_data_glob))
  } else if (os == "Windows") {
    grass_data_dir <- "C:/Users/philip.moffat/Documents/grassdata" # hard coded solution
  }
  
  return(grass_data_dir)
}

# Purpose: Find the full path to the GRASS GIS application and navigates to the
#          "Resources" folder within it.
# Parameters:
#   no parameters
# Returns:
#   The full file path to the GRASS application Resources folder.
# Notes:
#   This function is designed to work for windows, mac, and linux operating 
#   systems. That being said, the windows case is currently hard coded and 
#   assume GRASS GIS 8.2. For Linux and Darwin it is more flexible. 
get_grass_app_path <- function() {
  os <- Sys.info()["sysname"]
  
  # currently only works with linux and mac
  if (os == "Linux") {
    grass_path <- Sys.glob("program files/*grass*/Contents/Resources") # file.path("program files", "grass 8.2", "Contents", "Resources")
  } else if (os == "Darwin") {
    grass_path <- Sys.glob("/Applications/*GRASS*/Contents/Resources") # using a * instead of applications might be better --> works the same
  } else if (os == "Windows") {
    grass_path <- "C:/Program Files/GRASS GIS 8.2" # hard coded solution
  }
  
  return(grass_path)
}

# Purpose: Knits a markdown report and passes parameters to the report that 
#          facilitate the running of a specified SMR PERL model in a specified 
#          Location and Mapset inside GRASS GIS. 
# Parameters:
#   location: A character class name of the folder to use for location in the GRASS GIS database. Defaults to NULL.
#   mapset: A character class name of the folder to use for the mapset in the GRASS GIS database. Defaults to NULL.
#   perl_script: A character class name of the perl script to run SMR from in SMR_R/R/function_sets. Defaults to NULL.
#   simulation_note: A character class sentence describing the changes to the model for the current run. Defaults to NULL.
#   markdown_generic_outpath: A file path to the processed_data/smr_output location. Defaults to "processed_data/smr_output".
#   run_start_date: A character class yyyy-mm-dd object specifying the first day of the input data.
#   run_end_date: A character class yyyy-mm-dd object specifying the last day of the input data.
# Returns:
#   No returns. This function knits an R Markdown report resulting in an SMR 
#   model run. 
# Notes:
#   You can have the function written out with parameters set-up statically
#   if you think you will only be changing the internal setup of the model in PERL
#   for a while. Additionally, the function can be called with no parameters and 
#   it will prompt you for everything it needs. 
smr_wrapper <- 
  function(location = NULL, 
           mapset = NULL, 
           perl_script = NULL,
           simulation_note = NULL,
           markdown_generic_outpath = "processed_data/smr_output",
           run_start_date = NULL,
           run_end_date = NULL
           ) {
  
  if (is.null(run_start_date) || is.null(run_end_date)) {
    run_start_date <- readline(prompt = "Input a start date for the model run period: ")
    run_end_date <- readline(prompt = "Input an end date for the model run period: ")
  }
    
  # getting to the GRASS application
  grass_dir <- get_grass_app_path() # Sys.glob("/Applications/*GRASS*/Contents/Resources") 
  
  # find grassdata directory to prompt user for location and mapset
  grass_data_dir <- get_grassdata_dir() # Sys.glob("/Users/*/grassdata")
  
  
  if (is.null(location) || is.null(mapset) || is.null(perl_script)) {
    locations <- dir(grass_data_dir)
    location_list <- paste0(seq_along(locations), ". ", locations)
    
    full_string <- paste0("Locations in ", grass_data_dir, ":\n", paste(location_list, collapse = "\n"))
    cat(full_string)
    
    # prompt user for location and mapsets
    location <- readline(prompt = "Choose a location from the list of locations above: ")
    
    # make a list of the mapsets excluding PERMANENT in the location
    mapsets <- setdiff(dir(paste0(grass_data_dir, "/", location)), "PERMANENT")
    
    # print out the mapsets in that locations
    cat(paste0(cat("Mapsets in", location, ":\n "), seq_along(mapsets), ". ", mapsets, "\n"))
    
    # prompt user for a mapset name
    mapset <- readline(prompt = "Choose a mapset from the list of mapsets above: ")
    
    # prompt user for perl script name
    perl_script <- readline(prompt = "Write the name of the perl script you want to run: ")
    
  }
  
  if (is.null(simulation_note)) {
    simulation_note <- readline(prompt = "Write a simulation note describing changes to the model for this run: ")
  }
  
  params <- list(location = location, 
                 mapset = mapset,
                 perl_script_name = perl_script, # this is now an input variable with a user prompt
                 simulation_note = simulation_note,
                 date = Sys.Date(),
                 run_start_date = run_start_date,
                 run_end_date = run_end_date)
  
  print(paste0("markdown generic outpath: ", markdown_generic_outpath))
  output_path <- get_next_output_path(processed_path = markdown_generic_outpath)
  print(paste0("output_path: ", output_path))
  
  report_folder <- "R/function_sets" # could be made more flexible
  
  report_path <- file.path(report_folder, "model_report_v2.Rmd") # changed this but it should be made a parameter
  print(paste0("report path: ", report_path))
  
  markdown_report <- render(input = report_path,
         output_format = "pdf_document",
         params = params,
         output_dir = file.path(output_path, "model_report"))

}



