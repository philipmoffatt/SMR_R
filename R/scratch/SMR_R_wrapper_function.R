# this function will take readline inputs for location and mapset and 
# pass them as parameters to an rmarkdown which will then run the model, produce
# the outputs, and write a pdf to the correct directory
library(rmarkdown)
library(dplyr)

get_next_output_path <- function(processed_path="processed_data/smr_output") {
  
  outputDir <- processed_path
  
  SMR_log_file <- file.path(outputDir, "SMR_run_log.txt")
  print(paste0("SMR Log File Path: ", SMR_log_file))
  
  if(!file.exists(SMR_log_file)){
    runID <- "1"
  } else {
    table_temp <- (read.table(file= SMR_log_file) %>% row.names())
    runID <- table_temp[length(table_temp)] %>% as.integer() + 1
  }  
  
  runID <- paste0("MFC_", rep(0,(3-nchar(runID))) %>% paste0(., collapse = ""), runID)
  
  out_path <- file.path(outputDir, paste0(runID, "_", Sys.Date()))
  
  print(paste0("This is the final outpath: ", out_path))
  
  return(out_path)
  
}

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
  
  report_folder <- "R/scratch" # could be made more flexible
  
  report_path <- file.path(report_folder, "model_report_v2.Rmd") # changed this but it should be made a parameter
  print(paste0("report path: ", report_path))
  
  markdown_report <- render(input = report_path,
         output_format = "pdf_document",
         params = params,
         output_dir = file.path(output_path, "model_report"))

}



