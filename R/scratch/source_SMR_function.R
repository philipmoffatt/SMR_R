library(rgrass)

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

run_SMR <- 
  function(location = NULL, mapset = NULL, perl_script, 
           map_output_dir = NULL) {
    
    if (!is.null(map_output_dir)) {
      if (!file.exists(map_output_dir)) {
        dir.create(map_output_dir, recursive = TRUE)
      }
    }
    
  grass_dir <- get_grass_app_path() # Sys.glob("/Applications/*GRASS*/Contents/Resources") 
  #grass_data_dir <- Sys.glob("/Users/*/grassdata")
  grass_data_dir <- get_grassdata_dir() # testing operating system flexible grass data directory
  
  if (is.null(location) || is.null(mapset)) {
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
    
    perl_script <- readline(prompt = "Write the name of the perl script you want to run: ")

    
  }
  
  initGRASS(gisBase = grass_dir,
            gisDbase = grass_data_dir,
            location = location,
            mapset = mapset,
            override = T)
  
  execGRASS(cmd = "g.remove", flags = "f", type = "rast", pattern="*")
  print("grass maps removed")
  smr_subdir <- "R/scratch"
  current_dir <- file.path(getwd())
  
  if (!grepl(paste0("/", smr_subdir), current_dir, fixed = TRUE)) {
    current_dir <- file.path(current_dir, smr_subdir)
  }
  
  setwd(current_dir)
  
  #system("perl smr_buildup.pl") 
  
  system(paste0("perl ", perl_script))
  
}

import_files_into_GRASS <- 
  function(import_dir, location = NULL, mapset = NULL) {
  
  grass_dir <- get_grass_app_path()  # Sys.glob("/Applications/*GRASS*/Contents/Resources") 
  grass_data_dir <- get_grassdata_dir() # Sys.glob("/Users/*/grassdata") 
  
  if (is.null(location) || is.null(mapset)) {
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
    
  }

  initGRASS(gisBase = grass_dir,
            gisDbase = grass_data_dir,
            location = location,
            mapset = mapset,
            override = T)
  
  files_in_dir <- list.files(import_dir)

  #csv_ini_files <- grep("\\.(csv|ini)$", files_in_dir)
  #files_in_dir <- files_in_dir[-csv_ini_files]
  print("Files in Dir:")
  print(files_in_dir)
  
  # import each file into GRASS
  for (file in files_in_dir) {
    
    file_dir <- file.path(import_dir, file)
    
    execGRASS("r.import", 
              input = file_dir, 
              output = gsub("\\.[^.]+$", "", file), 
              flags=c("o", "overwrite"))
    
    print(paste0(file, " was imported"))
    
  }
}

