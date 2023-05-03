library(rgrass)

run_SMR <- 
  function(location = NULL, mapset = NULL) {
    
  grass_dir <- Sys.glob("/Applications/*GRASS*/Contents/Resources") 
  grass_data_dir <- Sys.glob("/Users/*/grassdata")
  
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
  
  execGRASS(cmd = "g.remove", flags = "f", type = "rast", pattern="*")
  print("grass maps removed")
  smr_subdir <- "R/scratch"
  current_dir <- file.path(getwd())
  
  if (!grepl(paste0("/", smr_subdir), current_dir, fixed = TRUE)) {
    current_dir <- file.path(current_dir, smr_subdir)
  }
  
  setwd(current_dir)
  
  system("perl smr_buildup.pl")
}

import_files_into_GRASS <- 
  function(import_dir, location = NULL, mapset = NULL) {
  
  grass_dir <- Sys.glob("/Applications/*GRASS*/Contents/Resources") 
  grass_data_dir <- Sys.glob("/Users/*/grassdata") 
  
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
  
  csv_ini_files <- grep("\\.(csv|ini)$", files_in_dir)
  files_in_dir <- files_in_dir[-csv_ini_files]
  
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


