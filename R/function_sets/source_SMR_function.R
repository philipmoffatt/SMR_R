library(rgrass)

### File Description:
#     This function-set uses the rgrass package to run a PERL script SMR model
#     from R. It contains two helper functions and 2 main functions.

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

# Purpose: Runs the SMR model in GRASS.
# Parameters:
#   location: the location in "grassdata" to use in GRASS GIS
#   mapset: the mapset to use from the specified location
#   perl_script: the name of the perl script to use (the SMR model file)
#   map_outpur_dir: the file path specifying where maps should be written out
#   to during model runs --> just builds it if it doesn't exist alredy.
# Returns:
#   no direct returns.
# Notes:
#   This function uses get_grass_app_path() and get_grassdata_dir(). This
#   this function also assumes that the perl script sits in the "scratch"
#   though this can be changed. 
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

# Purpose: Imports all the files in a given folder into a GRASS mapset and 
#          overwrites any existing identical files in the mapset. 
# Parameters:
#   import_dir: a file path to a folder containing the maps to import
#   location (defaults NULL prompts user iput): the GRASS location to use
#   mapset (defaults NULL prompts user input): the GRASS mapset in that location
#                                              to use
# Returns:
#   No returns but there will be new maps in the GRASS data directory.
# Notes:
#   This function uses helper functions (get_grass_app_path() & 
#   grass_data_dir()) that assume certain directory structures. Check those
#   assumptions if import_files_into_GRASS() can't find locations, mapsets, or
#   the GRASS application.
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

