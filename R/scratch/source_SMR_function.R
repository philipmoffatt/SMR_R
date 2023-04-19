library(rgrass)

run_SMR <- function(location, mapset) {
  grass_dir <- Sys.glob("/Applications/*GRASS*/Contents/Resources") 
  grass_data_dir <- Sys.glob("/Users/*/grassdata") 
  
  initGRASS(gisBase = grass_dir,
            gisDbase = grass_data_dir,
            location = location,
            mapset = mapset,
            override = T)
  
  execGRASS(cmd = "g.remove", flags = "f", type = "rast", pattern="*")
  
  smr_subdir <- "R/scratch"
  current_dir <- file.path(getwd())
  
  if (!grepl(paste0("/", smr_subdir), current_dir, fixed = TRUE)) {
    current_dir <- file.path(current_dir, smr_subdir)
  }
  
  setwd(current_dir)
  
  system("perl smr_buildup.pl")
}


import_files_into_GRASS <- function(import_dir, location_name, mapset_name) {
  
  grass_dir <- Sys.glob("/Applications/*GRASS*/Contents/Resources") 
  grass_data_dir <- Sys.glob("/Users/*/grassdata") 
  
  initGRASS(gisBase = grass_dir,
            gisDbase = grass_data_dir,
            location = location_name,
            mapset = mapset_name,
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
              flags="o")
    
    print(paste0(file, " was imported"))
    
  }
}


