library(rgrass)

run_SMR <- function(location, mapset) {
  grass_dir <- Sys.glob("/Applications/*GRASS*/Contents/Resources") 
  grass_data_dir <- Sys.glob("/Users/*/grassdata") 
  
  initGRASS(gisBase = grass_dir,
            gisDbase = grass_data_dir,
            location = location,
            mapset = mapset,
            override = T)
  
  # Remove all maps from the current mapset and location
  execGRASS(cmd = "g.remove", flags = "f", type = "rast", pattern="*")
  
  smr_subdir <- "R/scratch"
  current_dir <- file.path(getwd())
  
  if (!grepl(paste0("/", smr_subdir), current_dir, fixed = TRUE)) {
    current_dir <- file.path(current_dir, smr_subdir)
  }
  
  setwd(current_dir)
  
  # Execute script
  system("perl smr_buildup.pl")
}

