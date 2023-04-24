# script to open GRASS and run smr_buildup.pl from R
library(rgrass)

# if GRASS application name is slightly different (or different version) 
# it should still grab it
grass_dir <- Sys.glob("/Applications/*GRASS*/Contents/Resources") 


# this should replace 'duncanjurayj' with 'philipmoffatt'
grass_data_dir <- Sys.glob("/Users/*/grassdata") 
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


