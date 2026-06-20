library(yaml)

cfg      <- yaml::read_yaml("config.yaml")
root_dir <- cfg$root_dir
setwd(root_dir)

r_dir    <- file.path(dirname(root_dir), "R")   # ../R relative to scratch/

message("=== Stage 1: Spatial preprocessing ===")
source(file.path(r_dir, "SMR_data_setup.R"))
# STOP HERE after Stage 1 and verify config.yaml:basin_mask_value matches
# the basin ID shown in the dem_basins plot, then re-run from Stage 2.

message("=== Stage 2: Weather file ===")
source(file.path(r_dir, "wxData.R"))

message("=== Stage 3: Simulation ===")
message("Run from PowerShell after setting OSGeo4W env vars (see SETUP.md):")
message("  cd ", root_dir)
grass_py  <- "C:/OSGeo4W/bin/python3.exe"
grass_sc  <- "C:/OSGeo4W/apps/grass/grass85/etc/grass85.py"
location  <- file.path(root_dir, "grass_workspace/BullRun/PERMANENT")
wx_file   <- file.path(root_dir, "raw_data/weather/gridmet_weather.csv")
out_file  <- file.path(root_dir, "smr_output/MFC_mass_balance_1.csv")
perl_script <- file.path(r_dir, "smr_buildup.pl")
message(sprintf('  & "%s" "%s" --text "%s" --exec perl "%s" "%s" "%s"',
                grass_py, grass_sc, location, perl_script, wx_file, out_file))

message("=== Stage 4: Validation ===")
source(file.path(r_dir, "functional_PP.R"))
