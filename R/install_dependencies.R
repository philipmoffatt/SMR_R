# SMR_R Dependency Installer
# Run this once before using any other script in the project.
# Uses {pak} for fast parallel installs and unified CRAN/GitHub handling.

# ---- bootstrap pak ----
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak", repos = "https://r-lib.github.io/p/pak/stable/")
}
library(pak)

# ---- CRAN packages ----
cran_pkgs <- c(
  # config
  "yaml",
  # core data wrangling
  "tidyverse", "dplyr", "lubridate",
  # geospatial
  "sf", "terra", "raster", "mapview", "nhdplusTools",
  # federal data APIs
  "FedData",
  # DEM processing (installs R bindings; binary downloaded separately below)
  "whitebox",
  # weather / climate
  "chillR", "humidity", "bigleaf",
  # validation metrics
  "hydroGOF",
  # visualization
  "ggplot2", "RColorBrewer", "gridExtra", "cowplot"
)

pak::pkg_install(cran_pkgs)

# ---- GitHub-only packages ----
# AOI and climateR are no longer on CRAN; maintained by mikejohnson51
# rnoaa was archived from CRAN; use ropensci GitHub mirror
github_pkgs <- c(
  "mikejohnson51/AOI",
  "mikejohnson51/climateR",
  "ropensci/rnoaa"
)

pak::pkg_install(github_pkgs)

# ---- WhiteBox Tools binary ----
# The CRAN package is just R bindings; the actual binary must be downloaded once.
whitebox::install_whitebox()
whitebox::wbt_init()
cat("WhiteBox version:", whitebox::wbt_version(), "\n")

# ---- verify all packages load ----
all_pkgs <- c(
  "yaml", "tidyverse", "dplyr", "lubridate",
  "sf", "terra", "raster", "mapview", "nhdplusTools",
  "FedData", "whitebox",
  "AOI", "climateR", "rnoaa", "chillR", "humidity", "bigleaf",
  "hydroGOF",
  "ggplot2", "RColorBrewer", "gridExtra", "cowplot"
)

results <- sapply(all_pkgs, requireNamespace, quietly = TRUE)
failed  <- names(results[!results])

if (length(failed) == 0) {
  cat("\n--- All packages installed and loadable. ---\n")
} else {
  cat("\n--- FAILED to load: ---\n")
  print(failed)
}
