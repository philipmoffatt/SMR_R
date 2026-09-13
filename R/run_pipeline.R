library(yaml)

valid_stages <- c("spatial", "weather", "validate")
args <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args)) tolower(args[[1]]) else "spatial"
if (!stage %in% valid_stages) {
  stop("Usage: Rscript R/run_pipeline.R [spatial|weather|validate]")
}

config_path <- Sys.getenv("SMR_CONFIG", unset = "config.yaml")
cfg <- yaml::read_yaml(config_path)
root_dir <- Sys.getenv("SMR_ROOT_DIR", unset = cfg$root_dir)
setwd(root_dir)

file_arg <- grep("^--file=", commandArgs(), value = TRUE)
r_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
} else {
  file.path(dirname(root_dir), "R")
}

stage_scripts <- c(
  spatial = "SMR_data_setup.R",
  weather = "wxData.R",
  validate = "functional_PP.R"
)

message("=== SMR stage: ", stage, " ===")
source(file.path(r_dir, stage_scripts[[stage]]))
