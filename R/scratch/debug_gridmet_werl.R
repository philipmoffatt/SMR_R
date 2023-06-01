library(tidyverse)
library(dplyr)
library(lubridate)
library(hydroGOF)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(cowplot)
library(raster)
library(mapview)
source("./R/scratch/functional_PP.R")

werl_data <- read.csv("/Users/duncanjurayj/Dropbox/SMR_R/raw_data/validation_data/WERLdischarge.csv")
usgs_data <- read.csv("/Users/duncanjurayj/Dropbox/SMR_R/raw_data/validation_data/USGSdischarge.csv")

gridMET_mini <- read.csv("./raw_data/weather/gridMET_mini.csv", sep=' ', header = FALSE, col.names = column_headers) %>% 
  mutate(date = as.Date(date))
colnames(gridMET_mini)

noaa_mini <- read.csv("./raw_data/weather/noaa_pullman_mini.csv", sep=' ', header = FALSE, col.names = column_headers) %>% 
  mutate(date = as.Date(date))
colnames(gridMET_mini)