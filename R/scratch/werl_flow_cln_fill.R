# cleaning WERL flow data for model validation/check
library(tidyverse)
library(dplyr)
library(lubridate)

# read raw mfc flow data and find average daily discharge for recorded days
mfc_raw <- read.csv('raw_data/validation_data/WERL_flow_data_2018-2022.csv')
mfc_raw$datetime <- strptime(mfc_raw$date_time,format='%m/%d/%y %H:%M')
mfc_cln <- mfc_raw %>%
    mutate(date = as.Date(datetime),   # change date format
           day = day(date), 
           month = month(date), 
           year = year(date),
           cms = as.numeric(cms)) %>% # change data type for cms
  drop_na() %>%                       # drop na rows
  group_by(date) %>%                  # group the 10min measurements by day
  summarise(avgCMS = round(mean(cms),3))   # find average daily cms
plot(mfc_cln$date, mfc_cln$avgCMS)        # plot data
prgerr <- mfc_cln$date[mfc_cln$avgCMS>10]           # dates where flow device programing was wrong
mfc_cln <- data.frame(subset(mfc_cln, !(date %in% prgerr)))     # drop out the error date from Oct 2019

write.csv(mfc_cln,'raw_data/validation_data/WERL_flow_data_2018-2022_cln.csv')

# read in usgs data for Paradise Creek
pc_raw <- read.csv('raw_data/validation_data/usgs_pc_2018-2021.csv')
pc_cln <- pc_raw %>% 
  mutate(date = as.Date(date_pc, format = '%m/%d/%y'), # change date format
         pc_cms = round(pc_cfs * 0.028316847,3))          # convert cfs to cms

# find common record dates between mfc and pc
common_dates <- as.Date(intersect(mfc_cln$date, pc_cln$date), origin = '1970-01-01')
pc_cln_subset <- subset(pc_cln, date %in% common_dates)  # subset the pc data to common date records         
mfc_fill <- inner_join(mfc_cln, pc_cln_subset) %>%        # join common records
  dplyr::select(date,avgCMS,pc_cms)
plot(mfc_fill)                                          # plot check on flow data

# find a relationship between pc and mfc
fill_mdl <- lm(mfc_fill$avgCMS~mfc_fill$pc_cms)
summary(fill_mdl)

# find dates with no flow records from mfc
missing_records <- as.Date(setdiff(pc_cln$date, mfc_cln$date), origin = '1970-01-01')

# calculate an estimate for the flow in mfc based on pc data
mfc_mis_records <- pc_cln %>% 
  mutate(mfc_fill_cms = 0.7753186*pc_cln$pc_cms+0.0670547) %>% 
  dplyr::select(date, mfc_fill_cms)

# fill the empty records for mfc with the calculated values 
mfc_filled <- full_join(mfc_mis_records, mfc_cln) %>% 
  mutate(mean_discharge_cms = coalesce(avgCMS,mfc_fill_cms))

plot(mfc_filled) # check 

write.csv(mfc_filled,'raw_data/validation_data/mfc_cln_filled.csv')





