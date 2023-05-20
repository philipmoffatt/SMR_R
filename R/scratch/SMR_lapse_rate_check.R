# script to check that the lapse rate equations in SMR will behave as expected 
# across elevation gradients

library(ggplot2)

### Relevant SMR perl script code
  # lapse rate variables
  # $high_site_el = 1443;
  # $low_site_el = 784;
  # $mean_annual_onsite_precip = 43.921; # cm 
  # $mean_annual_sea_level_precip = 15.53944; # calculated based on pullman average and lapse rate
  # $precip_LR = 0.09195; # cm / m
  # $temp_LR = -3.20; # C / km
  
  # this equation may be wrong 
  # print `r.mapcalc 'precip = (el*$precip_LR+$mean_annual_sea_level_precip) / ($mean_annual_onsite_precip) * $precip' --o`;
  # print `r.mapcalc 'tavg = $tavg + ($temp_LR/1000) * (el - $low_site_el)' --o`;
  # print `r.mapcalc 'tdew = $tdew + ($temp_LR/1000) * (el - $low_site_el)' --o`;
  # print `r.mapcalc 'rain = if(tavg>tmax_rain,precip,if(tavg<tmin_snow,0.0,(tavg-tmin_snow)/(tmax_rain-tmin_snow)*precip))' --o`;
  # print `r.mapcalc 'snow = precip-rain' --o`;
  
  # THESE ARE STILL FROM MICA!!! NEED TO UPDATE -- Duncan 3/12/23 --> these have been updated to get MM from Pullman: Duncan 3/13/23
  # print `r.mapcalc 'pet_data = (1.8281 * $pet - 0.2486)' --o`;
  # print `r.mapcalc "pet = (pet_data - $pet) * (el - $low_site_el)/($high_site_el - $low_site_el)+pet_data" --o`;

# foreach $hrly_tmp (@hourly_tmp_array)
# {
  #			print `r.mapcalc 'temp = $hrly_tmp' --o`; # MZ 20170210
#  print `r.mapcalc 'temp = $hrly_tmp + ($temp_LR/1000)*($el - $low_site_el)' --o`;  # switched order of el and low_site_el as is in the initial lapse rate calculations above
#  print `r.mapcalc 'temp_sum = max(0.0,temp)+temp_sum' --o`;
  
# };
# print "\n \n";
# print "\n|----- 6 HOUR TEMPERATURE ARRAY -----|\n";
# print "\n \n";
# foreach $hrly_tmp (@hourly_tmp_array)
# {
  #			print `r.mapcalc 'temp = $hrly_tmp' --o`; # MZ 20170210
#  print `r.mapcalc 'temp = $hrly_tmp + ($temp_LR/1000)*(el - $low_site_el)' --o`; # switched order of el and low_site_el as is in the initial lapse rate calculations above
  
  #print "$hrly_tmp \n";

### Expectations
# 1. temperature (tavg and tdew) is expected to decrease with increasing elevation -- result: behaved as expected
# 2. precipitation is expected to increase with increasing elevation -- result: behaved as expected but should review intercept with philip
# 3. pet is expected to decrease with increasing elevation -- 

### Setting up constants from perl SMR in R
elevation_test_gradient <- seq(1,2000, 1) # meters
### Relevant SMR perl script code
# lapse rate variables
high_site_el = 1443;
low_site_el = 784;
mean_annual_onsite_precip = 43.921; # cm 
mean_annual_sea_level_precip = 15.53944; # calculated based on pullman average and lapse rate
precip_LR = 0.09195; # cm / m
temp_LR = -3.20; # C / km

### Testing temperature lapse rate behavior based on SMR equation
input_tavg <- 10 # celsius - the assumption here is that this is the temp at 784 (because that's the low_site_el where measurements are taken from)
input_tdew <- 8 # celsius - the assumption here is that this is the temp at 784 (because that's the low_site_el where measurements are taken from)

# tavg equation
tavg_gradient <- input_tavg + (temp_LR/1000) * (elevation_test_gradient - low_site_el)

# plot tavg behavior
ggplot() + 
  geom_point(aes(y=tavg_gradient, x=elevation_test_gradient)) + 
  geom_point(aes(y=input_tavg, x=low_site_el), size=5, color='red') + 
  geom_point(aes(y=tavg_gradient[1433], x=high_site_el), size=5, color='blue') +
  geom_text(aes(y = input_tavg, x = low_site_el + 300, label = "Temperature Recorded at Pullman"), color = "black", size=3) +
  geom_text(aes(y = tavg_gradient[high_site_el], x = high_site_el - 370, label = "Temperature Recorded at Moscow Mountain"), color = "black", size=3) +
  labs(title="Temperature as Function of Elevation based on an Empirical Lapse Rate",
       y="Average Temperature (C)",
       x="Elevation (m)") +
  theme(plot.title = element_text(hjust = 0.5))

# tdew equation
tdew_gradient <- input_tdew + (temp_LR/1000) * (elevation_test_gradient - low_site_el)

# plot tdew behavior
ggplot() + 
  geom_point(aes(y=tdew_gradient, x=elevation_test_gradient)) + 
  geom_point(aes(y=input_tdew, x=low_site_el), size=5, color='red') + 
  geom_point(aes(y=tdew_gradient[1433], x=high_site_el), size=5, color='blue') +
  geom_text(aes(y = input_tdew, x = low_site_el + 300, label = "Dewpoint Temperature Predicted at Pullman"), color = "black", size=3) +
  geom_text(aes(y = tdew_gradient[high_site_el], x = high_site_el - 370, label = "Dewpoint Temperature Predicted at Moscow Mountain"), color = "black", size=3) +
  labs(title="Dewpoint Temperature as Function of Elevation based on an Empirical Lapse Rate",
       y="Dewpoint Temperature (C)",
       x="Elevation (m)") +
  theme(plot.title = element_text(hjust = 0.5))

### Testing precipitation lapse rate behavior based on SMR equation
input_precip <- 2 # centimeters - this would be a precipitation value at Pullman (so at 784 meters)
precip_gradient <-  ((elevation_test_gradient*precip_LR+mean_annual_sea_level_precip) / (mean_annual_onsite_precip)) * input_precip

ggplot() + 
  geom_point(aes(y=precip_gradient, x=elevation_test_gradient)) + 
  geom_point(aes(y=input_precip, x=low_site_el), size=5, color='red') + 
  geom_point(aes(y=precip_gradient[1433], x=high_site_el), size=5, color='blue') +
  geom_text(aes(y = input_precip, x = low_site_el + 300, label = "Precipitation Recorded at Pullman"), color = "black", size=3) +
  geom_text(aes(y = precip_gradient[high_site_el], x = high_site_el - 370, label = "Precipitation Predicted at Moscow Mountain"), color = "black", size=3) +
  labs(title="Precipitation as Function of Elevation based on an Empirical Lapse Rate",
       y="Precipitation (cm)",
       x="Elevation (m)") +
  theme(plot.title = element_text(hjust = 0.5))


### Testing PET lapse rate behavior based on SMR equation
# note: pet_data is the high pet estimate based on a linear regression

input_pet <- 2 # centimeters

### ALTERNATIVES? ###

pet_lapse_rate <- -2.226955e-05 # centimeters per meter
# trying pet lapse rate in same way as temperature lapse rate
alternative_pet_gradient <- input_pet + ((pet_lapse_rate) * (elevation_test_gradient - low_site_el))

# plot of current pet gradient method

# plot alternative pet gradient with identical method to temperature
ggplot() +
  geom_point(aes(y=alternative_pet_gradient, x=elevation_test_gradient)) + 
  geom_point(aes(y=input_pet, x=low_site_el), size=5, color='red') + 
  geom_point(aes(y=alternative_pet_gradient[high_site_el], x=high_site_el), size=5, color='blue') +
  geom_text(aes(y = input_pet, x = low_site_el + 200, label = "PET at Pullman"), color = "black", size=3) +
  geom_text(aes(y = alternative_pet_gradient[high_site_el], x = high_site_el - 400, label = "PET Predicted at Moscow Mountain (alternative)"), color = "black", size=3) +
  labs(title="PET as Function of Elevation based on an Empirical Lapse Rate (alternative)",
       y="PET (cm)",
       x="Elevation (m)") +
  theme(plot.title = element_text(hjust = 0.5))

### Plot of alternative precip gradient method with simple 'temperature-like' lapse rate calculation
alternative_precip_LR <- 0.0002221415 # centimeters per meter
alternative_precip_gradient <- input_precip + ((alternative_precip_LR) * (elevation_test_gradient - low_site_el))

ggplot() +
  geom_point(aes(y=alternative_precip_gradient, x=elevation_test_gradient)) + 
  geom_point(aes(y=input_precip, x=low_site_el), size=5, color='red') + 
  geom_point(aes(y=alternative_precip_gradient[high_site_el], x=high_site_el), size=5, color='blue') +
  geom_text(aes(y = input_precip, x = low_site_el + 250, label = "Precipitation at Pullman"), color = "black", size=3) +
  geom_text(aes(y = alternative_precip_gradient[high_site_el], x = high_site_el - 570, label = "Precipitation Predicted at Moscow Mountain (alternative)"), color = "black", size=3) +
  labs(title="Precipitation as Function of Elevation based on an Empirical Lapse Rate (alternative)",
       y="Precipitation (cm)",
       x="Elevation (m)") +
  theme(plot.title = element_text(hjust = 0.5))

