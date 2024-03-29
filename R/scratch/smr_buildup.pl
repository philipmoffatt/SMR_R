#!/usr/bin/perl

print `g.remove -f type=raster name=MASK`;
print `g.region rast=watershed`;
print `r.mapcalc 'watershed = if(watershed==1.0,1,1)' --o`;
#print `rm Q_* `;
#print `rm M_* `;
#print `rm R_* `;
#print `rm tmp.txt`;
#print `rm runoff_total_* `;
#print `rm saturation_* `; 
#print `rm Psat_* `;

# Initialize variables
$gridsize = 30.0; # m

# The time step variable refers to the time step at which snowmelt is calculated
$time_step = 24.0; # hrs

# The temp time step determines how frequent the hydrology is simulated.  This requires
# temperature data in the weather file.  For example a 1 hr $temp_time_step requires hourly temperature data in the weather file
$temp_time_step = 6.0; # hrs

# lapse rate variables
$high_site_el = 1443;
$low_site_el = 784;
$mean_annual_onsite_precip = 43.921; # cm 
$mean_annual_sea_level_precip = 15.53944; # calculated based on pullman average and lapse rate
$precip_LR = 0.09195; # cm / m
$temp_LR = -3.20; # C / km

#print `r.mapcalc 'MASK = if(wsheds_all>0,1,0)' --o`;

# --------------------------------------------------------------------------------
#    SLOPE & ASPECT FOR SAM
# --------------------------------------------------------------------------------

system('r.slope.aspect elevation=el slope=slope aspect=aspect --o --quiet');

# --------------------------------------------------------------------------------
# READ IN WATERSHED AND RESERVOIR PROPERTIES
# --------------------------------------------------------------------------------

my $filename = 'wshed_res_properties.ini';

# Open properties file
open(my $WATERSHED, '<', $filename) || die "Cannot open the watershed properties file";

# Read each line of the file and assign values based on watershed ID
#print "\n \n";
#print "\n|----- READING WATERSHED PROPERTIES -----|\n";
#print "\n \n";
while (my $line = <$WATERSHED>) {
	chop($line);
	($wshed_id,$area_cells,$res_vol,$res_coeff) = split(' ', $line);
	
	if ($wshed_id > 0.0) {
		$area_{$wshed_id} = $area_cells * $gridsize * $gridsize; #square meters
		$res_vol_{$wshed_id} = $res_vol;
		$res_coeff_{$wshed_id} = $res_coeff;
		$base_flow_{$wshed_id} = 0.0;
		
		#print "	|---- Properties of watershed with ID: $wshed_id ----|\n";
		#print "		number of cells in watershed: $area_cells\n";
		#print "		cell size (m): $gridsize\n";
		#print "		total watershed area (m^2): $area_{$wshed_id}\n";
        #print "		res_vol: $res_vol_{$wshed_id}\n";
        #print "		res_coeff: $res_coeff_{$wshed_id}\n";
        #print "		base_flow: $base_flow_{$wshed_id}\n\n";
		#print "	|-------------------------------------------------|\n";
	};
}

# Close the watershed properties file
close($WATERSHED) || die "Cannot close the watershed properties file";


# OLD land use map legend:  1 = 100% forest cover, 2 = partial cut, 3 = clear cut
	#print `r.mapcalc 'landuse = 1.0' --o`;
# NEW land use map legend: 1 = water, 2 = urban/rock/barren/other, 3 = forest/woody wetlands, 4 = shrub, 5 = grass/grassy wetland/pasture, 6 = row crop

#	maximum canopy storage (cm) in descending land use number (6-1) order -- can range: 
	#	1. water: 0 cm
	#	2. urban/rock/barren/other: 0.1 cm
	#	3. forest/woody wetlands: .276 - 0.31 cm
	# 	4. shrub: 0.15 cm
	#	5. grass/grassy wetland/pasture: 0.15 cm
	#	6. row crop: 0.068 - 0.147 cm
print `r.mapcalc 'max_canopy_storage_amt = if(landuse==6.0,0.147,if(landuse==5.0,0.15,if(landuse==4.0,0.15,if(landuse==3.0,0.31,if(landuse==2.0,0.1,if(landuse==1.0,0.0,0.0))))))' --o`;

# values for our 6 land covers were mapped over from the existing forest, partial, clear, road classes --> #print `r.mapcalc 'canopy_cover = if(landuse==3.0,0.01,if(landuse==2.0,0.5,0.99))' --o`; # MZ 20200510 add canopy cover fraction for solar radiation calc
print `r.mapcalc 'canopy_cover = if(landuse==6.0,0.50,if(landuse==5.0,0.60,if(landuse==4.0,0.60,if(landuse==3.0,0.80,if(landuse==2.0,0,if(landuse==1.0,0,0.0))))))' --o`; ## DJ 3/14/23 --> changed the canopy_cover values to match those in email from Philip

#for now just mapped kfactor over from old landuse and assumed urban to be the same as clear cut (should probably be higher) --> #print `r.mapcalc 'kfactor = if(landuse==3.0,0.71,if(landuse==2.0,0.695,0.68))*$time_step/24.0' --o`;  #from output 68
#print `r.mapcalc 'kfactor = if(landuse==6.0, 0.695,if(landuse==5.0, 0.71,if(landuse==4.0,0.695,if(landuse==3.0,0.68,if(landuse==2.0,0.71,if(landuse==1.0,0.0,0.0))))))*$time_step/24.0' --o`;  #from output 68

# just using a walkover from old land use for now but we can improve on this --> #print `r.mapcalc 'tbase = if(landuse==3.0,1.7,if(landuse==2.0,1.9,2.1))' --o`;
#print `r.mapcalc 'tbase = if(landuse==6.0, 1.9,if(landuse==5.0, 1.7,if(landuse==4.0, 1.9,if(landuse==3.0, 2.1,if(landuse==2.0, 1.7,if(landuse==1.0, 1.7, 0.0))))))' --o`;

# land use values were mapped over from the old clear cut, partial, forest, and road system --> #print `r.mapcalc 'root_zone = if(landuse==3.0,soil_depth_A,soil_depth)' --o`;
print `r.mapcalc 'root_zone = if(landuse==6.0,soil_depth,if(landuse==5.0,soil_depth_A,if(landuse==4.0,soil_depth,if(landuse==3.0,soil_depth,if(landuse==2.0,0.0,if(landuse==1.0,0.0,0.0))))))' --o`;

#print `r.mapcalc 'wiltpt_amt = if(landuse==3.0,wiltpt_mc_A*soil_depth_A,wiltpt_mc_A*soil_depth_A+wiltpt_mc_B*soil_depth_B)' --o`;
print `r.mapcalc 'wiltpt_amt = if(landuse==6.0,wiltpt_mc_A*soil_depth_A+wiltpt_mc_B*soil_depth_B,if(landuse==5.0,wiltpt_mc_A*soil_depth_A,if(landuse==4.0,wiltpt_mc_A*soil_depth_A+wiltpt_mc_B*soil_depth_B,if(landuse==3.0,wiltpt_mc_A*soil_depth_A+wiltpt_mc_B*soil_depth_B,if(landuse==2.0,0.0,if(landuse==1.0,0.0,0.0))))))' --o`;

# ETreduction was done based on 80% for all classes except forest which is 70% and urban/barren/rock/water being 0 --> 	#print `r.mapcalc 'ETreduction_mc = if(landuse==3.0,fieldcap_amt_A*0.8/soil_depth_A,fieldcap_amt*0.8/soil_depth)' --o`;
print `r.mapcalc 'ETreduction_mc = if(landuse==6.0,fieldcap_amt*0.8/soil_depth,if(landuse==5.0,fieldcap_amt_A*0.8/soil_depth_A,if(landuse==4.0,fieldcap_amt*0.8/soil_depth,if(landuse==3.0,fieldcap_amt*0.7/soil_depth,if(landuse==2.0,0.0,if(landuse==1.0,0.0,0.0))))))' --o`;

# these values are just walked over from initial forest, partial, clear, road system assuming water,urban/barren are both clear cut (not perfect by any means) --> #print `r.mapcalc 'tmax_rain = if(landuse==3.0,3.1,if(landuse==2.0,1.8,-0.5))' --o`;
print `r.mapcalc 'tmax_rain = if(landuse==6.0,1.8,if(landuse==5.0,-0.5,if(landuse==4.0,1.8,if(landuse==3.0,3.1,if(landuse==2.0,-0.5,if(landuse==1.0,-0.5,0))))))' --o`;

# these were walked over from initial forest, partial, clear, road system where water,urban/barren are both clear cut values (could definitely be improved) --> #print `r.mapcalc 'tmin_snow = if(landuse==3.0,0.9,if(landuse==2.0,-1.0,-3.0))' --o`;
print `r.mapcalc 'tmin_snow = if(landuse==6.0,-1.0,if(landuse==5.0,-3.0,if(landuse==4.0,-1.0,if(landuse==3.0,0.9,if(landuse==2.0,-3.0,if(landuse==1.0,-3.0,0.0))))))' --o`;

# Initial conditions storage amount defined (% of saturation)
print `r.mapcalc 'storage_amt = sat_amt*0.4' --o`;
print `r.mapcalc 'storage_amt_A = sat_amt*0.4' --o`;
print `r.mapcalc 'storage_amt_A = (wiltpt_mc_A+(fieldcap_mc_A-wiltpt_mc_A)*0.2)*soil_depth_A' --o`;

#  storage_amt_ini was the storage amt of the previous run (9/31/1998)
#  storage_amt in stream cells is assumed at saturation
#  soil depth in stream cells is assumed at 1 cm A horizon, 19 cm B horizon
#  stream cells are defined as those cells having an upslope contributing area of 50 cells (using a 30m DEM)
print `r.mapcalc 'canopy_storage_amt = 0.0' --o`;
print `r.mapcalc 'swe = 0.0' --o`;

#  Initiation for SAM MZ 20190127
print `r.mapcalc 'snow.age = 1.0' --o`;
print `r.mapcalc 'swe.yesterday = 0.0' --o`;
print `r.mapcalc 'albedo = 0.2' --o`;
print `r.mapcalc 'liquid.water = 0.0' --o`;
print `r.mapcalc 'ice.content = 0.0' --o`;
print `r.mapcalc 'tsnow_surf = -2.0' --o`;
print `r.mapcalc 'u.surface = 24.0*tsnow_surf*(2.1*1000.0*min(swe/100.0,0.02)+1000.0*2.1*max(0.0,0.02-swe/100.0))' --o`; # convert to daily MZ 20190410 | should not be convert to daily becasue u.surface is in KJ/m^2. It is not a rate, it is the energy content of the surface layer. So the comment on 20190410 is wrong. MZ 20200330
print `r.mapcalc 'tsnow.pack = -0.5' --o`;
print `r.mapcalc 'u.total = 24.0*tsnow.pack*(2.1*1000.0*swe/100.0+2.1*1000.0*0.4)' --o`; # convert to daily MZ 20190410 | should not be convert to daily becasue u.surface is in KJ/m^2. It is not a rate, it is the energy content of the surface layer. So the comment on 20190410 is wrong. MZ 20200330

#  set the initial (t-1) time step, by MZ 2017.2.5
print `r.mapcalc 'mass_balance_total = 0.0' --o`;
print `r.mapcalc 'canopy_storage_amt_pre = 0.0' --o`;
print `r.mapcalc 'storage_amt_pre = storage_amt' --o`;

# added in for updating maximum storage map start values -- DJ, 2/20/23
print `r.mapcalc 'max_storage_A_horizon = 0.0' --o`;
print `r.mapcalc 'max_storage_B_horizon = 0.0' --o`;
print `r.mapcalc 'max_storage = 0.0' --o`;
print `r.mapcalc 'max_runoff = 0.0' --o`;

# annual and monthly map initialization
print `r.mapcalc 'runoff_1965 = 0.0' --o`;
print `r.mapcalc 'runoff_1966 = 0.0' --o`;
print `r.mapcalc 'runoff_1967 = 0.0' --o`;
print `r.mapcalc 'runoff_1968 = 0.0' --o`;
print `r.mapcalc 'runoff_1969 = 0.0' --o`;
print `r.mapcalc 'runoff_1970 = 0.0' --o`;
print `r.mapcalc 'runoff_1971 = 0.0' --o`;
print `r.mapcalc 'runoff_1972 = 0.0' --o`;
print `r.mapcalc 'runoff_1973 = 0.0' --o`;
print `r.mapcalc 'runoff_1974 = 0.0' --o`;
print `r.mapcalc 'runoff_1975 = 0.0' --o`;
print `r.mapcalc 'runoff_1976 = 0.0' --o`;
print `r.mapcalc 'runoff_1977 = 0.0' --o`;
print `r.mapcalc 'runoff_1978 = 0.0' --o`;
print `r.mapcalc 'runoff_1979 = 0.0' --o`;
print `r.mapcalc 'runoff_1980 = 0.0' --o`;
print `r.mapcalc 'runoff_1981 = 0.0' --o`;
print `r.mapcalc 'runoff_1982 = 0.0' --o`;
print `r.mapcalc 'runoff_1983 = 0.0' --o`;
print `r.mapcalc 'runoff_1984 = 0.0' --o`;
print `r.mapcalc 'runoff_1985 = 0.0' --o`;
print `r.mapcalc 'runoff_1986 = 0.0' --o`;
print `r.mapcalc 'runoff_1987 = 0.0' --o`;
print `r.mapcalc 'runoff_1988 = 0.0' --o`;

print `r.mapcalc 'runoff_feb_1965 = 0.0' --o`;
print `r.mapcalc 'runoff_feb_1966 = 0.0' --o`;
print `r.mapcalc 'runoff_feb_1967 = 0.0' --o`;
print `r.mapcalc 'runoff_feb_1968 = 0.0' --o`;
print `r.mapcalc 'runoff_feb_1969 = 0.0' --o`;
print `r.mapcalc 'runoff_feb_1970 = 0.0' --o`;
print `r.mapcalc 'runoff_feb_1971 = 0.0' --o`;
print `r.mapcalc 'runoff_feb_1972 = 0.0' --o`;
print `r.mapcalc 'runoff_feb_1973 = 0.0' --o`;
print `r.mapcalc 'runoff_feb_1974 = 0.0' --o`;
print `r.mapcalc 'runoff_feb_1975 = 0.0' --o`;
print `r.mapcalc 'runoff_feb_1976 = 0.0' --o`;
print `r.mapcalc 'runoff_feb_1977 = 0.0' --o`;
print `r.mapcalc 'runoff_feb_1978 = 0.0' --o`;
print `r.mapcalc 'runoff_feb_1979 = 0.0' --o`;
print `r.mapcalc 'runoff_feb_1980 = 0.0' --o`;
print `r.mapcalc 'runoff_feb_1981 = 0.0' --o`;
print `r.mapcalc 'runoff_feb_1982 = 0.0' --o`;
print `r.mapcalc 'runoff_feb_1983 = 0.0' --o`;
print `r.mapcalc 'runoff_feb_1984 = 0.0' --o`;
print `r.mapcalc 'runoff_feb_1985 = 0.0' --o`;
print `r.mapcalc 'runoff_feb_1986 = 0.0' --o`;
print `r.mapcalc 'runoff_feb_1987 = 0.0' --o`;
print `r.mapcalc 'runoff_feb_1988 = 0.0' --o`;

print `r.mapcalc 'Psat_1965 = 0.0' --o`;
print `r.mapcalc 'Psat_1966 = 0.0' --o`;
print `r.mapcalc 'Psat_1967 = 0.0' --o`;
print `r.mapcalc 'Psat_1968 = 0.0' --o`;
print `r.mapcalc 'Psat_1969 = 0.0' --o`;
print `r.mapcalc 'Psat_1970 = 0.0' --o`;
print `r.mapcalc 'Psat_1971 = 0.0' --o`;
print `r.mapcalc 'Psat_1972 = 0.0' --o`;
print `r.mapcalc 'Psat_1973 = 0.0' --o`;
print `r.mapcalc 'Psat_1974 = 0.0' --o`;
print `r.mapcalc 'Psat_1975 = 0.0' --o`;
print `r.mapcalc 'Psat_1976 = 0.0' --o`;
print `r.mapcalc 'Psat_1977 = 0.0' --o`;
print `r.mapcalc 'Psat_1978 = 0.0' --o`;
print `r.mapcalc 'Psat_1979 = 0.0' --o`;
print `r.mapcalc 'Psat_1980 = 0.0' --o`;
print `r.mapcalc 'Psat_1981 = 0.0' --o`;
print `r.mapcalc 'Psat_1982 = 0.0' --o`;
print `r.mapcalc 'Psat_1983 = 0.0' --o`;
print `r.mapcalc 'Psat_1984 = 0.0' --o`;
print `r.mapcalc 'Psat_1985 = 0.0' --o`;
print `r.mapcalc 'Psat_1986 = 0.0' --o`;
print `r.mapcalc 'Psat_1987 = 0.0' --o`;
print `r.mapcalc 'Psat_1988 = 0.0' --o`;

#____________________________________________________________________________________
#z
# START READING WEATHER DATA
#____________________________________________________________________________________

#  This set of commands splits a tab delimited array using a while loop
open ($WEATHER, '/Users/duncanjurayj/GrassWorkSpace/pullman_historical.csv') || die "Can't open file\n";


while (<$WEATHER>) {
	chop($_);
	($date,$year,$month,$day,$doy,$tmax,$tmin,$tavg,$tdew,$precip,$pet,$hour_1,$hour_6,$hour_12,$hour_18,$l_turb,$cloud,$cc_water,$cc_urban,$cc_forest,$cc_shrub,$cc_grass,$cc_row_crop,$rh_water,$rh_urban,$rh_forest,$rh_shrub,$rh_grass,$rh_row_crop,$output) = split(' ',$_);

	#print "\n \n";
	#print "\n|----- DAILY WEATHER READ FOR $date -----|\n";
	#print "\n \n";
	#print "\n	WEATHER SUMMARY:\n";
	#print "\n	| Date $date | Year = $year | Tair = $tavg °C | Precip = $precip cm | PET = $pet cm | Tdew = $tdew | Output = $output\n";
	#print "\n	|-------------------------------------------------------------------------------------------------------------------------| \n";

	@hourly_tmp_array = ($hour_1,$hour_6,$hour_12,$hour_18);

#   changed snowmelt so that clear cut melt is the same as partial cut melt
#	print `r.mapcalc 'kfactor = if(landuse==3.0,0.734,if(landuse==2.0,0.734,0.719))*$time_step/24.0' --o`;
#	print `r.mapcalc 'tbase = if(landuse==3.0,1.54,if(landuse==2.0,1.54,2.14))' --o`;
	#print `r.mapcalc 'kfactor = if(landuse==3.0,0.71,if(landuse==2.0,0.695,0.68))*$time_step/24.0' --o`; # Degree-day
	#print `r.mapcalc 'kfactor = if(landuse==6.0,0.0,if(landuse==5.0,0.0,if(landuse==4.0,0.0,if(landuse==3.0,0.71,if(landuse==2.0,0.695,if(landuse==1.0,0.68,0.0))))))' --o`;

	#print `r.mapcalc 'tbase = if(landuse==3.0,1.7,if(landuse==2.0,1.9,2.1))' --o`; # Degree-day
	#print `r.mapcalc 'tbase = if(landuse==6.0,0.0,if(landuse==5.0,0.0,if(landuse==4.0,0.0,if(landuse==3.0,1.7,if(landuse==2.0,1.9,if(landuse==1.0,2.1,0.0))))))' --o`;

    ##print `r.mapcalc 'max_canopy_storage_amt = if(landuse==6.0,0.147,if(landuse==5.0,0.15,if(landuse==4.0,0.15,if(landuse==3.0,0.31,if(landuse==2.0,0.1,if(landuse==1.0,0.0,0.0))))))' --o`;
	
	# this is redundant because it's the same as above but i'm not sure why it's here
	##print `r.mapcalc 'canopy_cover = if(landuse==6.0,0.50,if(landuse==5.0,0.60,if(landuse==4.0,0.60,if(landuse==3.0,0.80,if(landuse==2.0,0,if(landuse==1.0,0,0.0))))))' --o`; ## DJ 3/14/23 --> changed the canopy_cover values to match those in email from Philip

	# this is redundant because it's the same as above but i'm not sure why it's here
	##print `r.mapcalc 'root_zone = if(landuse==6.0,soil_depth,if(landuse==5.0,soil_depth_A,if(landuse==4.0,soil_depth,if(landuse==3.0,soil_depth,if(landuse==2.0,0.0,if(landuse==1.0,0.0,0.0))))))' --o`;

	# this is redundant because it's the same as above but i'm not sure why it's here
	##print `r.mapcalc 'wiltpt_amt = if(landuse==6.0,wiltpt_mc_A*soil_depth_A+wiltpt_mc_B*soil_depth_B,if(landuse==5.0,wiltpt_mc_A*soil_depth_A,if(landuse==4.0,wiltpt_mc_A*soil_depth_A+wiltpt_mc_B*soil_depth_B,if(landuse==3.0,wiltpt_mc_A*soil_depth_A+wiltpt_mc_B*soil_depth_B,if(landuse==2.0,0.0,if(landuse==1.0,0.0,0.0))))))' --o`;

	# ETreduction was done based on 80% for all classes except forest which is 70% and urban/barren/rock/water being 0
	##print `r.mapcalc 'ETreduction_mc = if(landuse==6.0,fieldcap_amt*0.8/soil_depth,if(landuse==5.0,fieldcap_amt_A*0.8/soil_depth_A,if(landuse==4.0,fieldcap_amt*0.8/soil_depth,if(landuse==3.0,fieldcap_amt*0.7/soil_depth,if(landuse==2.0,0.0,if(landuse==1.0,0.0,0.0))))))' --o`;


#  --------------------------------------------------------------------
#  --------------------------------------------------------------------
#____________________________________________________________________________________
#
#  1.  SNOWMELT AND INTERCEPTION ALGORITHM
#____________________________________________________________________________________

#  rain and snow are in (cm)
#  Mica snotel elevation = 1447.8 m, Mica snotel mean annual precip = 145.7 cm
#  Precip lapse rate = 0.08136 cm / meter of elevation
#  Mean Annual Precip at sea level (0 m elevation or intercept) = 27.9 cm
#  Temperature lapse rate is -5.31 C/km
#  Potential evapotranspiration warm calculated using the Hargreaves method using 
#  Mica snotel temperature.  A PET lapse rate was determined by calculating PET
#  for temperatures extrapolated to St Maries elevation (675.4 m) and fitting

# rain and snow are in (cm)
# Moscow Mountain SNOTEL elevation = 1443 m, Pullman Station elevation = 784 m
# precipitation lapse rate between sites = 0.09195 cm / meter of elevation --> seems low but move forward with this for now
# temperature lapse rate between sites = - -3.20 C/km --> also on the lower end it seems, but fine to use for now we can always modify both values --> divide by 1000 to work with meters
# precipitation throughout watershed is calculated with: (elevation * precip_lapse_rate + mean_annual_sea_level_precip) / (higher_station_mean_annual_precip) * (actual_precip_data)
# temperature throughout watershed is calculated with: temp - (temp_lapse_rate * (el - low_site_elevation)) --> this assumes data is taken from the lower site (Pullman site)
#  an equation between PET and elevation.  PET_moscow_mountain = 1.8281 * PET_pullman + -0.2486 (cm)


# this equation may be wrong 
print `r.mapcalc 'precip = (el*$precip_LR+$mean_annual_sea_level_precip) / ($mean_annual_onsite_precip) * $precip' --o`;
print `r.mapcalc 'tavg = $tavg + ($temp_LR/1000) * (el - $low_site_el)' --o`;
print `r.mapcalc 'tdew = $tdew + ($temp_LR/1000) * (el - $low_site_el - el)' --o`;
print `r.mapcalc 'rain = if(tavg>tmax_rain,precip,if(tavg<tmin_snow,0.0,(tavg-tmin_snow)/(tmax_rain-tmin_snow)*precip))' --o`;
print `r.mapcalc 'snow = precip-rain' --o`;

# THESE ARE STILL FROM MICA!!! NEED TO UPDATE -- Duncan 3/12/23 --> these have been updated to get MM from Pullman: Duncan 3/13/23
print `r.mapcalc 'pet_data = (1.8281*$pet-0.2486)' --o`;
print `r.mapcalc "pet = (pet_data - $pet)*(el - $low_site_el)/($low_site_el - $high_site_el)+pet_data" --o`;
#print `r.mapcalc "pet = (pet_data-$pet)*($low_site_el-el)/($high_site_el-$low_site_el)+pet_data" --o`;
# interception is calculated for Spruce trees based on work by Lankreijer et al. 1999, Agricultural and Forest
# Meteorology 98-99:595.  Max Storage of Canopy was taken as 2.0 mm.  During rain, evaporation is 0.04 mm/hr.
# When there is no rain, for the time being, it is assumed that evaporation = 50% of PET.

#  Assume canopy storage is 3 mm for full canopy, 1.5 mm for partial canopy
#  Assume evaporation of canopy storage is at the potential rate

print `r.mapcalc 'canopy_storage_amt = canopy_storage_amt + rain' --o`;
print `r.mapcalc 'throughfall = if(canopy_storage_amt>max_canopy_storage_amt,canopy_storage_amt-max_canopy_storage_amt,0.0)' --o`;
print `r.mapcalc 'canopy_storage_amt = canopy_storage_amt-throughfall' --o`;
print `r.mapcalc 'canopy_evap = if(rain>0.0,min(canopy_storage_amt,pet),min(canopy_storage_amt,pet))' --o`;
print `r.mapcalc 'pet = max(0.0,pet-canopy_evap)' --o`;
print `r.mapcalc 'canopy_storage_amt = canopy_storage_amt-canopy_evap' --o`;

#  Degree Day snowmelt
#  USACE recommends forest = 0.23 cm/C/day, 0 C, open = 0.27 cm/C/day, Tb=0
#  Fitting Mica Creek snotel data  Ksnow = 0.719 cm/C/day, Tb = 2.14
#  For open areas assume USACE difference between forest and open melt accurate
#  Using solver in excel fit open area snow parameters by preserving USACE difference
#  Assumed melt in partial cut increased by 1/3 of the difference between open and forest areas 
#  Fitted open area parameters Ksnow_open = 0.764 cm/C/day, Tb_open = 0.42
#  Fitted partial cut parameters Ksnow_partial = 0.734 cm/C/day, Tb_partial = 1.54

#print `r.mapcalc 'snowmelt = if(swe+snow-max(0.0,kfactor*(tavg-tbase))<0.0,swe+snow,max(0.0,kfactor*(tavg-tbase)))' --o`;
#print `r.mapcalc 'swe = swe+snow-snowmelt' --o`;

#print "\n \n";
print "\n|----- SNOW ACCUMULATION AND MELT MODEL -----|\n";
#  ------------------------------- SAM --------------------------------------
#  Snow accumulation and melt (SAM) model
#  Modified by MZ on 20190126
#  albedo approximated from DHSVM relationship
#  assumes the albedo of the soil is 0.2
#  albedo (100 x %)

#print `r.mapcalc 'rh = if(swe>0,if(landuse==1 || landuse==2,$rh_veg,$rh_snow),$rh_veg)' --o`; # modified MZ 20190205
#print `r.mapcalc 'rh = if(swe>0,$rh_snow,$rh_veg)' --o`; # modified MZ 20190205
# LAZY FIX RIGHT NOW WHERE RH_SNOW IS RH_URBAN BUT THIS WILL NEED TO CHANGE (THOUGH THEY WON'T BE TOO DIFFERENT) --> ADDITIONALL YJUST USING ROW CROP RIGHT NOW
print `r.mapcalc 'rh = if(swe>0,$rh_urban/(1.0-(canopy_cover/3.0)),$rh_row_crop/(1.0-(canopy_cover/3.0)))' --o`; # add canopy cover effects on rh MZ 20200601
print `r.mapcalc 'snow.age = if(snow>0.0 && throughfall==0.0,1.0,snow.age+1.0)' --o`; # modified MZ 20190210
print `r.mapcalc 'albedo = if(swe.yesterday+snow>0.0,min(0.95,0.7383*snow.age^(-0.1908)),0.2)' --o`;

#  Calculate the diffuse transmissivity following Bristow and Campbell (1985)
#  The clear sky transmissivity for Troy, ID is 0.75
print `r.mapcalc 't_diff = if($cloud==1.0,0.1,0.75*(1-$cloud)*(1-exp(-0.6*$cloud/((0.75-0.4)*(1-$cloud))))) * watershed' --o`;

#  Calculate the real-sky diffuse radiation coefficient
#  The clear sky diffuse transmissivity for the Palouse is 0.10 (Bristow and Campbell, 1985)
print `r.mapcalc 'coefdh = t_diff/0.10' --o`;

#  Calculate the real-sky beam radiation coefficient
#  The clear sky transmissivity for Troy, ID is 0.75
print `r.mapcalc 'coefbh = max(0.0,(1-$cloud)-t_diff/0.75)' --o`;

#  Calculate the beam, diffuse, and reflected irradiance at time $s_time
#  Change to daily time step so the output of beam, diffuse and reflected radiation is in Wh/m2/day, MZ 20190127
print `r.sun elevation=el aspect=aspect slope=slope linke_value=$l_turb albedo=albedo coeff_bh=coefbh coeff_dh=coefdh beam_rad=beam_rad diff_rad=diff_rad refl_rad=refl_rad day=$doy nprocs=12 --o --quiet`;

#  Rad_tot is the total radiation in units of KJ/m^2
#  Convert to daily value by multiply 24, MZ 20170127; Update: Do not multiply 24 because the output of beam etc is already daily in Wh/m2/day. 1Wh/m2/day = 3.6KJ/m2/day
print `r.mapcalc 'q.srad = (1.0-canopy_cover)*(1.0-albedo)*(beam_rad+diff_rad+refl_rad) * 3600.0/1000.0' --o`; # add canopy cover impact MZ 20200510

#  q.lw is the net heat energy added to the snowpack by longwave rad. 
#  assume the emissivity of snow is 0.98
#  Stefan Boltzmann constant = 5.67 x 10^-8 W/m^2/K^4 
#  (NOTE: THIS CONSTANT MUST BE CONVERTED TO J/m^2/K^4/hr or /day ACCORDING TO THE DESIRED TIMESTEP)
#  Cloudiness determined in a spreadsheet by the 1 - measured srad divided by predicted clear sky srad 
#  Emissivity of cloudy skies calculated from the Unsworth and Montieth (1975)
#  relationship e=(0.72+0.005*Tavg)*(1-0.84*Cloud fraction)+0.84*Cloud fraction
#  q.lw (KJ/m^2)
#  Convert to daily value by multiply 24, MZ 20190127

#print `r.mapcalc 'q.lw = 5.67*(10.0^-8.0)*3600.0*24.0/1000.0*(max(0.72,(0.72+0.005*tavg)*(1.0-0.84*$cloud)+0.84*$cloud)*(tavg+273.15)^4.0-0.98*(tsnow_surf+273.15)^4.0)' --o`;
#print `r.mapcalc 'q.lw = 5.67*(10.0^-8.0)*3600.0*24.0/1000.0*(max(0.72,(9.2*(10^-6.0)*(tavg+273.15)^2.0)*(1.0-0.84*$cloud)+0.84*$cloud)*(tavg+273.15)^4.0-0.98*(tsnow_surf+273.15)^4.0)' --o`; # MZ 20190408 from Campbell an introductino to environmental biophysics eqn (10.11) to calculate clear sky emissivity e(0)=9.2*10^-6*Tavg^2, where Tavg is in Kelvin
print `r.mapcalc 'q.lw = 5.67*(10.0^-8.0)*3600.0*24.0/1000.0*(max(0.72,((9.2*(10^-6.0)*(tavg+273.15)^2.0)*(1.0-0.84*$cloud)+0.84*$cloud)*(1-canopy_cover)+0.92*canopy_cover)*(tavg+273.15)^4.0-0.98*(tsnow_surf+273.15)^4.0)' --o`; # add canopy cover impact MZ 20200510

#  act. air vapor density calculated using tdew 
#  esat = e((16.78*tdew-116.9)/(tdew+237.3))                     units = (kPa)
#  act. air vapor density (kg/m^3) = esat/(273 + tdew)/0.4615    units = (kg/m^3)
#  Thermodynamic vapor constant = 0.4615                         units =  KJ/(kg K)
#  Vapor density of the snow is assumed saturated if snow surface temperature = 0
#  Otherwise the snow surface vapor density is calculated using the freezing point depression equation (Flerchinger)
#  vap.d.air, vap.d.snow_sat, vap.d.snow (kg/m^3)

print `r.mapcalc 'vap.d.snow_sat = exp((16.78*tsnow_surf-116.9)/(tsnow_surf+237.3))/(273.15+tsnow_surf)/0.4615' --o`;
print `r.mapcalc 'vap.d.snow = if(tsnow_surf<0.0,vap.d.snow_sat*exp(335.0*tsnow_surf/(0.4615*(273.15+tsnow_surf)^2.0)),vap.d.snow_sat)' --o`;
print `r.mapcalc 'vap.d.air = exp((16.78*tdew-116.9)/(tdew+237.3))/(tdew+273.15)/0.4615' --o`;

# q.vap is the heat from convective vapor exchange (condensation/evaporation)
# latent heat of vaporization = 2500 KJ/kg
# latent heat of fusion = 335 KJ/kg
# It is assumed initially that there is a vapor/solid or solid/vapor change 
# The refreeze calculation will correct this assumption if there is only a vapor/liquid or liquid/vapor change
# the wind roughness $rh is inputed in units of s/m, 
#  NOTE:  THIS NEEDS TO BE CONVERTED TO DAY/M OR HR/M ACORDING TO THE DESIRED TIMESTEP
# q.latent (KJ/m^2)

print `r.mapcalc 'q.latent = (2500.0+335.0)*(vap.d.air-vap.d.snow)/(rh/(3600.0*24.0))' --o`; # Duncan Jurayj - divided by 1000 --> q.latent maps went up to 50,000 (maybe in J rather than KJ?) - removed this to test something
print `r.mapcalc 'condens = if((swe.yesterday+snow)>0.0,q.latent/((2500.0+335.0)*1000.0)*100.0,0.0)' --o`; # in cm MZ 20170128; get daily average value by dividing 24 MZ 20190410
#print `r.mapcalc 'condens = if((swe.yesterday+snow)>0.0,if(q.latent/((2500.0+335.0)*1000.0)*100.0>0,q.latent/((2500.0+335.0)*1000.0)*100.0, max(q.latent/((2500.0+335.0)*1000.0)*100.0,-(swe.yesterday+snow))),0.0)' --o`; # in cm MZ 20170128; get daily average value by dividing 24 MZ 20190410

# q.cond is the heat from convective temp. gradients between air and snow
# density of air 1.29 kg/m^3, heat capacity of air 1 kJ/kg/C
# the wind roughness $rh is inputed in units of s/m
#  NOTE:  THIS NEEDS TO BE CONVERTED TO DAY/M OR HR/M ACORDING TO THE DESIRED TIMESTEP
# q.cond (KJ/m^2)

print `r.mapcalc 'q.sensible = 1.0*1.29*(tavg-tsnow_surf)/(rh/(3600.0*24.0))' --o`; # convert to KJ/(m^2 day) MZ 20200329
#print `r.mapcalc 'q.sensible = 1.0*1.29*(tavg-tsnow_surf)/(rh/3600)' --o`; #

# q.rain.ground is the combined conductive heat from rain, snow, and ground 
# assumes ground melt constant at 0.0 KJ/m^2/day
# NOTE:  THIS CONSTANT NEEDS TO BE CONVERTED TO A DAILY OR HOURLY TIMESTEP
# assumes heat capacity of water is 4.2 KJ/(kg C)
# assumes heat capacity of snow is 2.1 KJ/(kg C)
# assumes the density of water is 1000 kg/m^3
# units of rain are cm (remember to convert cm to m)
# q.rain.ground (KJ/m^2)

print `r.mapcalc 'q.rain.ground = 4.2*1000.0*(throughfall/100.0)*(max(0.0,tavg)-tsnow_surf)+2.1*1000.0*(snow/100.0)*(min(0,tavg)-tsnow_surf)+0.0' --o`; # unit convert to daily done MZ 20190409

# Total incoming energy flux at snow surface.  
# q.* are in units of (KJ/m^2)
# q.total (KJ/m^2)

print `r.mapcalc 'q.total = q.srad+q.lw+q.sensible+q.latent+q.rain.ground' --o`; 

# refreeze
# If positive it is the amount of ice added to the snowpack by the freezing of liquid water 
#     in order to raise the snow pack temperature to 0 C 
# If negative it is the amount of ice which must thaw to cool the temperature of the pack to 0 C
#  refreeze (cm) 

print `r.mapcalc 'refreeze = min((throughfall+liquid.water),max(0.0,-1.0*100.0*(u.total+q.total)/(335.0*1000.0)))-min((snow+ice.content+condens),max(0.0,100.0*(u.total+q.total)/(335.0*1000.0)))' --o`;

#  Internal energy of the snow + soil layer (KJ/m^2) relative to a snowpack at 0 degrees C
#  A snowpack at 0 degrees C has an internal energy of 0 KJ/m^2
#  U = U + Qin + refreeze 

print `r.mapcalc 'u.total = u.total+q.total+refreeze/100.0*335.0*1000.0' --o`;

# ice.content
# if U > 0 then ice.content = 0
# if U < 0 then everything frozen
# cm

print `r.mapcalc 'ice.content = if(u.total>0.000001,0.0,ice.content+snow+refreeze+condens)' --o`;

# liquid.water
# if U > 0 or U < 0 then liquid water = 0
# Otherwise it is assumed that the water holding capacity is 3% of the SWE
# cm

print `r.mapcalc 'liquid.water = if(u.total>0.000001||u.total<-0.000001,0.0,min(0.03/(1.0-0.03)*ice.content,liquid.water+throughfall-refreeze))' --o`; # Modified based on Erin SAM paper MZ 20190410

# snow water equivalent (cm)

print `r.mapcalc 'swe = liquid.water+ice.content' --o`;

#  melt.water  
#  the amount of liquid water which leaves the pack and enters the soil
#  If a pack is not present it is equal to the rainfall

print `r.mapcalc 'snowmelt = (swe.yesterday - swe) + snow+condens' --o`;

print `r.mapcalc 'swe.y.save = swe.yesterday' --o`;
print `r.mapcalc 'swe.yesterday = swe' --o`;

# snow pack temperature
# It is not required to calculate snow pack temperature in the model
# Rather the energy content of the snow pack is tracked (u.total)
# If it is desired to output snow pack temperature then use the following equation
# Assuming the soil damping depth is 0.4 m, the bulk density of the soil is 1000 Kg/m^3
# Assuming the specific heat of the soil is 2.1 KJ/kg/C
# Assumes the specific heat of the snow is 2.1 KJ/kg/C and the density of water is 1000 Kg/m^3
#
#r.mapcalc tsnow.pack='u.total/(2.1*1000.0*swe/100.0+2.1*1000.0*0.4)'

# u.surface (KJ/m^2)
# The energy content of the surface layer
# Assumes the bulk density of the soil is 1000.0 Kg/m^3
# Latent heat of fusion = 335 KJ/kg, density of water = 1000.0 Kg/m^3
# Assumes the snow skin thickness is 0.05 m
# Assumes the snow skin damping depth is 0.15 m (related to the thermal conductivity of the snow)
# Assumes the soil skin damping depth is 0.3 m (related to the thermal conductivity of the soil)
# If u.total indicates melting (=0) then u.surface = 0
# If swe > 0 then u.surface must be <= 0 C or 0 KJ/m^2 

print `r.mapcalc 'u.surface = if(u.total==0.0,0.0,if(swe>0.0,min(0.0,u.surface+(q.total+refreeze/100.0*335.0*1000.0)*(min(swe/100.0,0.02)/0.05+max(0.02-swe/100.0,0.0)/0.3)),u.surface+(q.total+refreeze/100.0*335.0*1000.0)*(min(swe/100.0,0.02)/0.05+max(0.02-swe/100.0,0.0)/0.3)))' --o`;

# snow surface temperature (C)
# Assumes the specific heat of the soil is 2.1 KJ/kg/C
# Assumes the bulk density of the soil is 1000.0 Kg/m^3
# Assumes the density of water is 1000.0 Kg/m^3
# Assumes the specific heat of the snow is 2.1 KJ/kg/C
# Assumes the snow skin thickness is 0.02 m
# Assumes the snow skin damping depth is 0.15 m
# Assumes the soil skin damping depth is 0.3 m

print `r.mapcalc 'tsnow_surf = u.surface/24.0/(2.1*1000.0*min(swe/100.0,0.02)+1000.0*2.1*max(0.0,0.02-swe/100.0))' --o`;

#  -----------------------------End of SAM ----------------------------------
#print "\n|----- END OF SNOW ACCUMULATION AND MELT MODEL -----|\n";
#print "\n \n";

# assume road width is 5 m and assume road surface storage is 1.5 cm 
# assume once surface storage is exceeded all runoff reaches stream
# print `echo "road_runoff=if(roads==1.0,max(0.0,snowmelt+throughfall-1.5)*5.0/$gridsize,0.0)" | r.mapcalc`;
#print `r.mapcalc 'road_runoff = if(roads==1.0,max(0.0,snowmelt+throughfall-1.5)*5.0/$gridsize,0.0)' --o`;
print `r.mapcalc 'road_runoff = if(landuse==2.0,max(0.0,snowmelt+throughfall-1.5)/$gridsize,0.0)' --o`; #removed the x5 part of the road runoff (becuase grid size is already 30)
print `r.mapcalc 'water_input = snowmelt+throughfall' --o`;
print `r.mapcalc 'snowmelt = if(water_input>0.0,snowmelt*(water_input-road_runoff)/water_input,snowmelt)' --o`;
print `r.mapcalc 'throughfall = if(water_input>0.0,throughfall*(water_input-road_runoff)/water_input,throughfall)' --o`;

print `r.mapcalc 'temp_sum = 0.0' --o`;
print `r.mapcalc 'actualET_daily_flow = 0.0' --o`;
print `r.mapcalc 'perc_daily_flow = 0.0' --o`;
print `r.mapcalc 'runoff_daily_flow = 0.0' --o`;
print `r.mapcalc 'lateral_daily_out = 0.0' --o`;
print `r.mapcalc 'lateral_daily_in = 0.0' --o`;
print `r.mapcalc 'input_daily = 0.0' --o`;
print `r.mapcalc 'SM_test = 0.0' --o`; 

print `r.mapcalc 'mass_daily_balance = 0.0' --o`; 
print `r.mapcalc 'input_daily_balance = 0.0' --o`;


	foreach $hrly_tmp (@hourly_tmp_array)
		{
#			print `r.mapcalc 'temp = $hrly_tmp' --o`; # MZ 20170210
			print `r.mapcalc 'temp = $hrly_tmp-($temp_LR/1000)*($low_site_el - el)' --o`;
			print `r.mapcalc 'temp_sum = max(0.0,temp)+temp_sum' --o`;

		};
	#print "\n \n";
	#print "\n|----- 6 HOUR TEMPERATURE ARRAY -----|\n";
	#print "\n \n";
	foreach $hrly_tmp (@hourly_tmp_array)
		{
#			print `r.mapcalc 'temp = $hrly_tmp' --o`; # MZ 20170210
			print `r.mapcalc 'temp = $hrly_tmp-($temp_LR/1000)*($low_site_el - el)' --o`;

			#print "$hrly_tmp \n";


	#____________________________________________________________________________________
	#
	#  2.  DISTRIBUTE WATER TO EACH SOIL LAYER
	#____________________________________________________________________________________

	#    Storage is distributed assuming that all water infiltrates.  When total field capacity of
	#    the entire soil profile is not exceeded, moisture is evenly distributed.
	#    When total field capacity is exceeded, saturation begins

	print `r.mapcalc 'input = snowmelt*max(temp,0.0)/if(temp_sum==0.0,1.0,temp_sum)+throughfall*$temp_time_step/24.0' --o`;

	print `r.mapcalc 'input_daily = input_daily + input' --o`;

	print `r.mapcalc 'storage_amt = storage_amt+input' --o`;

	print `r.mapcalc 'storage_amt_A_tmp = if(storage_amt<fieldcap_amt,min(input+storage_amt_A,fieldcap_amt_A),storage_amt_A)' --o`;

	print `r.mapcalc 'storage_amt_B = if(storage_amt<fieldcap_amt,storage_amt-storage_amt_A_tmp,if((storage_amt-fieldcap_amt_A)<sat_mc_B*soil_depth_B,storage_amt-fieldcap_amt_A,sat_mc_B*soil_depth_B))' --o`;

	print `r.mapcalc 'storage_amt_A = if(storage_amt<fieldcap_amt,storage_amt_A_tmp,if((storage_amt-storage_amt_B)<sat_mc_A*soil_depth_A,storage_amt-storage_amt_B,sat_mc_A*soil_depth_A))' --o`;

	#____________________________________________________________________________________
	#
	#    3. SUBSURFACE LATERAL FLOW
	#____________________________________________________________________________________

	#    Lateral flow is based on Darcy's Law, with gradient
	#    equal to land slope, and direction maps (north, northeast, etc)
	#    calculated from the elevation map in the program smr.setup.

	#    Effective conductivity effK is calculated with Bresler's formula for
	#    unsaturated conductivity.  Kfc is conductivity at field capacity.
	#    Effective conductivity combines saturated and unsaturated.
	#    Assumes that in unsaturated flow is unaffected by macro-pores

	print `r.mapcalc 'effK_A = if(storage_amt_A<fieldcap_amt_A,Ksat_matrix_A*exp((-13.0/sat_mc_A)*(sat_mc_A-storage_amt_A/soil_depth_A)),if(storage_amt_A>=sat_mc_A*soil_depth_A,Ksat_mpore_A,(Ksat_mpore_A-Kfc_A)*(storage_amt_A/soil_depth_A-fieldcap_amt_A/soil_depth_A)/(sat_mc_A-fieldcap_amt_A/soil_depth_A)+Kfc_A))*$temp_time_step/24.0' --o`;
	print `r.mapcalc 'effK_B = if(storage_amt_B<fieldcap_amt_B,Ksat_matrix_B*exp((-13.0/sat_mc_B)*(sat_mc_B-storage_amt_B/soil_depth_B)),if(storage_amt_B>=sat_mc_B*soil_depth_B,Ksat_mpore_B,(Ksat_mpore_B-Kfc_B)*(storage_amt_B/soil_depth_B-fieldcap_amt_B/soil_depth_B)/(sat_mc_B-fieldcap_amt_B/soil_depth_B)+Kfc_B))*$temp_time_step/24.0' --o`;


	#******************************** Layer A B FIX ***************************************
	# MZ 20171205

	print `r.mapcalc 'lateral_flow_A = min(effK_A*slope/100.0*soil_depth_A/$gridsize/100.0,storage_amt_A)' --o`;
	print `r.mapcalc 'lateral_flow_A = if(flowunits==0.0,0.0,lateral_flow_A)' --o`;
	print `r.mapcalc 'lateral_flow_B = min(effK_B*slope/100.0*soil_depth_B/$gridsize/100.0,storage_amt_B)' --o`;
	print `r.mapcalc 'lateral_flow_B = if(flowunits==0.0,0.0,lateral_flow_B)' --o`;
	print `r.mapcalc 'lateral_flow = lateral_flow_A + lateral_flow_B' --o`;

	print `r.mapcalc 'storage_amt_A = storage_amt_A - lateral_flow_A  + (if(isnull(lateral_flow_A[-1,0]),0.0,lateral_flow_A[-1,0])*north + if(isnull(lateral_flow_A[-1,1]),0.0,lateral_flow_A[-1,1])*northeast+if(isnull(lateral_flow_A[0,1]),0.0,lateral_flow_A[0,1])*east+if(isnull(lateral_flow_A[1,1]),0.0,lateral_flow_A[1,1])*southeast+if(isnull(lateral_flow_A[1,0]),0.0,lateral_flow_A[1,0])*south+if(isnull(lateral_flow_A[1,-1]),0.0,lateral_flow_A[1,-1])*southwest + if(isnull(lateral_flow_A[0,-1]),0.0,lateral_flow_A[0,-1])*west +if(isnull(lateral_flow_A[-1,-1]),0.0,lateral_flow_A[-1,-1])*northwest)' --o`;

	print `r.mapcalc 'storage_amt_B = storage_amt_B - lateral_flow_B  + (if(isnull(lateral_flow_B[-1,0]),0.0,lateral_flow_B[-1,0])*north + if(isnull(lateral_flow_B[-1,1]),0.0,lateral_flow_B[-1,1])*northeast+if(isnull(lateral_flow_B[0,1]),0.0,lateral_flow_B[0,1])*east+if(isnull(lateral_flow_B[1,1]),0.0,lateral_flow_B[1,1])*southeast+if(isnull(lateral_flow_B[1,0]),0.0,lateral_flow_B[1,0])*south+if(isnull(lateral_flow_B[1,-1]),0.0,lateral_flow_B[1,-1])*southwest + if(isnull(lateral_flow_B[0,-1]),0.0,lateral_flow_B[0,-1])*west +if(isnull(lateral_flow_B[-1,-1]),0.0,lateral_flow_B[-1,-1])*northwest)' --o`;

	print `r.mapcalc 'storage_amt = storage_amt - lateral_flow  + (if(isnull(lateral_flow[-1,0]),0.0,lateral_flow[-1,0])*north + if(isnull(lateral_flow[-1,1]),0.0,lateral_flow[-1,1])*northeast+if(isnull(lateral_flow[0,1]),0.0,lateral_flow[0,1])*east+if(isnull(lateral_flow[1,1]),0.0,lateral_flow[1,1])*southeast+if(isnull(lateral_flow[1,0]),0.0,lateral_flow[1,0])*south+if(isnull(lateral_flow[1,-1]),0.0,lateral_flow[1,-1])*southwest + if(isnull(lateral_flow[0,-1]),0.0,lateral_flow[0,-1])*west +if(isnull(lateral_flow[-1,-1]),0.0,lateral_flow[-1,-1])*northwest)' --o`;

	# water distribution
	print `r.mapcalc 'storage_amt_A_tmp = if(storage_amt<fieldcap_amt,min(storage_amt_A,fieldcap_amt_A),storage_amt_A)' --o`;

	print `r.mapcalc 'storage_amt_B = if(storage_amt<fieldcap_amt,storage_amt-storage_amt_A_tmp,if((storage_amt-fieldcap_amt_A)<sat_mc_B*soil_depth_B,storage_amt-fieldcap_amt_A,sat_mc_B*soil_depth_B))' --o`;

	print `r.mapcalc 'storage_amt_A = if(storage_amt<fieldcap_amt,storage_amt_A_tmp,if((storage_amt-storage_amt_B)<sat_mc_A*soil_depth_A,storage_amt-storage_amt_B,sat_mc_A*soil_depth_A))' --o`;

	print `r.mapcalc 'storage_diff = storage_amt - storage_amt_A - storage_amt_B' --o`;
		
	#____________________________________________________________________________________
	#
	#    4. EVAPOTRANSPIRATION CALCULATION
	#____________________________________________________________________________________

	#    Actual ET is based on moisture content, ETcoeff, and potential ET.
	#    ET takes place at potential rate at ETreduction_mc or higher
	#    moisture contents, stops at wilting point or below, and is linearly
	#    related to moisture content between (Thornthwaite-Mather assumption)

	# 	Landuse classes:
	#	1. water: 0 cm
	#	2. urban/rock/barren/other: 0.1 cm
	#	3. forest/woody wetlands: .276 - 0.31 cm
	# 	4. shrub: 0.15 cm
	#	5. grass/grassy wetland/pasture: 0.15 cm
	#	6. row crop: 0.068 - 0.147 cm


	print `r.mapcalc 'ET_coeff = if(landuse==6.0,$cc_row_crop,if(landuse==5.0,$cc_grass,if(landuse==4.0,$cc_shrub,if(landuse==3.0,$cc_forest,if(landuse==2.0,$cc_urban,if(landuse==1.0,$cc_water,$cc_water))))))' --o`; # DJ - 3/14/23 --> taking out the -10 and removing the /100 below because our Kc values are small already

	print `r.mapcalc 'root_storage_amt = if(landuse==6.0,storage_amt,if(landuse==5.0,storage_amt,if(landuse==4.0,storage_amt,if(landuse==3.0,storage_amt_A,if(landuse==2.0,storage_amt,if(landuse==1.0,storage_amt,storage_amt))))))' --o`;

	print `r.mapcalc 'actualET_flow = if(root_storage_amt>(root_zone*ETreduction_mc),min(max(0.0,(root_storage_amt-wiltpt_amt)),pet*ET_coeff*$temp_time_step/24.0),if(root_storage_amt>wiltpt_amt,min(max(0.0,(root_storage_amt-wiltpt_amt)),max(0.0,pet*((root_storage_amt/root_zone-wiltpt_amt/root_zone)/(ETreduction_mc-wiltpt_amt/root_zone))*ET_coeff*$temp_time_step/24.0)),0.0))' --o`; ### FIX  to avoid negative values - MZ 20171030 Fixed # DJ - 3/14/23 --> removing

	print `r.mapcalc 'actualET_daily_flow = actualET_daily_flow+actualET_flow' --o`;
	print `r.mapcalc 'storage_amt = storage_amt-actualET_flow' --o`;

	print `r.mapcalc 'storage_amt_A = if(landuse==6.0,if(storage_amt<fieldcap_amt,fieldcap_amt_A*storage_amt/fieldcap_amt,fieldcap_amt_A),if(landuse==5.0,if(storage_amt<fieldcap_amt,fieldcap_amt_A*storage_amt/fieldcap_amt,fieldcap_amt_A),if(landuse==4.0,if(storage_amt<fieldcap_amt,fieldcap_amt_A*storage_amt/fieldcap_amt,fieldcap_amt_A),if(landuse==3.0,storage_amt_A-actualET_flow,if(landuse==2.0,if(storage_amt<fieldcap_amt,fieldcap_amt_A*storage_amt/fieldcap_amt,fieldcap_amt_A),if(landuse==1.0,if(storage_amt<fieldcap_amt,fieldcap_amt_A*storage_amt/fieldcap_amt,fieldcap_amt_A),0))))))' --o`;

	print `r.mapcalc 'storage_amt_B = storage_amt-storage_amt_A' --o`; # 20171207 MZ


	#____________________________________________________________________________________
	#
	#    5. PERCOLATION
	#____________________________________________________________________________________

	#    Moisture above fieldcap.mc can percolate to the restricting layer
	#    if the hydraulic conductivity of the restricting layer is adequate.
	#    The amount of percolation is substracted from the storage.

	print `r.mapcalc 'perc = min(min(max(storage_amt-fieldcap_amt,0.0),storage_amt_B),Ksubsurface*$temp_time_step/24.0)' --o`;
	print `r.mapcalc 'storage_amt = storage_amt - perc' --o`;
	print `r.mapcalc 'storage_amt_B = storage_amt_B - perc' --o`;

	# soil layer water distribution
	print `r.mapcalc 'storage_amt_A_tmp = if(storage_amt<fieldcap_amt,min(storage_amt_A,fieldcap_amt_A),storage_amt_A)' --o`;

	print `r.mapcalc 'storage_amt_B = if(storage_amt<fieldcap_amt,storage_amt-storage_amt_A_tmp,if((storage_amt-fieldcap_amt_A)<sat_mc_B*soil_depth_B,storage_amt-fieldcap_amt_A,sat_mc_B*soil_depth_B))' --o`;

	print `r.mapcalc 'storage_amt_A = if(storage_amt<fieldcap_amt,storage_amt_A_tmp,if((storage_amt-storage_amt_B)<sat_mc_A*soil_depth_A,storage_amt-storage_amt_B,sat_mc_A*soil_depth_A))' --o`;

	print `r.mapcalc 'perc_daily_flow = perc_daily_flow + perc' --o`;

	#____________________________________________________________________________________
	#
	#    6. RUNOFF
	#____________________________________________________________________________________

	#    Moisture above saturation in each cell becomes runoff.

	print `r.mapcalc 'runoff_flow = if(storage_amt>sat_amt,storage_amt-sat_amt,0.0)' --o`;

	print `r.mapcalc 'storage_amt = storage_amt-runoff_flow' --o`;
	print `r.mapcalc 'runoff_daily_flow = runoff_daily_flow+runoff_flow' --o`;

	print `r.mapcalc 'saturation = storage_amt/sat_amt*100.0' --o`;
	print `r.mapcalc 'runoff_$year = runoff_flow+runoff_$year' --o`;

	}; # end of the 6 hours loop

	#print "\n \n";
	#print "\n|----- END OF 6 HOUR TIME LOOPING -----|\n";
	#print "\n \n";

	print `r.mapcalc 'runoff_daily_flow = runoff_daily_flow+road_runoff' --o`;
	print `r.mapcalc 'input_daily_balance = input_daily - (snowmelt + throughfall)' --o`;

	print `r.mapcalc 'Psat_$year = Psat_$year + if(saturation>=100,1,0)' --o`;


#____________________________________________________________________________________
#
#  7. Output maps
#____________________________________________________________________________________

	# the idea here is too sum up the daily runoff if it's february (month == 2) and then 
	# if the day is 60, which will be right at the end of february write out the file to 
	# an ascii (it will be hard to say if this works until we start trying to run the model)
	if ($month == 2) {
		`r.mapcalc 'runoff_feb_$year = runoff_feb_$year + runoff_daily_flow' --o`;
		}

	# this outputs accumulated and average february runoff but doesn't account for shortened 
	# february --> this may need to be changed.
	if ($doy == 61) {
		print `r.mapcalc 'average_feb_runoff_$year = runoff_feb_$year / 29' --o`;
		print `r.out.ascii input=average_feb_runoff output=feb_avg_runoff_$year`;
		}

	#print "\n \n";
	#print "\n|----- OUTPUT ZONAL STATS -----|\n";
	#print "\n \n";


	open($WATERSHED, '<', 'wshed_res_properties.ini') || "Can't open file\n";
	while (<$WATERSHED>) {
		chop($_);
		($wshed_id,$area_cells,$res_vol,$res_coeff) = split;
			#print `r.mapcalc 'MASK = if(wsheds_all==$wshed_id,1,0)' --o`;
		
		# **********************************  1  ***************************************
		# Runoff output
			print `r.stats.zonal base=watershed cover=runoff_daily_flow out=temp1 method=sum --o --quiet`;
			print $runoff_cm_{$wshed_id} = `r.stats -A -n -N input=temp1`;
			$runoff_cm_{$wshed_id} = $runoff_cm_{$wshed_id}*1; # sum of all cells in cm
			$runoff_cms_{$wshed_id} = ($runoff_cm_{$wshed_id}/100)*($gridsize*$gridsize)/($time_step*3600.0); # Calculate the runoff in m3/s (cubic meter per second) for sum of all cells
			$runoff_mm_{$wshed_id} = $runoff_cm_{$wshed_id}*10/$area_cells; # Convert runoff into area average depth (mm)


		# **********************************  3  ***************************************
		# Precip output
			print `r.stats.zonal base=watershed cover=precip out=temp3 method=sum --o --quiet`;
			print $precip_cm_{$wshed_id} = `r.stats -A -n -N input=temp3`;
			$precip_cm_{$wshed_id} = $precip_cm_{$wshed_id}*1;


		# **********************************  4  ***************************************
		# Rain output
			print `r.stats.zonal base=watershed cover=rain out=temp4 method=sum --o --quiet`;
			print $rain_cm_{$wshed_id} = `r.stats -A input=temp4 nv= `;
			$rain_cm_{$wshed_id} = $rain_cm_{$wshed_id}*1;


		# **********************************  5  ***************************************
		# Actual evaporation output
			print `r.stats.zonal base=watershed cover=actualET_daily_flow out=temp5 method=sum --o --quiet`;
			print $actualET_flow_cm_{$wshed_id} = `r.stats -A input=temp5 nv= `;
			$actualET_flow_cm_{$wshed_id} = $actualET_flow_cm_{$wshed_id}*1;


		# **********************************  6  ***************************************
		# Canopy ET output
			print `r.stats.zonal base=watershed cover=canopy_evap out=temp6 method=sum --o --quiet`;
			print $canopy_evap_cm_{$wshed_id} = `r.stats -A input=temp6 nv= `;
			$canopy_evap_cm_{$wshed_id} = $canopy_evap_cm_{$wshed_id}*1;


		# **********************************  7  ***************************************
		# Snowmelt output
			print `r.stats.zonal base=watershed cover=snowmelt out=temp7 method=sum --o --quiet`;
			print $snowmelt_cm_{$wshed_id} = `r.stats -A input=temp7 nv= `;
			$snowmelt_cm_{$wshed_id} = $snowmelt_cm_{$wshed_id}*1;


		# **********************************  8  ***************************************
		# Storage amount output
			print `r.stats.zonal base=watershed cover=storage_amt out=temp8 method=sum --o --quiet`;
			print $storage_amt_cm_{$wshed_id} = `r.stats -A input=temp8 nv= `;
			$storage_amt_cm_{$wshed_id} = $storage_amt_cm_{$wshed_id}*1;


		# **********************************  9  ***************************************
		# Throughfall output
			print `r.stats.zonal base=watershed cover=throughfall out=temp9 method=sum --o --quiet`;
			print $throughfall_cm_{$wshed_id} = `r.stats -A input=temp9 nv= `;
			$throughfall_cm_{$wshed_id} = $throughfall_cm_{$wshed_id}*1;		


		# **********************************  10  ***************************************
		# Canopy Storage Amount output
			print `r.stats.zonal base=watershed cover=canopy_storage_amt out=temp10 method=sum --o --quiet`;
			print $canopy_storage_amt_cm_{$wshed_id} = `r.stats -A input=temp10 nv= `;
			$canopy_storage_amt_cm_{$wshed_id} = $canopy_storage_amt_cm_{$wshed_id}*1;


		# **********************************  11  ***************************************
		# Percolation output
			print `r.stats.zonal base=watershed cover=perc_daily_flow out=temp11 method=sum --o --quiet`;
			print $perc_cm_{$wshed_id} = `r.stats -A input=temp11 nv= `;
			$perc_cm_{$wshed_id} = $perc_cm_{$wshed_id}*1;
			$perc_cms_{$wshed_id} = $perc_cm_{$wshed_id}/100.0/($time_step*3600.0); # Calculate the percolation in m3/s for sum of all cells # Duncan Jurayj - 3/21/23: removed this: "*$gridsize*$gridsize"
			$perc_mm_{$wshed_id} = $perc_cm_{$wshed_id}*10/$area_cells; # Convert percolation into area average depth (mm)

		# **********************************  12  ***************************************
		# Streamflow calculations
			$res_vol_{$wshed_id} = $res_vol_{$wshed_id} + $perc_cms_{$wshed_id} - $base_flow_{$wshed_id};
			$base_flow_{$wshed_id} = $res_coeff_{$wshed_id} * $res_vol_{$wshed_id};
			$Q_{$wshed_id}=$base_flow_{$wshed_id}+$runoff_cms_{$wshed_id}; # m3/s for sum of all cells
			$Q_mm_{$wshed_id} = $Q_{$wshed_id}*86400*1000/($area_cells*$gridsize*$gridsize); #convert Q into average depth in mm

		# **********************************  17 swe temporary test  ***************************************
		# Snow Water Equivalent (average) output
			print `r.stats.zonal base=watershed cover=swe out=temp17 method=sum --o --quiet`;
			print $swe_cm_{$wshed_id} = `r.stats -A input=temp17 nv= `;
			$swe_cm_{$wshed_id} = $swe_cm_{$wshed_id}*1;

		# **********************************  22 condens temporary test  ***************************************
		# Condens (average) output
			print `r.stats.zonal base=watershed cover=condens out=temp22 method=sum --o --quiet`;
			print $condens_cm_{$wshed_id} = `r.stats -A input=temp22 nv= `;
			$condens_cm_{$wshed_id} = $condens_cm_{$wshed_id}*1;

		# **********************************  23 snow temporary test  ***************************************
		# Snowfall (average) output
			print `r.stats.zonal base=watershed cover=snow out=temp23 method=sum --o --quiet`;
			print $snow_cm_{$wshed_id} = `r.stats -A input=temp23`;
			$snow_cm_{$wshed_id} = $snow_cm_{$wshed_id}*1;

				# **********************************  24  ***************************************
		# q.srad output
			print `r.stats.zonal base=watershed cover=q.srad out=temp24 method=sum --o --quiet`;
			print $srad_{$wshed_id} = `r.stats -A -n -N input=temp24`;
			$srad_{$wshed_id} = $srad_{$wshed_id}*1;

		# **********************************  25  ***************************************
		# q.latent output
			print `r.stats.zonal base=watershed cover=q.latent out=temp25 method=sum --o --quiet`;
			print $latent_{$wshed_id} = `r.stats -A -n -N input=temp25`;
			$latent_{$wshed_id} = $latent_{$wshed_id}*1;


		# **********************************  26  ***************************************
		# q.sensible output
			print `r.stats.zonal base=watershed cover=q.sensible out=temp26 method=sum --o --quiet`;
			print $sensible_{$wshed_id} = `r.stats -A input=temp26`;
			$sensible_{$wshed_id} = $sensible_{$wshed_id}*1;


		# **********************************  27  ***************************************
		# q.lw output
			print `r.stats.zonal base=watershed cover=q.lw out=temp27 method=sum --o --quiet`;
			print $lw_{$wshed_id} = `r.stats -A input=temp27`;
			$lw_{$wshed_id} = $lw_{$wshed_id}*1;


		# **********************************  28  ***************************************
		# q.rain.ground output
			print `r.stats.zonal base=watershed cover=q.rain.ground out=temp28 method=sum --o --quiet`;
			print $q_rain_ground_cm_{$wshed_id} = `r.stats -A input=temp28`;
			$q_rain_ground_cm_{$wshed_id} = $q_rain_ground_{$wshed_id} * 1;

		# **********************************  29 ***************************************
		# q.total output
			print `r.stats.zonal base=watershed cover=q.total out=temp29 method=sum --o --quiet`;
			print $q_total_{$wshed_id} = `r.stats -A input=temp29`;
			$q_total_{$wshed_id} = $q_total_{$wshed_id}*1;

		#  Create mass balance output file
		open(OUT, ">>/Users/duncanjurayj/Documents/SMR_R/raw_data/smr_output/MFC_mass_balance_$wshed_id.csv") || die("Cannot Open File");
		print OUT "$wshed_id $date $year $runoff_cm_{$wshed_id} $precip_cm_{$wshed_id} $rain_cm_{$wshed_id} $actualET_flow_cm_{$wshed_id} $canopy_evap_cm_{$wshed_id} $snowmelt_cm_{$wshed_id} $storage_amt_cm_{$wshed_id} $throughfall_cm_{$wshed_id} $canopy_storage_amt_cm_{$wshed_id} $perc_cm_{$wshed_id} $Q_{$wshed_id} $swe_cm_{$wshed_id} $condens_cm_{$wshed_id} $snow_cm_{$wshed_id} $base_flow_{$wshed_id} $srad_{$wshed_id} $latent_{$wshed_id} $sensible_{$wshed_id} $lw_{$wshed_id} $q_rain_ground_cm_{$wshed_id} $q_total_{$wshed_id}\n";
		close(OUT); 


		#print "\n \n";
		#print "\n|----- END OF OUTPUTS ZONAL STATS -----|\n";
		#print "\n \n";
		
	}
	close ($WATERSHED);


}
close ($WEATHER);

#print "\n \n";
#print "\n|---------- SMR-SAM HAS FINISHED RUNNING ----------|\n";
#print "\n \n";

