#!/usr/bin/perl
# this is a test smr perl script where outputs are dependent on the user --> 
#   this method should allow multiple users on different - smr_multiple_users.pl

print `g.remove -f type=raster name=MASK`;
print `g.region rast=watershed`;
print `r.mapcalc 'watershed = if(watershed==1,1,0)' --o`;

# setting user specific variables: user name, operating system, and home 
# directory home directory name

my $operating_system = $^O;
my $home_dir_name;
my $user;

if ($operating_system eq 'msys') {
    $user = $ENV{'USERNAME'};
} elsif ($operating_system eq 'darwin') {
    $user = $ENV{'USER'};
} elsif ($operating_system eq 'linux') { # this may also need to be changed
    $user = $ENV{'USER'}; # this may also need to be changed 
}

if ($operating_system eq 'darwin') {
    $home_dir_name = "Users";
} elsif ($operating_system eq 'linux') {
    $home_dir_name = "home";
} elsif ($operating_system eq 'msys') {
    $home_dir_name = "C:/Users"; # this may need to become #C:/Users -- 
}

print "The current user is: $user\n";
print "The current operating system is: $operating_system\n";
print "Based on the operating system the home directory is: $home_dir_name\n";
print "Weather file directory is: /$home_dir_name/$user/Dropbox/SMR_R/raw_data/weather/noaa_pullman_mini_sn.csv\n";

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
#$mean_annual_sea_level_precip = 15.53944; # calculated based on pullman average and lapse rate
$precip_LR = 0.08; # cm / m based on 2000-2020 mean annual precip diff betwn pullman & Moscow Mountain
$temp_LR = -2.20; # C / km based on 2000-2020 mean annual temp diff betwn pullman & Moscow Mountain

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

# Read each line of the watershed file and assign values based on watershed ID
print "\n \n";
print "\n|----- READING WATERSHED PROPERTIES -----|\n";
print "\n \n";
while (my $line = <$WATERSHED>) {
	chop($line);
	($wshed_id,$area_cells,$res_vol,$res_coeff) = split(' ', $line);
	
	if ($wshed_id > 0.0) {
		$area_{$wshed_id} = $area_cells * $gridsize * $gridsize; #square meters
		$res_vol_{$wshed_id} = 10;
		$res_coeff_{$wshed_id} = $res_coeff;
		$base_flow_{$wshed_id} = 0.0;
		
		print "	|---- Properties of watershed with ID: $wshed_id ----|\n";
		print "		number of cells in watershed: $area_cells\n";
		print "		cell size (m): $gridsize\n";
		print "		total watershed area (m^2): $area_{$wshed_id}\n";
        print "		res_vol: $res_vol_{$wshed_id}\n";
        print "		res_coeff: $res_coeff_{$wshed_id}\n";
        print "		base_flow: $base_flow_{$wshed_id}\n\n";
		print "	|-------------------------------------------------|\n";
	};
}

# Close the watershed properties file
close($WATERSHED) || die "Cannot close the watershed properties file";

#____________________________________________________________________________________
#
# PERMANENT MAPS REFERENCED FOR MODEL CALCULATIONS
#____________________________________________________________________________________

#         el (The DEM in meters)
#         watershed (A map defining the extent of the watershed, 1=inside watershed)
#         fieldcap_mc =  field capacity moisture content (m^3/m^3), water drains at moisture
#           contents above field capacity
#         soil_depth =  Depth to a hydraulic restrictive layer (cm), includes all soil layers
#         wiltpt_amt =  Wilting point moisture depth (cm)
#         ETreduction_mc = Moisture content at which actual ET becomes limited by soil moisture
#					  (assumed to be 80% of field capacity) (m^3/m^3)
#         Ksat_matix =  The lateral matrix Ksat (ignoring macropores) of the soil layer (cm/day)
#	        Ksat_mpore =  The lateral Ksat of the soil layer when macropores are active (cm/day)
#         Ksubsurface = The vertical Ksat through the hydraulic restrictive layer (cm/day)

#____________________________________________________________________________________
#
# MODEL INITIAL CONDITIONS & INTIAL MAPS
#____________________________________________________________________________________
# land use map legend: 1 = water, 2 = urban/rock/barren/other, 3 = forest/woody wetlands, 4 = shrub, 5 = grass/grassy wetland/pasture, 6 = row crop
#	maximum canopy storage (cm) in descending land use number (6-1) order -- can range: 
	#	1. water: 0 cm
	#	2. urban/rock/barren/other: 0.1 cm (Li, X. et al. 2005) 
	#	3. forest/woody wetlands: 0.276 - 0.31 cm (Soto‐Schönherr, S., & Iroumé, A. 2016), spruce (Floriancic et al. 2022)
	# 4. shrub: 0.15 cm (Tromble, J. M. 1988)
	#	5. grass/grassy wetland/pasture: 0.15 cm (Gifford, G. F. 1976)
	#	6. row crop: 0.068 - 0.147 cm (Di, W., Jiusheng, L., & Minjie, R. 2006)
print `r.mapcalc 'max_canopy_storage_amt = if(landuse==6.0,0.147,if(landuse==5.0,0.15,if(landuse==4.0,0.15,if(landuse==3.0,0.31,if(landuse==2.0,0.1,if(landuse==1.0,0.0,0.0))))))' --o`;

# initial water stored in the canopies is zero
print `r.mapcalc 'canopy_storage_amt = 0.0' --o`;

# Canopy cover is defined by landuse and is used to modify the calculation of resistence to heat transfer in the SAM Module; values range from 0.1-0.99; larger values increase resistence which slows snowmelt. 
print `r.mapcalc 'canopy_cover = if(landuse==6.0,0.1,if(landuse==5.0,0.2,if(landuse==4.0,0.50,if(landuse==3.0,0.99,if(landuse==2.0,0.01,if(landuse==1.0,0.001,0.001))))))' --o`;

# Root zone defines the portion of the soil that is affected by actual evaporation.
print `r.mapcalc 'root_zone = if(landuse==6.0,soil_depth,if(landuse==5.0,soil_depth_A,if(landuse==4.0,soil_depth,if(landuse==3.0,soil_depth,if(landuse==2.0,1.0,if(landuse==1.0,1.0,1.0))))))' --o`;

# Wiltpt amount is the amount of water that causes plants to wilt. This is a base map already made in the setup.
#print `r.mapcalc 'wiltpt_amt = if(landuse==6.0,wiltpt_mc_A*soil_depth_A+wiltpt_mc_B*soil_depth_B,if(landuse==5.0,wiltpt_mc_A*soil_depth_A,if(landuse==4.0,wiltpt_mc_A*soil_depth_A+wiltpt_mc_B*soil_depth_B,if(landuse==3.0,wiltpt_mc_A*soil_depth_A+wiltpt_mc_B*soil_depth_B,if(landuse==2.0,wiltpt_mc_A*soil_depth_A+wiltpt_mc_B*soil_depth_B,if(landuse==1.0,wiltpt_mc_A*soil_depth_A+wiltpt_mc_B*soil_depth_B,wiltpt_mc_A*soil_depth_A+wiltpt_mc_B*soil_depth_B))))))' --o`;

# ETreduction is the amount of water that can evaporate from the soil; 80% for all classes except forest which is 70% and urban/barren/rock/water is 0. 
print `r.mapcalc 'ETreduction_mc = if(landuse==6.0,fieldcap_amt_A*0.8/soil_depth_A,if(landuse==5.0,fieldcap_amt_A*0.8/soil_depth_A,if(landuse==4.0,fieldcap_amt_A*0.8/soil_depth_A,if(landuse==3.0,fieldcap_amt_A*0.8/soil_depth_A,if(landuse==2.0,fieldcap_amt_A*0.8/soil_depth_A,if(landuse==1.0,fieldcap_amt_A*0.8/soil_depth_A,fieldcap_amt_A*0.8/soil_depth_A))))))' --o`;

# The temperature above which all precip is considered to fall as rain.
print `r.mapcalc 'tmax_rain = if(landuse==6.0,1.0,if(landuse==5.0,1.0,if(landuse==4.0,1.0,if(landuse==3.0,1.0,if(landuse==2.0,1.0,if(landuse==1.0,1.0,0))))))' --o`;

# The temperature below which all precip is considered to fall as snow.
print `r.mapcalc 'tmin_snow = if(landuse==6.0,0.0,if(landuse==5.0,0.0,if(landuse==4.0,0.0,if(landuse==3.0,0.0,if(landuse==2.0,0.0,if(landuse==1.0,0.0,0.0))))))' --o`;

# Initial conditions storage amount defined (% of saturation)
#  in setup - storage_amt in stream cells is assumed at saturation
#  in setp up - soil depth in stream cells is assumed at 1 cm A horizon, 19 cm B horizon
print `r.mapcalc 'storage_amt_A = sat_mc_A * 0.4' --o`;
print `r.mapcalc 'storage_amt_B = sat_mc_B * 0.4' --o`;
#print `r.mapcalc 'storage_amt_A = (wiltpt_mc_A+(fieldcap_mc_A-wiltpt_mc_A)*0.2)*soil_depth_A' --o`;
print `r.mapcalc 'storage_amt = storage_amt_A + storage_amt_B' --o`;

#  Initiation for SAM
print `r.mapcalc 'swe = 0.0' --o`;
print `r.mapcalc 'snow.age = 1.0' --o`;
print `r.mapcalc 'swe.yesterday = 0.0' --o`;
print `r.mapcalc 'albedo = 0.2' --o`;
print `r.mapcalc 'liquid.water = 0.0' --o`;
print `r.mapcalc 'ice.content = 0.0' --o`;
print `r.mapcalc 'tsnow_surf = -2.0' --o`;
print `r.mapcalc 'u.surface = 24.0*tsnow_surf*(2.1*1000.0*min(swe/100.0,0.02)+1000.0*2.1*max(0.0,0.02-swe/100.0))' --o`; # convert to daily MZ 20190410 | should not be convert to daily becasue u.surface is in KJ/m^2. It is not a rate, it is the energy content of the surface layer. So the comment on 20190410 is wrong. MZ 20200330
print `r.mapcalc 'tsnow.pack = -0.5' --o`;
print `r.mapcalc 'u.total = 24.0*tsnow.pack*(2.1*1000.0*swe/100.0+2.1*1000.0*0.4)' --o`; # convert to daily MZ 20190410 | should not be convert to daily becasue u.surface is in KJ/m^2. It is not a rate, it is the energy content of the surface layer. So the comment on 20190410 is wrong. MZ 20200330

#  set the initial (t-1) time step
print `r.mapcalc 'mass_balance_total = 0.0' --o`;
print `r.mapcalc 'canopy_storage_amt_pre = 0.0' --o`;
print `r.mapcalc 'storage_amt_pre = storage_amt' --o`;

# initial maps for each year of weather
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

print `r.mapcalc 'A_amt_feb_1965 = 0.0' --o`;
print `r.mapcalc 'A_amt_feb_1966 = 0.0' --o`;
print `r.mapcalc 'A_amt_feb_1967 = 0.0' --o`;
print `r.mapcalc 'A_amt_feb_1968 = 0.0' --o`;
print `r.mapcalc 'A_amt_feb_1969 = 0.0' --o`;
print `r.mapcalc 'A_amt_feb_1970 = 0.0' --o`;
print `r.mapcalc 'A_amt_feb_1971 = 0.0' --o`;
print `r.mapcalc 'A_amt_feb_1972 = 0.0' --o`;
print `r.mapcalc 'A_amt_feb_1973 = 0.0' --o`;
print `r.mapcalc 'A_amt_feb_1974 = 0.0' --o`;
print `r.mapcalc 'A_amt_feb_1975 = 0.0' --o`;
print `r.mapcalc 'A_amt_feb_1976 = 0.0' --o`;
print `r.mapcalc 'A_amt_feb_1977 = 0.0' --o`;
print `r.mapcalc 'A_amt_feb_1978 = 0.0' --o`;
print `r.mapcalc 'A_amt_feb_1979 = 0.0' --o`;
print `r.mapcalc 'A_amt_feb_1980 = 0.0' --o`;
print `r.mapcalc 'A_amt_feb_1981 = 0.0' --o`;
print `r.mapcalc 'A_amt_feb_1982 = 0.0' --o`;
print `r.mapcalc 'A_amt_feb_1983 = 0.0' --o`;
print `r.mapcalc 'A_amt_feb_1984 = 0.0' --o`;
print `r.mapcalc 'A_amt_feb_1985 = 0.0' --o`;
print `r.mapcalc 'A_amt_feb_1986 = 0.0' --o`;
print `r.mapcalc 'A_amt_feb_1987 = 0.0' --o`;
print `r.mapcalc 'A_amt_feb_1988 = 0.0' --o`;

print `r.mapcalc 'B_amt_feb_1965 = 0.0' --o`;
print `r.mapcalc 'B_amt_feb_1966 = 0.0' --o`;
print `r.mapcalc 'B_amt_feb_1967 = 0.0' --o`;
print `r.mapcalc 'B_amt_feb_1968 = 0.0' --o`;
print `r.mapcalc 'B_amt_feb_1969 = 0.0' --o`;
print `r.mapcalc 'B_amt_feb_1970 = 0.0' --o`;
print `r.mapcalc 'B_amt_feb_1971 = 0.0' --o`;
print `r.mapcalc 'B_amt_feb_1972 = 0.0' --o`;
print `r.mapcalc 'B_amt_feb_1973 = 0.0' --o`;
print `r.mapcalc 'B_amt_feb_1974 = 0.0' --o`;
print `r.mapcalc 'B_amt_feb_1975 = 0.0' --o`;
print `r.mapcalc 'B_amt_feb_1976 = 0.0' --o`;
print `r.mapcalc 'B_amt_feb_1977 = 0.0' --o`;
print `r.mapcalc 'B_amt_feb_1978 = 0.0' --o`;
print `r.mapcalc 'B_amt_feb_1979 = 0.0' --o`;
print `r.mapcalc 'B_amt_feb_1980 = 0.0' --o`;
print `r.mapcalc 'B_amt_feb_1981 = 0.0' --o`;
print `r.mapcalc 'B_amt_feb_1982 = 0.0' --o`;
print `r.mapcalc 'B_amt_feb_1983 = 0.0' --o`;
print `r.mapcalc 'B_amt_feb_1984 = 0.0' --o`;
print `r.mapcalc 'B_amt_feb_1985 = 0.0' --o`;
print `r.mapcalc 'B_amt_feb_1986 = 0.0' --o`;
print `r.mapcalc 'B_amt_feb_1987 = 0.0' --o`;
print `r.mapcalc 'B_amt_feb_1988 = 0.0' --o`;

print `r.mapcalc 'precip_1965 = 0.0' --o`;
print `r.mapcalc 'precip_1966 = 0.0' --o`;
print `r.mapcalc 'precip_1967 = 0.0' --o`;
print `r.mapcalc 'precip_1968 = 0.0' --o`;
print `r.mapcalc 'precip_1969 = 0.0' --o`;
print `r.mapcalc 'precip_1970 = 0.0' --o`;
print `r.mapcalc 'precip_1971 = 0.0' --o`;
print `r.mapcalc 'precip_1972 = 0.0' --o`;
print `r.mapcalc 'precip_1973 = 0.0' --o`;
print `r.mapcalc 'precip_1974 = 0.0' --o`;
print `r.mapcalc 'precip_1975 = 0.0' --o`;
print `r.mapcalc 'precip_1976 = 0.0' --o`;
print `r.mapcalc 'precip_1977 = 0.0' --o`;
print `r.mapcalc 'precip_1978 = 0.0' --o`;
print `r.mapcalc 'precip_1979 = 0.0' --o`;
print `r.mapcalc 'precip_1980 = 0.0' --o`;
print `r.mapcalc 'precip_1981 = 0.0' --o`;
print `r.mapcalc 'precip_1982 = 0.0' --o`;
print `r.mapcalc 'precip_1983 = 0.0' --o`;
print `r.mapcalc 'precip_1984 = 0.0' --o`;
print `r.mapcalc 'precip_1985 = 0.0' --o`;
print `r.mapcalc 'precip_1986 = 0.0' --o`;
print `r.mapcalc 'precip_1987 = 0.0' --o`;
print `r.mapcalc 'precip_1988 = 0.0' --o`;

print `r.mapcalc 'pet_1965 = 0.0' --o`;
print `r.mapcalc 'pet_1966 = 0.0' --o`;
print `r.mapcalc 'pet_1967 = 0.0' --o`;
print `r.mapcalc 'pet_1968 = 0.0' --o`;
print `r.mapcalc 'pet_1969 = 0.0' --o`;
print `r.mapcalc 'pet_1970 = 0.0' --o`;
print `r.mapcalc 'pet_1971 = 0.0' --o`;
print `r.mapcalc 'pet_1972 = 0.0' --o`;
print `r.mapcalc 'pet_1973 = 0.0' --o`;
print `r.mapcalc 'pet_1974 = 0.0' --o`;
print `r.mapcalc 'pet_1975 = 0.0' --o`;
print `r.mapcalc 'pet_1976 = 0.0' --o`;
print `r.mapcalc 'pet_1977 = 0.0' --o`;
print `r.mapcalc 'pet_1978 = 0.0' --o`;
print `r.mapcalc 'pet_1979 = 0.0' --o`;
print `r.mapcalc 'pet_1980 = 0.0' --o`;
print `r.mapcalc 'pet_1981 = 0.0' --o`;
print `r.mapcalc 'pet_1982 = 0.0' --o`;
print `r.mapcalc 'pet_1983 = 0.0' --o`;
print `r.mapcalc 'pet_1984 = 0.0' --o`;
print `r.mapcalc 'pet_1985 = 0.0' --o`;
print `r.mapcalc 'pet_1986 = 0.0' --o`;
print `r.mapcalc 'pet_1987 = 0.0' --o`;
print `r.mapcalc 'pet_1988 = 0.0' --o`;

print `r.mapcalc 'psat_1965 = 0.0' --o`;
print `r.mapcalc 'psat_1966 = 0.0' --o`;
print `r.mapcalc 'psat_1967 = 0.0' --o`;
print `r.mapcalc 'psat_1968 = 0.0' --o`;
print `r.mapcalc 'psat_1969 = 0.0' --o`;
print `r.mapcalc 'psat_1970 = 0.0' --o`;
print `r.mapcalc 'psat_1971 = 0.0' --o`;
print `r.mapcalc 'psat_1972 = 0.0' --o`;
print `r.mapcalc 'psat_1973 = 0.0' --o`;
print `r.mapcalc 'psat_1974 = 0.0' --o`;
print `r.mapcalc 'psat_1975 = 0.0' --o`;
print `r.mapcalc 'psat_1976 = 0.0' --o`;
print `r.mapcalc 'psat_1977 = 0.0' --o`;
print `r.mapcalc 'psat_1978 = 0.0' --o`;
print `r.mapcalc 'psat_1979 = 0.0' --o`;
print `r.mapcalc 'psat_1980 = 0.0' --o`;
print `r.mapcalc 'psat_1981 = 0.0' --o`;
print `r.mapcalc 'psat_1982 = 0.0' --o`;
print `r.mapcalc 'psat_1983 = 0.0' --o`;
print `r.mapcalc 'psat_1984 = 0.0' --o`;
print `r.mapcalc 'psat_1985 = 0.0' --o`;
print `r.mapcalc 'psat_1986 = 0.0' --o`;
print `r.mapcalc 'psat_1987 = 0.0' --o`;
print `r.mapcalc 'psat_1988 = 0.0' --o`;

#____________________________________________________________________________________
#
# START READING WEATHER DATA
#____________________________________________________________________________________

print "\n|----- READING Weather -----|\n";

#  This set of commands splits a tab delimited array using a while loop; it runs until line 1013
open ($WEATHER, '<', "/$home_dir_name/$user/Dropbox/SMR_R/raw_data/weather/noaa_pullman_mini_sn.csv") || die "Can't open weather file\n"; 

while (<$WEATHER>) {
	chop($_);
	($date,$year,$month,$day,$doy,$tmax,$tmin,$tavg,$tdew,$precip,$pet,$hour_1,$hour_6,$hour_12,$hour_18,$l_turb,$cloud,$cc_water,$cc_urban,$cc_forest,$cc_shrub,$cc_grass,$cc_row_crop,$rh_snow) = split(' ',$_);

	print "\n \n";
	print "\n|----- DAILY WEATHER READ FOR $date -----|\n";
	print "\n \n";
	print "\n	WEATHER SUMMARY:\n";
	print "\n	| Date $date | Year = $year | Julian Day = $doy | Tair = $tavg °C | Precip = $precip cm | PET = $pet cm | Tdew = $tdew";
	print "\n	|-------------------------------------------------------------------------------------------------------------------------| \n";

	@hourly_tmp_array = ($hour_1,$hour_6,$hour_12,$hour_18);

#  --------------------------------------------------------------------
#  --------------------------------------------------------------------
#____________________________________________________________________________________
#
#  1.  LAPSE RATE CALCS & PRECIPITATION PARTITIONING
#____________________________________________________________________________________
print "\n|----- Partitioning Precipitation for $date -----|\n";

# rain and snow are in (cm)
# Moscow Mountain SNOTEL elevation = 1443 m, Pullman Station elevation = 784 m
# precipitation lapse rate between sites = 0.08 cm / meter of elevation
# temperature lapse rate between sites = - -3.20 C/km --> also on the lower end it seems, but fine to use for now we can always modify both values --> divide by 1000 to work with meters
# precipitation throughout watershed is calculated with: precip_at_pullman * (elevation * lapse_rate / MAPpullman); the LR is from MAP differences btwn elevations of 784 and 1443 m; it predicts a MAP for a cell based on its elevation and divides that by the MAP at the precip measurement elevation.
# temperature throughout watershed is calculated with: temp - (temp_lapse_rate * (el - low_site_elevation)) --> this assumes data is taken from the lower site (Pullman site)
# an equation between PET and elevation.  PET_moscow_mountain = 1.8281 * PET_pullman + -0.2486 (cm)

print `r.mapcalc 'precip = $precip * (el * $precip_LR / $mean_annual_onsite_precip)' --o`;
print `r.mapcalc 'tavg = $tavg + ($temp_LR/1000) * (el - $low_site_el)' --o`;
print `r.mapcalc 'tdew = $tdew + ($temp_LR/1000) * (el - $low_site_el)' --o`;
print `r.mapcalc 'rain = if(tavg>tmax_rain,precip,if(tavg<tmin_snow,0.0,(tavg-tmin_snow)/(tmax_rain-tmin_snow)*precip))' --o`;
print `r.mapcalc 'snow = precip-rain' --o`;
print `r.mapcalc 'pet_data = (1.8281 * $pet - 0.2486)' --o`;
print `r.mapcalc "pet = (pet_data - $pet) * ($low_site_el - el)/($high_site_el - $low_site_el)+pet_data" --o`;

# interception is calculated for Spruce trees based on work by Lankreijer et al. 1999, Agricultural and Forest
# Meteorology 98-99:595.  Max Storage of Canopy was taken as 2.0 mm.  During rain, evaporation is 0.04 mm/hr.
# When there is no rain, for the time being, it is assumed that evaporation = 50% of PET.
# Assume evaporation of canopy storage is at the potential rate

print `r.mapcalc 'canopy_storage_amt = canopy_storage_amt + rain' --o`;
print `r.mapcalc 'throughfall = if(canopy_storage_amt>max_canopy_storage_amt,canopy_storage_amt-max_canopy_storage_amt,0.0)' --o`;
print `r.mapcalc 'canopy_storage_amt = canopy_storage_amt-throughfall' --o`;
print `r.mapcalc 'canopy_evap = if(rain>0.0,min(canopy_storage_amt,pet),min(canopy_storage_amt,pet))' --o`;
print `r.mapcalc 'pet = max(0.0,pet-canopy_evap)' --o`;
print `r.mapcalc 'canopy_storage_amt = canopy_storage_amt-canopy_evap' --o`;

#____________________________________________________________________________________
#
#  2.  SNOWMELT AND INTERCEPTION ALGORITHM
#____________________________________________________________________________________
print "\n|----- SNOW ACCUMULATION AND MELT MODEL -----|\n";
#  ------------------------------- SAM --------------------------------------
#  Snow accumulation and melt (SAM) model
#  Modified by MZ on 20190126
#  albedo approximated from DHSVM relationship
#  assumes the albedo of the soil is 0.2
#  albedo (100 x %)

print `r.mapcalc 'rh = $rh_snow/(1.0-(canopy_cover/3.0))' --o`; # duncan jurayj 2023051 ## testing for why soil storage is 0
print `r.mapcalc 'snow.age = if(snow>0.0 && throughfall==0.0,1.0,snow.age+1.0)' --o`; # modified MZ 20190210
print `r.mapcalc 'albedo = if(swe.yesterday+snow>0.0,min(0.95,0.7383*snow.age^(-0.1908)),0.2)' --o`;

#  Calculate the diffuse transmissivity following Bristow and Campbell (1985)
#  The clear sky transmissivity for Troy, ID is 0.75
print `r.mapcalc 't_diff = if($cloud==1.0,0.1,0.75*(1-$cloud)*(1-exp(-0.6*$cloud/((0.75-0.4)*(1-$cloud)))))' --o`; #* watershed' --o`; # removed *watershed 

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
#print `r.mapcalc 'q.srad = (1.0-canopy_cover)*(1.0-albedo)*(beam_rad+diff_rad+refl_rad) * 3600.0/1000.0' --o`; # add canopy cover impact MZ 20200510
print `r.mapcalc 'q.srad = (1.0-albedo)*(beam_rad+diff_rad+refl_rad) * 3600.0/1000.0' --o`; # remove canopy cover plm 20230426

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
print `r.mapcalc 'q.lw = 5.67*(10.0^-8.0)*3600.0*24.0/1000.0*(max(0.72,(9.2*(10^-6.0)*(tavg+273.15)^2.0)*(1.0-0.84*$cloud)+0.84*$cloud)*(tavg+273.15)^4.0-0.98*(tsnow_surf+273.15)^4.0)' --o`; # MZ 20190408 from Campbell an introductino to environmental biophysics eqn (10.11) to calculate clear sky emissivity e(0)=9.2*10^-6*Tavg^2, where Tavg is in Kelvin
#print `r.mapcalc 'q.lw = 5.67*(10.0^-8.0)*3600.0*24.0/1000.0*(max(0.72,((9.2*(10^-6.0)*(tavg+273.15)^2.0)*(1.0-0.84*$cloud)+0.84*$cloud)*(1-canopy_cover)+0.92*canopy_cover)*(tavg+273.15)^4.0-0.98*(tsnow_surf+273.15)^4.0)' --o`; # add canopy cover impact MZ 20200510

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

print `r.mapcalc 'q.latent = (2500.0+335.0)*(vap.d.air-vap.d.snow)/(rh)' --o`; # Duncan Jurayj - divided by 1000 --> q.latent maps went up to 50,000 (maybe in J rather than KJ?) - removed this to test something
print `r.mapcalc 'condens = if((swe.yesterday+snow)>0.0,q.latent/((2500.0+335.0)*1000.0)*100.0,0.0)' --o`; # in cm MZ 20170128; get daily average value by dividing 24 MZ 20190410
#print `r.mapcalc 'condens = if((swe.yesterday+snow)>0.0,if(q.latent/((2500.0+335.0)*1000.0)*100.0>0,q.latent/((2500.0+335.0)*1000.0)*100.0, max(q.latent/((2500.0+335.0)*1000.0)*100.0,-(swe.yesterday+snow))),0.0)' --o`; # in cm MZ 20170128; get daily average value by dividing 24 MZ 20190410

# q.cond is the heat from convective temp. gradients between air and snow
# density of air 1.29 kg/m^3, heat capacity of air 1 kJ/kg/C
# the wind roughness $rh is inputed in units of s/m
#  NOTE:  THIS NEEDS TO BE CONVERTED TO DAY/M OR HR/M ACORDING TO THE DESIRED TIMESTEP
# q.cond (KJ/m^2)

print `r.mapcalc 'q.sensible = 1.0*1.29*(tavg-tsnow_surf)/(rh)' --o`; # convert to KJ/(m^2 day) MZ 20200329
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
print "\n|----- END OF SNOW ACCUMULATION AND MELT MODEL -----|\n";
print "\n \n";

print `r.mapcalc 'water_input = snowmelt+throughfall' --o`;
print `r.mapcalc 'snowmelt = if(water_input>0.0,snowmelt*(water_input-road_runoff)/water_input,snowmelt)' --o`;
print `r.mapcalc 'throughfall = if(water_input>0.0,throughfall*(water_input-road_runoff)/water_input,throughfall)' --o`;

#  -----------------------------Initiate 6-hr Hydrology ----------------------------------
print "\n|----- Initial Conditions for 6-hr Temp Array -----|\n";
print "\n \n";
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
			print `r.mapcalc 'temp = $hrly_tmp + ($temp_LR/1000)*($el - $low_site_el)' --o`;
			print `r.mapcalc 'temp_sum = max(0.0,temp)+temp_sum' --o`;

		};
	print "\n \n";
	print "\n|----- 6 HOUR TEMPERATURE ARRAY -----|\n";
	print "\n \n";
	foreach $hrly_tmp (@hourly_tmp_array)
		{
			print `r.mapcalc 'temp = $hrly_tmp + ($temp_LR/1000)*(el - $low_site_el)' --o`;


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

	print `r.mapcalc 'ET_coeff = if(landuse==6.0,$cc_row_crop,if(landuse==5.0,$cc_grass,if(landuse==4.0,$cc_shrub,if(landuse==3.0,$cc_forest,if(landuse==2.0,$cc_urban,if(landuse==1.0,$cc_water,$cc_water))))))' --o`; 
	
	print `r.mapcalc 'root_storage_amt = if(landuse==6.0,storage_amt,if(landuse==5.0,storage_amt,if(landuse==4.0,storage_amt,if(landuse==3.0,storage_amt,if(landuse==2.0,storage_amt,if(landuse==1.0,storage_amt,storage_amt))))))' --o`;
	
	print `r.mapcalc 'actualET_flow = if(root_storage_amt>(root_zone*ETreduction_mc),min(max(0.0,(root_storage_amt-wiltpt_amt)),pet*ET_coeff*$temp_time_step/24.0),if(root_storage_amt>wiltpt_amt,min(max(0.0,(root_storage_amt-wiltpt_amt)),max(0.0,pet*((root_storage_amt/root_zone-wiltpt_amt/root_zone)/(ETreduction_mc-wiltpt_amt/root_zone))*ET_coeff*$temp_time_step/24.0)),0.0))' --o`; 
	
	print `r.mapcalc 'actualET_daily_flow = actualET_daily_flow+actualET_flow' --o`;
	
	print `r.mapcalc 'storage_amt = storage_amt-actualET_flow' --o`;
	
	print `r.mapcalc 'storage_amt_A = if(landuse==6.0,if(storage_amt<fieldcap_amt,fieldcap_amt_A*storage_amt/fieldcap_amt,fieldcap_amt_A),if(landuse==5.0,if(storage_amt<fieldcap_amt,fieldcap_amt_A*storage_amt/fieldcap_amt,fieldcap_amt_A),if(landuse==4.0,if(storage_amt<fieldcap_amt,fieldcap_amt_A*storage_amt/fieldcap_amt,fieldcap_amt_A),if(landuse==3.0,storage_amt_A-actualET_flow,if(landuse==2.0,if(storage_amt<fieldcap_amt,fieldcap_amt_A*storage_amt/fieldcap_amt,fieldcap_amt_A),if(landuse==1.0,if(storage_amt<fieldcap_amt,fieldcap_amt_A*storage_amt/fieldcap_amt,fieldcap_amt_A),0))))))' --o`;
	
	print `r.mapcalc 'storage_amt_B = storage_amt-storage_amt_A' --o`;


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

	print `r.mapcalc 'storage_amt = storage_amt-runoff_flow' --o`; # needs stream effect
	print `r.mapcalc 'runoff_daily_flow = runoff_daily_flow+runoff_flow' --o`;

	print `r.mapcalc 'saturation = storage_amt/sat_amt*100' --o`; # DJ 04/12/23-- made it so saturation is always 100$ at stream cells 

	}; # end of the 6 hours loop

	print "\n \n";
	print "\n|----- END OF 6 HOUR TIME LOOPING -----|\n";
	print "\n \n";

	print `r.mapcalc 'runoff_daily_flow = runoff_daily_flow' --o`; # dropped road runoff
	print `r.mapcalc 'input_daily_balance = input_daily - (snowmelt + throughfall)' --o`;
	#print `r.mapcalc 'psat_$year = psat_$year + if(saturation>=100,1,0)' --o`;


#____________________________________________________________________________________
#
#  7. Output maps
#____________________________________________________________________________________

print "\n|----- Febrauary runoff & saturation -----|\n";  
	# the idea here is too sum up the daily stats if it's february (month == 2) and then 
	# if the day is 60, which will be right at the end of february write out the file to 
	# an ascii (it will be hard to say if this works until we start trying to run the model)
	if ($month == 2) {
	  
	  #print `g.rename raster=runoff_feb_$year,runoff_feb --o`;
		`r.mapcalc 'runoff_feb_$year = runoff_feb_$year + runoff_daily_flow' --o`;
		
		#print `g.rename raster=A_amt_feb_$year,A_amt_feb_$year_new --o`;
		`r.mapcalc 'A_amt_feb_$year = A_amt_feb_$year + storage_amt_A' --o`;
		
		#print `g.rename raster=B_amt_feb_$year,B_amt_feb_$year_new --o`;
		`r.mapcalc 'B_amt_feb_$year = B_amt_feb_$year + storage_amt_B' --o`;
		
		}

	# this outputs accumulated and average february runoff but doesn't account for shortened 
	# february --> this may need to be changed.
	if ($doy == 62) {
		print `r.mapcalc 'avg_feb_runoff_$year = runoff_feb_$year / 29' --o`;
		print `r.out.gdal input=avg_feb_runoff_$year output=/$home_dir_name/$user/Dropbox/SMR_R/raw_data/smr_output/maps_$user/feb_avg_runoff_$year.tif --o`;
		
		print `r.mapcalc 'avg_A_amt_feb_$year = A_amt_feb_$year / 29' --o`;
		print `r.out.gdal input=avg_A_amt_feb_$year output=/$home_dir_name/$user/Dropbox/SMR_R/raw_data/smr_output/maps_$user/avg_A_amt_feb_$year.tif --o`;

		print `r.mapcalc 'avg_B_amt_feb_$year = B_amt_feb_$year / 29' --o`;
		print `r.out.gdal input=avg_B_amt_feb_$year output=/$home_dir_name/$user/Dropbox/SMR_R/raw_data/smr_output/maps_$user/avg_B_amt_feb_$year.tif --o`;
		}

print "\n|----- annual precip and pet -----|\n";  
  # maps of annual precipitation and pet to confirm distribution is correct
  # with respect to elevation
  if ($doy <= 365) {
    `r.mapcalc 'precip_$year = precip_$year + precip' --o`;
    `r.mapcalc 'pet_$year = pet_$year + pet' --o`;
  }
  
  if ($doy == 365) { # will need to change to 365
    print `r.out.gdal input=precip_$year output=/$home_dir_name/$user/Dropbox/SMR_R/raw_data/smr_output/maps_$user/precip_$year.tif --o`;
    print `r.out.gdal input=pet_$year output=/$home_dir_name/$user/Dropbox/SMR_R/raw_data/smr_output/maps_$user/pet_$year.tif --o`; 
  }
 
 print "\n|----- annual saturation -----|\n";   
  # count up number of days in a year each cell is at saturation
  if ($doy <= 365) { # need to make 365 eventually once markdown is set up
    `r.mapcalc 'psat_$year = psat_$year + if(saturation >= 100, 1, 0)' --o`; # make this 100 but want it to be showing up in ten days for now
  }
  
  # at the end of the year divide by number of days in the year and multiply by 100 to get the percentage of the year each cell spends at saturation
  if ($doy == 365) { # need to make 365 eventually once markdown is set up
    `r.mapcalc 'psat_$year = (psat_$year / 365) * 100' --o`; 
    print `r.out.gdal input=psat_$year output=/$home_dir_name/$user/Dropbox/SMR_R/raw_data/smr_output/maps_$user/psat_$year.tif --o`;
  }

	open($WATERSHED, '<', 'wshed_res_properties.ini') || "Can't open watershed properties file\n";
	while (<$WATERSHED>) {
		chop($_);
		($wshed_id,$area_cells,$res_vol,$res_coeff) = split;
			#print `r.mapcalc 'MASK = if(wsheds_all==$wshed_id,1,0)' --o`;
		
		# **********************************  1  ***************************************
		# Runoff output
			print `r.stats.zonal base=watershed cover=runoff_daily_flow out=temp1 method=sum --o --quiet`; # Calculates the sum of runoff_daily_flow for each cell in watershed 
      print "\n|-----The runoff_daily_flow is $temp1 -----|\n";			
			print $runoff_cm_{$wshed_id} = `r.stats -A -n -N input=temp1`; # A:  Print averaged values, n:  Do not report no data value, N:  Do not report cells where all maps have no data
			$runoff_cm_{$wshed_id} = $runoff_cm_{$wshed_id}*1; # sum of all cells in cm
			$runoff_cms_{$wshed_id} = ($runoff_cm_{$wshed_id}/100)*($gridsize*$gridsize)/($time_step*3600.0); # Calculate the runoff in m3/s (cubic meter per second) for sum of all cells
			$runoff_mm_{$wshed_id} = $runoff_cm_{$wshed_id}*10/$area_cells; # Convert runoff into area average depth (mm)


		# **********************************  3  ***************************************
		print "\n|----- annual average watershed precip -----|\n";	
		# Precip output
			print `r.stats.zonal base=watershed cover=precip out=temp3 method=sum --o --quiet`;
			print $precip_cm_{$wshed_id} = `r.stats -A -n -N input=temp3`;
			$precip_cm_{$wshed_id} = $precip_cm_{$wshed_id}*1;


		# **********************************  4  ***************************************
			print "\n|----- annual average watershed rain -----|\n";
		# Rain output
			print `r.stats.zonal base=watershed cover=rain out=temp4 method=sum --o --quiet`;
			print $rain_cm_{$wshed_id} = `r.stats -A input=temp4 null_value="null" `;
			$rain_cm_{$wshed_id} = $rain_cm_{$wshed_id}*1;


		# **********************************  5  ***************************************
			print "\n|----- annual average AET -----|\n";
		# Actual evaporation output
			print `r.stats.zonal base=watershed cover=actualET_daily_flow out=temp5 method=sum --o --quiet`;
			print $actualET_flow_cm_{$wshed_id} = `r.stats -A input=temp5 null_value="AET is null"`;
			$actualET_flow_cm_{$wshed_id} = $actualET_flow_cm_{$wshed_id}*1;


		# **********************************  6  ***************************************
		print "\n|----- annual average canopy ET -----|\n";
		# Canopy ET output
			print `r.stats.zonal base=watershed cover=canopy_evap out=temp6 method=sum --o --quiet`;
			print $canopy_evap_cm_{$wshed_id} = `r.stats -A input=temp6 nv= `;
			$canopy_evap_cm_{$wshed_id} = $canopy_evap_cm_{$wshed_id}*1;


		# **********************************  7  ***************************************
		print "\n|----- annual average snowmelt -----|\n";
		# Snowmelt output
			print `r.stats.zonal base=watershed cover=snowmelt out=temp7 method=sum --o --quiet`;
			print $snowmelt_cm_{$wshed_id} = `r.stats -A input=temp7 nv= `;
			$snowmelt_cm_{$wshed_id} = $snowmelt_cm_{$wshed_id}*1;


		# **********************************  8  ***************************************
		print "\n|----- annual average Storage -----|\n";
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
			$perc_cms_{$wshed_id} = $perc_cm_{$wshed_id}/100.0*($gridsize*$gridsize)/($time_step*3600.0);
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
			
		# **********************************  30 ***************************************
		# ice.content
			print `r.stats.zonal base=watershed cover=ice.content out=temp30 method=sum --o --quiet`;
			print $ice_content_{$wshed_id} = `r.stats -A input=temp30`;
			$ice_content_{$wshed_id} = $ice_content_{$wshed_id}*1;

		# ********************************** 31  ***************************************
		# liquid.water
			print `r.stats.zonal base=watershed cover=liquid.water out=temp31 method=sum --o --quiet`;
			print $liquid_water_{$wshed_id} = `r.stats -A input=temp31`;
			$liquid_water_{$wshed_id} = $liquid_water_{$wshed_id}*1;
			
		# **********************************  32 ***************************************
		# refreeze
			print `r.stats.zonal base=watershed cover=refreeze out=temp32 method=sum --o --quiet`;
			print $refreeze_{$wshed_id} = `r.stats -A input=temp32`;
			$refreeze_{$wshed_id} = $refreeze_{$wshed_id}*1;
			
		# **********************************  33 ***************************************
		# vap.d.air
			print `r.stats.zonal base=watershed cover=vap.d.air out=temp33 method=sum --o --quiet`;
			print $vap_d_air_{$wshed_id} = `r.stats -A input=temp33`;
			$vap_d_air_{$wshed_id} = $vap_d_air_{$wshed_id}*1;
	
		# **********************************  34 ***************************************
		# vap.d.snow
			print `r.stats.zonal base=watershed cover=vap.d.snow out=temp34 method=sum --o --quiet`;
			print $vap_d_snow_{$wshed_id} = `r.stats -A input=temp34`;
			$vap_d_snow_{$wshed_id} = $vap_d_snow_{$wshed_id}*1;
			
		# **********************************  35 ***************************************
		# u.surface
			print `r.stats.zonal base=watershed cover=u.surface out=temp35 method=sum --o --quiet`;
			print $u_surface_{$wshed_id} = `r.stats -A input=temp35`;
			$u_surface_{$wshed_id} = $u_surface_{$wshed_id}*1;
			
		# **********************************  36 ***************************************		
    # actualET_daily_flow
    	print `r.stats.zonal base=watershed cover=actualET_daily_flow out=temp36 method=sum --o --quiet`;
			print $actual_ET_daily_{$wshed_id} = `r.stats -A input=temp36`;
			$actual_ET_daily_{$wshed_id} = $actual_ET_daily_{$wshed_id}*1;

		#  Create mass balance output file
    open(OUT, ">>", "/$home_dir_name/$user/Dropbox/SMR_R/raw_data/smr_output/MFC_mass_balance_$user.csv") || die "Can't open weather file\n";
		print OUT "$wshed_id $date $year $runoff_cm_{$wshed_id} $precip_cm_{$wshed_id} $rain_cm_{$wshed_id} $actualET_flow_cm_{$wshed_id} $canopy_evap_cm_{$wshed_id} $snowmelt_cm_{$wshed_id} $storage_amt_cm_{$wshed_id} $throughfall_cm_{$wshed_id} $canopy_storage_amt_cm_{$wshed_id} $perc_cm_{$wshed_id} $Q_{$wshed_id} $swe_cm_{$wshed_id} $condens_cm_{$wshed_id} $snow_cm_{$wshed_id} $base_flow_{$wshed_id} $srad_{$wshed_id} $latent_{$wshed_id} $sensible_{$wshed_id} $lw_{$wshed_id} $q_rain_ground_cm_{$wshed_id} $q_total_{$wshed_id} $ice_content_{$wshed_id} $liquid_water_{$wshed_id} $refreeze_{$wshed_id} $vap_d_air_{$wshed_id} $vap_d_snow_{$wshed_id} $u_surface_{$wshed_id} $actual_ET_daily_{$wshed_id}\n";
		close(OUT); 

		print "\n \n";
		print "\n|----- END OF OUTPUTS ZONAL STATS -----|\n";
		print "\n \n";
		
	}
	close ($WATERSHED);


}
close ($WEATHER);

print "\n \n";
print "\n|---------- SMR-SAM HAS FINISHED RUNNING ----------|\n";
print "\n \n";

