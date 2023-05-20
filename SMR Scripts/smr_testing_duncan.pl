#!/usr/bin/perl

print `g.remove -f type=raster name=MASK`;
print `g.region rast=watershed`;
print `rm Q_* `;
print `rm M_* `;
print `rm R_* `;
print `rm tmp.txt`;
print `rm runoff_total_* `; # 03/19/2017 MZ
print `rm saturation_* `; # 03/19/2017 MZ
print `rm Psat_* `; # 04/09/2019 MZ

# Initialize variables
$gridsize = 30.0;		#m

# The time step variable refers to the time step at which snowmelt is calculated
$time_step = 24.0;	#hrs

# The temp time step determines how frequent the hydrology is simulated.  This requires
# temperature data in the weather file.  For example a 1 hr $temp_time_step requires hourly temperature data in the weather file
$temp_time_step = 6.0; #hrs

print `r.mapcalc 'MASK = if(wsheds_all>0,1,0)' --o`;

# --------------------------------------------------------------------------------
#    SLOPE & ASPECT FOR SAM
# --------------------------------------------------------------------------------

print 'r.slope.aspect elevation=el slope=slope aspect=aspect'

# --------------------------------------------------------------------------------
#    FLOW DIRECTION PROGRAM
# --------------------------------------------------------------------------------
#    A "flow unit" is an elevation difference of one unit for an adjacent
#    cell.  An elevation difference of one unit for a diagonal cell would
#    be 1/(square root of 2)=0.707  flow units.
# --------------------------------------------------------------------------------

print `r.mapcalc 'flowunits = max((el-if(isnull(el[-1,0]),el,el[-1,0])),0.0) + max((el-if(isnull(el[-1,1]),el,el[-1,1])),0.0)*0.707 + max((el-if(isnull(el[0,1]),el,el[0,1])),0.0) + max((el-if(isnull(el[1,1]),el,el[1,1])),0.0)*0.707 + max((el-if(isnull(el[1,0]),el,el[1,0])),0.0) + max((el-if(isnull(el[1,-1]),el,el[1,-1])),0.0)*0.707 + max((el-if(isnull(el[-1,-1]),el,el[-1,-1])),0.0)*0.707 + max((el-if(isnull(el[0,-1]),el,el[0,-1])),0.0)' --o`;

#    The following maps are the percent of the cell's neighbors in a direction
#    that will flow to the cell.  For example, "north", is the percent of the
#    lateral flow out of the cell to the north that will flow to the cell.

print `r.mapcalc 'north = if(isnull(flowunits[-1,0]),0,max((if(isnull(el[-1,0]),0,el[-1,0])-el),0.0)/flowunits[-1,0])' --o`;
print `r.mapcalc 'northeast = if(isnull(flowunits[-1,1]),0,max((if(isnull(el[-1,1]),0,el[-1,1])-el),0.0)/flowunits[-1,1]*0.707)' --o`;
print `r.mapcalc 'east = if(isnull(flowunits[0,1]),0,max((if(isnull(el[0,1]),0,el[0,1])-el),0.0)/flowunits[0,1])' --o`;
print `r.mapcalc 'southeast = if(isnull(flowunits[1,1]),0,max((if(isnull(el[1,1]),0,el[1,1])-el),0.0)/flowunits[1,1]*0.707)' --o`;
print `r.mapcalc 'south = if(isnull(flowunits[1,0]),0,max((if(isnull(el[1,0]),0,el[1,0])-el),0.0)/flowunits[1,0])' --o`;
print `r.mapcalc 'southwest = if(isnull(flowunits[1,-1]),0,max((if(isnull(el[1,-1]),0,el[1,-1])-el),0.0)/flowunits[1,-1]*0.707)' --o`;
print `r.mapcalc 'west = if(isnull(flowunits[0,-1]),0,max((if(isnull(el[0,-1]),0,el[0,-1])-el),0.0)/flowunits[0,-1])' --o`;
print `r.mapcalc 'northwest = if(isnull(flowunits[-1,-1]),0,max((if(isnull(el[-1,-1]),0,el[-1,-1])-el),0.0)/flowunits[-1,-1]*0.707)' --o`;

print `r.mapcalc 'north = if(isnull(north),0.0,north)' --o`;
print `r.mapcalc 'northeast = if(isnull(northeast),0.0,northeast)' --o`;
print `r.mapcalc 'east = if(isnull(east),0.0,east)' --o`;
print `r.mapcalc 'southeast = if(isnull(southeast),0.0,southeast)' --o`;
print `r.mapcalc 'south = if(isnull(south),0.0,south)' --o`;
print `r.mapcalc 'southwest = if(isnull(southwest),0.0,southwest)' --o`;
print `r.mapcalc 'west = if(isnull(west),0.0,west)' --o`;
print `r.mapcalc 'northwest = if(isnull(northwest),0.0,northwest)' --o`;

# --------------------------------------------------------------------------------
#     RESIDUAL MOISTURE CONTENT INITIALIZATION
# --------------------------------------------------------------------------------

print `r.mapcalc 'residual_mc_A = 0.02' --o`;
print `r.mapcalc 'residual_mc_B = 0.02' --o`;

# --------------------------------------------------------------------------------
# READ IN WATERSHED AND RESERVOIR PROPERTIES
# --------------------------------------------------------------------------------

# Open properties file
open($WATERSHED, '<', 'wshed_res_properties.ini') || die "Cannot open the watershed properties file: $!";

# Read each line of the file and assign values based on watershed ID
while ($line = <$WATERSHED>) {
	chop($line);
	($wshed_id,$area_cells,$res_vol,$res_coeff) = split(' ', $line);
	
	if ($wshed_id > 0.0) {
		$area_{$wshed_id} = $area_cells * $gridsize * $gridsize; #square meters
		$res_vol_{$wshed_id} = $res_vol;
		$res_coeff_{$wshed_id} = $res_coeff;
		$base_flow_{$wshed_id} = 0.0;
		
		print "|---- Properties of watershed with ID: $wshed_id ----|\n";
		print "number of cells in watershed: $area_cells\n";
		print "cell size (m): $gridsize\n";
		print "total watershed area (m^2): $area_{$wshed_id}\n";
        print "res_vol: $res_vol_{$wshed_id}\n";
        print "res_coeff: $res_coeff_{$wshed_id}\n";
        print "base_flow: $base_flow_{$wshed_id}\n\n";
		print "|----------------------------------------------------|\n"
	};
}

# Close the watershed properties file
close($WATERSHED) || die "Cannot close the watershed properties file: $!";


# initialize landuse and roads maps and associated parameters
	print `r.mapcalc 'roads = 0.0' --o`;

#  Landuse map legend:  1 = 100% forest cover, 2 = partial cut, 3 = clear cut
	print `r.mapcalc 'landuse = 1.0' --o`;

#  assume max_canopy storage is 0.1 cm for partial cut and 0.2 cm for 100% forest cover
	#print `r.mapcalc 'max_canopy_storage_amt = if(landuse==3.0,0.0,if(landuse==2.0,0.1,0.2))' --o`;
    print `r.mapcalc 'max_canopy_storage_amt = if(landuse==6.0, 0.0,
												if(landuse==5.0, 0.0,
												if(landuse==4.0, 0.0,
												if(landuse==3.0, 0.0,
												if(landuse==2.0, 0.1,
												if(landuse==1.0, 0.2,
												 0.0))))))' --o`;

	#print `r.mapcalc 'canopy_cover = if(landuse==3.0,0.01,if(landuse==2.0,0.5,0.99))' --o`; # MZ 20200510 add canopy cover fraction for solar radiation calc
	print `r.mapcalc 'canopy_cover = if(landuse==6.0, 0.0,
												if(landuse==5.0, 0.0,
												if(landuse==4.0, 0.0,
												if(landuse==3.0, 0.01,
												if(landuse==2.0, 0.5,
												if(landuse==1.0, 0.99,
												0.0))))))' 
												--o`;

#	print `r.mapcalc 'kfactor = if(landuse==3.0,0.734,if(landuse==2.0,0.734,0.719))*$time_step/24.0' --o`;
#	print `r.mapcalc 'tbase = if(landuse==3.0,1.54,if(landuse==2.0,1.54,2.14))' --o`;

	#print `r.mapcalc 'kfactor = if(landuse==3.0,0.71,if(landuse==2.0,0.695,0.68))*$time_step/24.0' --o`;  #from output 68
	print `r.mapcalc 'kfactor = if(landuse==6.0, 0.0,
								if(landuse==5.0, 0.0,
								if(landuse==4.0, 0.0,
								if(landuse==3.0, 0.71,
								if(landuse==2.0, 0.695,
								if(landuse==1.0, 0.68,
								0.0))))))
								*$time_step/24.0' 
								--o`;  #from output 68


	#print `r.mapcalc 'tbase = if(landuse==3.0,1.7,if(landuse==2.0,1.9,2.1))' --o`;
	print `r.mapcalc 'tbase = if(landuse==6.0, 0.0,
							  if(landuse==5.0, 0.0,
							  if(landuse==4.0, 0.0,
							  if(landuse==3.0, 1.7,
							  if(landuse==2.0, 1.9,
							  if(landuse==1.0, 2.1, 
							  0.0))))))' 
							  --o`;


#	print `r.mapcalc 'kfactor = if(landuse==3.0,0.61,if(landuse==2.0,0.595,0.58))*$time_step/24.0' --o`;   # Change 1991-2001 to the values that fit cutting years MZ 20171109
#	print `r.mapcalc 'tbase = if(landuse==3.0,1.7,if(landuse==2.0,1.9,2.1))' --o`;


	#print `r.mapcalc 'root_zone = if(landuse==3.0,soil_depth_A,soil_depth)' --o`;
	print `r.mapcalc 'root_zone = if(landuse==6.0, 0.0,
								  if(landuse==5.0, 0.0,
								  if(landuse==4.0, 0.0,
								  if(landuse==3.0, soil_depth_A,
								  if(landuse==2.0, soil_depth,
								  if(landuse==1.0, soil_depth, 
								  0.0))))))' 
								  --o`;


	#print `r.mapcalc 'wiltpt_amt = if(landuse==3.0,wiltpt_mc_A*soil_depth_A,wiltpt_mc_A*soil_depth_A+wiltpt_mc_B*soil_depth_B)' --o`;
	print `r.mapcalc 'wiltpt_amt = if(landuse==6.0,wiltpt_mc_A*soil_depth_A+wiltpt_mc_B*soil_depth_B,
								   if(landuse==5.0,wiltpt_mc_A*soil_depth_A+wiltpt_mc_B*soil_depth_B,
								   if(landuse==4.0,wiltpt_mc_A*soil_depth_A+wiltpt_mc_B*soil_depth_B,
								   if(landuse==3.0,wiltpt_mc_A*soil_depth_A,
								   if(landuse==2.0,wiltpt_mc_A*soil_depth_A+wiltpt_mc_B*soil_depth_B,
								   if(landuse==1.0,wiltpt_mc_A*soil_depth_A+wiltpt_mc_B*soil_depth_B,
								   0.0))))))' 
								   --o`;


	#print `r.mapcalc 'ETreduction_mc = if(landuse==3.0,fieldcap_amt_A*0.8/soil_depth_A,fieldcap_amt*0.8/soil_depth)' --o`;
	print `r.mapcalc 'ETreduction_mc = if(landuse==6.0,fieldcap_amt*0.8/soil_depth,
								       if(landuse==5.0,fieldcap_amt*0.8/soil_depth,
								   	   if(landuse==4.0,fieldcap_amt*0.8/soil_depth,
								   	   if(landuse==3.0,fieldcap_amt_A*0.8/soil_depth_A,
								   	   if(landuse==2.0,fieldcap_amt*0.8/soil_depth,
								       if(landuse==1.0,fieldcap_amt*0.8/soil_depth,
								       0.0))))))' 
								       --o`;


#	$tmax_rain = 1.3;
#	$tmin_snow = -0.6;
#	$tmax_rain = 2.0; # higher temperature to move streamflow late. MZ 20190607
#	$tmin_snow = 0.0;

	#print `r.mapcalc 'tmax_rain = if(landuse==3.0,3.1,if(landuse==2.0,1.8,-0.5))' --o`;
	print `r.mapcalc 'tmax_rain = if(landuse==6.0,0,
								  if(landuse==5.0,0,
								  if(landuse==4.0,0,
								  if(landuse==3.0,3.1,
								  if(landuse==2.0,1.8,
								  if(landuse==1.0,-0.5,
								  0))))))' 
								  --o`;


	#print `r.mapcalc 'tmin_snow = if(landuse==3.0,0.9,if(landuse==2.0,-1.0,-3.0))' --o`;
	print `r.mapcalc 'tmin_snow = if(landuse==6.0,0,
							  	  if(landuse==5.0,0,
							      if(landuse==4.0,0,
							      if(landuse==3.0,0.9,
							      if(landuse==2.0,-1.0,
							      if(landuse==1.0,-3.0,
							      0.0))))))' 
							      --o`;


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

print `r.mapcalc 'runoff_1990 = 0.0' --o`;
print `r.mapcalc 'runoff_1991 = 0.0' --o`;
print `r.mapcalc 'runoff_1992 = 0.0' --o`;
print `r.mapcalc 'runoff_1993 = 0.0' --o`;
print `r.mapcalc 'runoff_1994 = 0.0' --o`;
print `r.mapcalc 'runoff_1995 = 0.0' --o`;
print `r.mapcalc 'runoff_1996 = 0.0' --o`;
print `r.mapcalc 'runoff_1997 = 0.0' --o`;
print `r.mapcalc 'runoff_1998 = 0.0' --o`;
print `r.mapcalc 'runoff_1999 = 0.0' --o`;
print `r.mapcalc 'runoff_2000 = 0.0' --o`; 
print `r.mapcalc 'runoff_2001 = 0.0' --o`;
print `r.mapcalc 'runoff_2002 = 0.0' --o`;
print `r.mapcalc 'runoff_2003 = 0.0' --o`;
print `r.mapcalc 'runoff_2004 = 0.0' --o`;
print `r.mapcalc 'runoff_2005 = 0.0' --o`;
print `r.mapcalc 'runoff_2006 = 0.0' --o`;
print `r.mapcalc 'runoff_2007 = 0.0' --o`;
print `r.mapcalc 'runoff_2008 = 0.0' --o`;
print `r.mapcalc 'runoff_2009 = 0.0' --o`;
print `r.mapcalc 'runoff_2010 = 0.0' --o`;
print `r.mapcalc 'runoff_2011 = 0.0' --o`;
print `r.mapcalc 'runoff_2012 = 0.0' --o`;
print `r.mapcalc 'runoff_2013 = 0.0' --o`;

print `r.mapcalc 'Psat_1990 = 0.0' --o`;
print `r.mapcalc 'Psat_1991 = 0.0' --o`;
print `r.mapcalc 'Psat_1992 = 0.0' --o`;
print `r.mapcalc 'Psat_1993 = 0.0' --o`;
print `r.mapcalc 'Psat_1994 = 0.0' --o`;
print `r.mapcalc 'Psat_1995 = 0.0' --o`;
print `r.mapcalc 'Psat_1996 = 0.0' --o`;
print `r.mapcalc 'Psat_1997 = 0.0' --o`;
print `r.mapcalc 'Psat_1998 = 0.0' --o`;
print `r.mapcalc 'Psat_1999 = 0.0' --o`;
print `r.mapcalc 'Psat_2000 = 0.0' --o`;
print `r.mapcalc 'Psat_2001 = 0.0' --o`;
print `r.mapcalc 'Psat_2002 = 0.0' --o`;
print `r.mapcalc 'Psat_2003 = 0.0' --o`;
print `r.mapcalc 'Psat_2004 = 0.0' --o`;
print `r.mapcalc 'Psat_2005 = 0.0' --o`;
print `r.mapcalc 'Psat_2006 = 0.0' --o`;
print `r.mapcalc 'Psat_2007 = 0.0' --o`;
print `r.mapcalc 'Psat_2008 = 0.0' --o`;
print `r.mapcalc 'Psat_2009 = 0.0' --o`;
print `r.mapcalc 'Psat_2010 = 0.0' --o`;
print `r.mapcalc 'Psat_2011 = 0.0' --o`;
print `r.mapcalc 'Psat_2012 = 0.0' --o`;
print `r.mapcalc 'Psat_2013 = 0.0' --o`;


#____________________________________________________________________________________
#
# START READING WEATHER DATA
#____________________________________________________________________________________

#  This set of commands splits a tab delimited array using a while loop
open (WEATHER, "<wea_mica_24hr_1990_2013_output_wy_ccopen_SAM_cloud_rh_2006.txt") || die "Can't open file\n"; # 20170912 MZ  ###from 20171020, use 'wea_mica_24hr_1990_2013_output_wy.txt', all the $year means wateryear

while (<WEATHER>) {
	chop($_);
	($date,$year,$tmax,$tmin,$tavg,$precip,$pet_snotel,$cc_forest,$cc_partial,$cc_open,$output,$t0,$t1,$t2,$t3,$tdew,$cloud,$rh_snow,$rh_veg,$l_turb) = split(/\t/,$_);
	print " $date | $year | Tair = $tavg °C | Precip = $precip cm | PET = $pet_snotel cm | Output = $output\n";
	@hourly_tmp_array = ($t0,$t1,$t2,$t3);


#  changed snowmelt so that clear cut melt is the same as partial cut melt
#	print `r.mapcalc 'kfactor = if(landuse==3.0,0.734,if(landuse==2.0,0.734,0.719))*$time_step/24.0' --o`;
#	print `r.mapcalc 'tbase = if(landuse==3.0,1.54,if(landuse==2.0,1.54,2.14))' --o`;
	#print `r.mapcalc 'kfactor = if(landuse==3.0,0.71,if(landuse==2.0,0.695,0.68))*$time_step/24.0' --o`; # Degree-day
	print `r.mapcalc 'kfactor = if(landuse==6.0,0.0,
								if(landuse==5.0,0.0,
								if(landuse==4.0,0.0,
								if(landuse==3.0,0.71,
								if(landuse==2.0,0.695,
								if(landuse==1.0,0.68,
								0.0))))))' 
								--o`

	#print `r.mapcalc 'tbase = if(landuse==3.0,1.7,if(landuse==2.0,1.9,2.1))' --o`; # Degree-day
	print `r.mapcalc 'tbase = if(landuse==6.0,0.0,
							  if(landuse==5.0,0.0,
							  if(landuse==4.0,0.0,
							  if(landuse==3.0,1.7,
							  if(landuse==2.0,1.9,
							  if(landuse==1.0,2.1,
							  0.0))))))' 
							  --o`;

#  Increased the max_canopy_storage_amt in complete forest cover to 0.3
	#print `r.mapcalc 'max_canopy_storage_amt = if(landuse==3.0,0.0,if(landuse==2.0,0.1,0.3))' --o`;
	print `r.mapcalc 'max_canopy_storage_amt = if(landuse==6.0,0.0,
							  if(landuse==5.0,0.0,
							  if(landuse==4.0,0.0,
							  if(landuse==3.0,0.0,
							  if(landuse==2.0,2.0,
							  if(landuse==1.0,0.3,
							  0.0))))))' 
							  --o`;

	#print `r.mapcalc 'canopy_cover = if(landuse==3.0,0.01,if(landuse==2.0,0.5,0.99))' --o`; # MZ 20200510 add canopy cover fraction for solar radiation calc
	print `r.mapcalc 'canopy_cover = if(landuse==6.0,0.0,
									 if(landuse==5.0,0.0,
									 if(landuse==4.0,0.0,
									 if(landuse==3.0,0.01,
									 if(landuse==2.0,0.5,
									 if(landuse==1.0,0.99,
									 0.0))))))'
									 --o`


	print `r.mapcalc 'root_zone = if(landuse==3.0,soil_depth_A,soil_depth)' --o`;
	print `r.mapcalc 'root_zone = if(landuse==6.0,soil_depth,
								  if(landuse==5.0,soil_depth,
								  if(landuse==4.0,soil_depth,
								  if(landuse==3.0,soil_depth_A,
								  if(landuse==2.0,soil_depth,
								  if(landuse==1.0,soil_depth,
								  0.0))))))' 
								  --o`

	#print `r.mapcalc 'wiltpt_amt = if(landuse==3.0,wiltpt_mc_A*soil_depth_A,wiltpt_mc_A*soil_depth_A+wiltpt_mc_B*soil_depth_B)' --o`;
	print `r.mapcalc 'wiltpt_amt = if(landuse==6.0, wiltpt_mc_A*soil_depth_A + wiltpt_mc_B*soil_depth_B,
								   if(landuse==5.0, wiltpt_mc_A*soil_depth_A + wiltpt_mc_B*soil_depth_B,
      							   if(landuse==4.0, wiltpt_mc_A*soil_depth_A + wiltpt_mc_B*soil_depth_B,
         						   if(landuse==3.0, wiltpt_mc_A*soil_depth_A,
            					   if(landuse==2.0, wiltpt_mc_A*soil_depth_A + wiltpt_mc_B*soil_depth_B,
               					   if(landuse==1.0, wiltpt_mc_A*soil_depth_A + wiltpt_mc_B*soil_depth_B,
                 				   0.0)))))))' 
								   --o`;

	#print `r.mapcalc 'ETreduction_mc = if(landuse==3.0,fieldcap_amt_A*0.8/soil_depth_A,fieldcap_amt*0.8/soil_depth)' --o`;

	print `r.mapcalc 'ETreduction_mc = if(landuse==6.0, fieldcap_amt*0.8/soil_depth,
   									   if(landuse==5.0, fieldcap_amt*0.8/soil_depth,
      								   if(landuse==4.0, fieldcap_amt*0.8/soil_depth,
        							   if(landuse==3.0, fieldcap_amt_A*0.8/soil_depth_A,
          							   if(landuse==2.0, fieldcap_amt*0.8/soil_depth,
               						   if(landuse==1.0, fieldcap_amt*0.8/soil_depth, 
									   0.0)))))))'
									    --o`;

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
#  an equation between PET and elevation.  PET_st_maries = 1.1088 * PET_snotel + 0.0097 (cm)

# print `r.mapcalc 'precip = $precip' --o`; # MZ 20170210
# print `r.mapcalc 'tavg = $tavg' --o`;
# print `r.mapcalc 'pet = $pet_snotel' --o`;

print `r.mapcalc 'precip = (el*0.08136+27.9)/(145.7)*($precip)' --o`;
print `r.mapcalc 'tavg = $tavg-0.00576*(el-1447.8)' --o`;
print `r.mapcalc 'tdew = $tdew-0.00576*(el-1447.8)' --o`; # temperature distribution MZ 20190205
print `r.mapcalc 'rain = if(tavg>tmax_rain,precip,if(tavg<tmin_snow,0.0,(tavg-tmin_snow)/(tmax_rain-tmin_snow)*precip))' --o`;
print `r.mapcalc 'snow = precip-rain' --o`;
print `r.mapcalc 'pet_st_maries = 1.1088*$pet_snotel+0.0097' --o`;
print `r.mapcalc 'pet = (pet_st_maries-$pet_snotel)*(el-675.4)/(1447.8-675.4)+pet_st_maries' --o`;

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

#  ------------------------------- SAM --------------------------------------
#  Snow accumulation and melt (SAM) model
#  Modified by MZ on 20190126
#  albedo approximated from DHSVM relationship
#  assumes the albedo of the soil is 0.2
#  albedo (100 x %)

#print `r.mapcalc 'rh = if(swe>0,if(landuse==1 || landuse==2,$rh_veg,$rh_snow),$rh_veg)' --o`; # modified MZ 20190205
#print `r.mapcalc 'rh = if(swe>0,$rh_snow,$rh_veg)' --o`; # modified MZ 20190205
print `r.mapcalc 'rh = if(swe>0,$rh_snow/(1.0-(canopy_cover/3.0)),$rh_veg/(1.0-(canopy_cover/3.0)))' --o`; # add canopy cover effects on rh MZ 20200601
print `r.mapcalc 'snow.age = if(snow>0.0 && throughfall==0.0,1.0,snow.age+1.0)' --o`; # modified MZ 20190210
print `r.mapcalc 'albedo = if(swe.yesterday+snow>0.0,min(0.95,0.7383*snow.age^(-0.1908)),0.2)' --o`;

#  Calculate the diffuse transmissivity following Bristow and Campbell (1985)
#  The clear sky transmissivity for Troy, ID is 0.75
print `r.mapcalc 't_diff = if($cloud==1.0,0.1,0.75*(1-$cloud)*(1-exp(-0.6*$cloud/((0.75-0.4)*(1-$cloud)))))' --o`;

#  Calculate the real-sky diffuse radiation coefficient
#  The clear sky diffuse transmissivity for the Palouse is 0.10 (Bristow and Campbell, 1985)
print `r.mapcalc 'coefdh = t_diff/0.10' --o`;
# print `r.mapcalc 'coefdh = t_diff' --o`; # 20190404 MZ Test

#  Calculate the real-sky beam radiation coefficient
#  The clear sky transmissivity for Troy, ID is 0.75
print `r.mapcalc 'coefbh = max(0.0,(1-$cloud)-t_diff/0.75)' --o`;
# print `r.mapcalc 'coefbh = max(0.0,0.75*(1-$cloud)-t_diff)' --o`; # MZ 20190409 t_total = t_diff + t_beam, t_total = B*(1-cloud), B is maximum clear-sky transmissivity (0.75)


#  Calculate the beam, diffuse, and reflected irradiance at time $s_time
#  Change to daily time step so the output of beam, diffuse and reflected radiation is in Wh/m2/day, MZ 20190127
print `r.sun elevation=el aspect=aspect slope=slope linke_value=$l_turb albedo=albedo coeff_bh=coefbh coeff_dh=coefdh beam_rad=beam_rad diff_rad=diff_rad refl_rad=refl_rad day=$date nprocs=12 --o`;

#  Rad_tot is the total radiation in units of KJ/m^2
#  Convert to daily value by multiply 24, MZ 20170127; Update: Do not multiply 24 because the output of beam etc is already daily in Wh/m2/day. 1Wh/m2/day = 3.6KJ/m2/day
#print `r.mapcalc 'q.srad = (1.0-albedo)*(beam_rad+diff_rad+refl_rad)*3600.0/1000.0' --o`;
print `r.mapcalc 'q.srad = (1.0-canopy_cover)*(1.0-albedo)*(beam_rad+diff_rad+refl_rad)*3600.0/1000.0' --o`; # add canopy cover impact MZ 20200510

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

print `r.mapcalc 'q.latent = (2500.0+335.0)*(vap.d.air-vap.d.snow)/(rh/(3600.0*24.0))' --o`;
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

#print `r.mapcalc 'snowmelt = (swe.yesterday-swe)+rain+snow+condens' --o`;
print `r.mapcalc 'snowmelt = (swe.yesterday-swe)+snow+condens' --o`;

print `r.mapcalc 'swe.y.save = swe.yesterday' --o`; # testing MA 20200329
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
#if ($date == 279){die};
#  -----------------------------End of SAM ----------------------------------

# assume road width is 5 m and assume road surface storage is 1.5 cm 
# assume once surface storage is exceeded all runoff reaches stream
# print `echo "road_runoff=if(roads==1.0,max(0.0,snowmelt+throughfall-1.5)*5.0/$gridsize,0.0)" | r.mapcalc`;
print `r.mapcalc 'road_runoff = if(roads==1.0,max(0.0,snowmelt+throughfall-1.5)*5.0/$gridsize,0.0)' --o`;
print `r.mapcalc 'water_input = snowmelt+throughfall' --o`;
print `r.mapcalc 'snowmelt = if(water_input>0.0,snowmelt*(water_input-road_runoff)/water_input,snowmelt)' --o`;
print `r.mapcalc 'throughfall = if(water_input>0.0,throughfall*(water_input-road_runoff)/water_input,throughfall)' --o`;

print `r.mapcalc 'temp_sum = 0.0' --o`;
print `r.mapcalc 'actualET_daily_flow = 0.0' --o`;
print `r.mapcalc 'perc_daily_flow = 0.0' --o`;
print `r.mapcalc 'runoff_daily_flow = 0.0' --o`;
print `r.mapcalc 'lateral_daily_out = 0.0' --o`;
print `r.mapcalc 'lateral_daily_in = 0.0' --o`;
print `r.mapcalc 'input_daily = 0.0' --o`; # daily input MZ 2017.2.16
print `r.mapcalc 'SM_test = 0.0' --o`; # test MZ 2017.2.21

print `r.mapcalc 'mass_daily_balance = 0.0' --o`; # daily mass balance MZ 2017.2.5
print `r.mapcalc 'input_daily_balance = 0.0' --o`; # daily input balance MZ 2017.2.16


	foreach $hrly_tmp (@hourly_tmp_array)
		{
#			print `r.mapcalc 'temp = $hrly_tmp' --o`; # MZ 20170210
			print `r.mapcalc 'temp = $hrly_tmp-0.00531*(el-1447.8)' --o`;
			print `r.mapcalc 'temp_sum = max(0.0,temp)+temp_sum' --o`;

		};

	foreach $hrly_tmp (@hourly_tmp_array)
		{
#			print `r.mapcalc 'temp = $hrly_tmp' --o`; # MZ 20170210
			print `r.mapcalc 'temp = $hrly_tmp-0.00531*(el-1447.8)' --o`;

			print " $hrly_tmp \n";
#  -----------------------------------------------------------------------
#  -----------------------------------------------------------------------

#____________________________________________________________________________________
#
#  2.  DISTRIBUTE WATER TO EACH SOIL LAYER
#____________________________________________________________________________________

#    Storage is distributed assuming that all water infiltrates.  When total field capacity of
#    the entire soil profile is not exceeded, moisture is evenly distributed.
#    When total field capacity is exceeded, saturation begins

print `r.mapcalc 'input = snowmelt*max(temp,0.0)/if(temp_sum==0.0,1.0,temp_sum)+throughfall*$temp_time_step/24.0' --o`;
print `r.mapcalc 'input_daily = input_daily + input' --o`;
# print `r.mapcalc 'SM_test = SM_test + max(temp,0.0)/if(temp_sum==0.0,1.0,temp_sum)' --o`;


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

# water distribution MZ 20171205
print `r.mapcalc 'storage_amt_A_tmp = if(storage_amt<fieldcap_amt,min(storage_amt_A,fieldcap_amt_A),storage_amt_A)' --o`;

print `r.mapcalc 'storage_amt_B = if(storage_amt<fieldcap_amt,storage_amt-storage_amt_A_tmp,if((storage_amt-fieldcap_amt_A)<sat_mc_B*soil_depth_B,storage_amt-fieldcap_amt_A,sat_mc_B*soil_depth_B))' --o`;

print `r.mapcalc 'storage_amt_A = if(storage_amt<fieldcap_amt,storage_amt_A_tmp,if((storage_amt-storage_amt_B)<sat_mc_A*soil_depth_A,storage_amt-storage_amt_B,sat_mc_A*soil_depth_A))' --o`;

print `r.mapcalc 'storage_diff = storage_amt - storage_amt_A - storage_amt_B' --o`;

#******************************** END FIX ***************************************


#____________________________________________________________________________________
#
#    4. EVAPOTRANSPIRATION CALCULATION
#____________________________________________________________________________________

#    Actual ET is based on moisture content, ETcoeff, and potential ET.
#    ET takes place at potential rate at ETreduction_mc or higher
#    moisture contents, stops at wilting point or below, and is linearly
#    related to moisture content between (Thornthwaite-Mather assumption)

# print `r.mapcalc 'ET_coeff = if(landuse==3.0,$cc_open+20.0,if(landuse==2.0,100.0,130.0))' --o`;
#print `r.mapcalc 'ET_coeff = if(landuse==3.0,$cc_open-10,if(landuse==2.0,80.0,120.0))' --o`;  # test!!!decrease crop coefficient to decrease ET, so the runoff increases MZ 1/8/2017

#print `r.mapcalc 'ET_coeff = if(landuse==3.0,$cc_open-10,if(landuse==2.0,$cc_partial-10,$cc_forest-10))' --o`;  # test!!!decrease crop coefficient to decrease ET, so the runoff increases MZ 12/14/2017
print `r.mapcalc 'r.mapcalc = if(landuse==6.0,$cc_open-10,
						if(landuse==5.0,$cc_open-10,
						if(landuse==4.0,$cc_open-10,
						if(landuse==3.0,$cc_open-10,
						if(landuse==2.0,$cc_partial-10,
						if(landuse==1.0,$cc_forest-10,
						$cc_open-10))))))' 
						--o`;

#print `r.mapcalc 'root_storage_amt = if(landuse==3.0,storage_amt_A,storage_amt)' --o`;
print `r.mapcalc 'root_storage_amt = if(landuse==6.0,storage_amt,
									 if(landuse==5.0,storage_amt,
									 if(landuse==4.0,storage_amt,
									 if(landuse==3.0,storage_amt_A,
									 if(landuse==2.0,storage_amt,
									 if(landuse==1.0,storage_amt,
									 storage_amt))))))'
								     --o`;


#print `r.mapcalc 'actualET_flow = if(root_storage_amt>(root_zone*ETreduction_mc),min((root_storage_amt-wiltpt_amt),pet*ET_coeff/100.0*$temp_time_step/24.0),if(root_storage_amt>wiltpt_amt,min((root_storage_amt-wiltpt_amt),pet*((root_storage_amt/root_zone-wiltpt_amt/root_zone)/(ETreduction_mc-wiltpt_amt/root_zone))*ET_coeff/100.0*$temp_time_step/24.0),0.0))' --o`;
print `r.mapcalc 'actualET_flow = if(root_storage_amt>(root_zone*ETreduction_mc),min(max(0.0,(root_storage_amt-wiltpt_amt)),pet*ET_coeff/100.0*$temp_time_step/24.0),if(root_storage_amt>wiltpt_amt,min(max(0.0,(root_storage_amt-wiltpt_amt)),max(0.0,pet*((root_storage_amt/root_zone-wiltpt_amt/root_zone)/(ETreduction_mc-wiltpt_amt/root_zone))*ET_coeff/100.0*$temp_time_step/24.0)),0.0))' --o`; ### FIX  to avoid negative values - MZ 20171030 Fixed

print `r.mapcalc 'actualET_daily_flow = actualET_daily_flow+actualET_flow' --o`;
print `r.mapcalc 'storage_amt = storage_amt-actualET_flow' --o`;

#print `r.mapcalc 'storage_amt_A = if(landuse==3.0,storage_amt_A-actualET_flow,if(storage_amt<fieldcap_amt,fieldcap_amt_A*storage_amt/fieldcap_amt,fieldcap_amt_A))' --o`;
print `r.mapcalc 'storage_amt_A = if(landuse==6.0,if(storage_amt<fieldcap_amt,fieldcap_amt_A*storage_amt/fieldcap_amt,fieldcap_amt_A),
								  if(landuse==5.0,if(storage_amt<fieldcap_amt,fieldcap_amt_A*storage_amt/fieldcap_amt,fieldcap_amt_A),
								  if(landuse==4.0,if(storage_amt<fieldcap_amt,fieldcap_amt_A*storage_amt/fieldcap_amt,fieldcap_amt_A),
								  if(landuse==3.0,storage_amt_A-actualET_flow,
								  if(landuse==2.0,if(storage_amt<fieldcap_amt,fieldcap_amt_A*storage_amt/fieldcap_amt,fieldcap_amt_A),
								  if(landuse==1.0,if(storage_amt<fieldcap_amt,fieldcap_amt_A*storage_amt/fieldcap_amt,fieldcap_amt_A),
								  0))))))'
								  --o`;

print `r.mapcalc 'storage_amt_B = storage_amt-storage_amt_A' --o`; # 20171207 MZ

#____________________________________________________________________________________
#
#    5. PERCOLATION
#____________________________________________________________________________________

#    Moisture above fieldcap.mc can percolate to the restricting layer
#    if the hydraulic conductivity of the restricting layer is adequate.
#    The amount of percolation is substracted from the storage.

#print `r.mapcalc 'perc = min(max(storage_amt-fieldcap_amt,0.0),Ksubsurface*$temp_time_step/24.0)' --o`; # commented on 20171205 MZ
print `r.mapcalc 'perc = min(min(max(storage_amt-fieldcap_amt,0.0),storage_amt_B),Ksubsurface*$temp_time_step/24.0)' --o`;
print `r.mapcalc 'storage_amt = storage_amt - perc' --o`; # MZ 20171205
print `r.mapcalc 'storage_amt_B = storage_amt_B - perc' --o`; # MZ 20171205

# water distribution MZ 20171207
print `r.mapcalc 'storage_amt_A_tmp = if(storage_amt<fieldcap_amt,min(storage_amt_A,fieldcap_amt_A),storage_amt_A)' --o`;

print `r.mapcalc 'storage_amt_B = if(storage_amt<fieldcap_amt,storage_amt-storage_amt_A_tmp,if((storage_amt-fieldcap_amt_A)<sat_mc_B*soil_depth_B,storage_amt-fieldcap_amt_A,sat_mc_B*soil_depth_B))' --o`;

print `r.mapcalc 'storage_amt_A = if(storage_amt<fieldcap_amt,storage_amt_A_tmp,if((storage_amt-storage_amt_B)<sat_mc_A*soil_depth_A,storage_amt-storage_amt_B,sat_mc_A*soil_depth_A))' --o`;


print `r.mapcalc 'perc_daily_flow = perc_daily_flow + perc' --o`;
#print `r.mapcalc 'storage_amt = storage_amt-perc' --o`; # commented on 20171205 MZ


#____________________________________________________________________________________
#
#    6. RUNOFF
#____________________________________________________________________________________

#    Moisture above saturation in each cell becomes runoff.

#print `r.mapcalc 'runoff_flow = if(storage_amt>sat_amt,min(storage_amt-sat_amt,storage_amt_A),0.0)' --o`; # MZ 20171205; comment 20180103 MZ
print `r.mapcalc 'runoff_flow = if(storage_amt>sat_amt,storage_amt-sat_amt,0.0)' --o`; # commented on 20171205 MZ; Uncomment 20180103 MZ

#print `r.mapcalc 'storage_amt_A = storage_amt_A-runoff_flow' --o`; # MZ 20171205 comment on 20171207 storage_amt_A should never have runoff coming out because it is always below sat_A
#print `r.mapcalc 'storage_amt_B = storage_amt - storage_amt_A' --o`; # MZ 20171205, comment on 20171207

print `r.mapcalc 'storage_amt = storage_amt-runoff_flow' --o`; # commented on 20171205 MZ, uncomment on 20171207
print `r.mapcalc 'runoff_daily_flow = runoff_daily_flow+runoff_flow' --o`;

print `r.mapcalc 'saturation = storage_amt/sat_amt*100.0' --o`;
print `r.mapcalc 'runoff_$year = runoff_flow+runoff_$year' --o`;

print `r.reclass input=saturation output=sat_index rules=sat_rcl.txt --o`;

}; # end of the 6 hours loop

print `r.mapcalc 'runoff_daily_flow = runoff_daily_flow+road_runoff' --o`;
# print `r.mapcalc 'input_daily_balance = input_daily -(snowmelt + throughfall)' --o`; # MZ 20170216

print `r.mapcalc 'Psat_$year = Psat_$year + if(saturation>=100,1,0)' --o`; # changed to if saturation >=100 just in case 20180103 MZ


#____________________________________________________________________________________
#
#  7. Output maps
#____________________________________________________________________________________

# will add to this as it makes sense

};
close (WEATHER);

