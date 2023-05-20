#!/usr/bin/perl

#    program smr_setup_mica.pl
#    One soil layer program for Mica creek watershed
#    This program creates the maps needed in the Soil Moisture Routing model

#    INPUT REQUIREMENTS
#    1. Input maps and soil parameters
#         el (The DEM in meters)
#         watershed (A map defining the extent of the watershed 1=inside watershed)
#         residual_mc  =  residual moisture content (m^3/m^3)
#         porosity = soil porosity neglecting rock content (m^3/m^3)
#         fieldcap_mc =  field capacity moisture content (m^3/m^3)
#         rock_percent =  Volumetric rock content (%)
#         soil_depth =  Depth to a hydraulic restrictive layer (cm)
#         wiltpt_mc =  Wilting point moisture content (m^3/m^3)
#         ETreduction_mc = Moisture content at which actual ET becomes limited by soil moisture
#					(assumed to be 80% of field capacity) (m^3/m^3)
#         Ksat_matix =  The lateral matrix Ksat (ignoring macropores) of the soil layer (cm/day)
#	  Ksat_mpore =  The lateral Ksat of the soil layer when macropores are active (cm/day)
#         Ksubsurface = The vertical Ksat through the hydraulic restrictive layer (cm/day)

#    BEGIN PROGRAM
#    Any mask is removed so that values fill the entire region
print `g.remove -f type=raster name=MASK`;
print `r.slope.aspect elevation=el slope=slope --o`;   # Slope is in degree (default)

# -----------------------------------------------------------------------------
#    MOISTURE CONTENTS
#         residual.mc  =  residual moisture content (m^3/m^3)
#         porosity = soil porosity neglecting rock content (m^3/m^3)
#         fieldcap.mc =  field capacity moisture content (m^3/m^3)
#         rock.percent =  Volumetric rock content (%)
#         soil.depth =  Depth to a hydraulic restrictive layer (cm)
#         wiltpt.mc =  Wilting point moisture content (m^3/m^3)
#	  max_canopy_storage.amt = maximum storage of water in Spruce canopy (cm)
#         tree elevation is assumed at 20 meters

print `r.mapcalc 'residual_mc_A = 0.02' --o`;
print `r.mapcalc 'residual_mc_B = 0.02' --o`;

# Using wiltpt_mc_a and wiltpt_mc_b divided by 100 to calculate A and B --MZ Feb. 16 2016
print `r.mapcalc 'wiltpt_mc_A = wp_mc_a/100.0' --o`;
print `r.mapcalc 'wiltpt_mc_B = wp_mc_b/100.0' --o`;
#print `r.mapcalc 'wiltpt_mc_A = if(wp_mc_a>4.0,5.0,3.8)/100.0' --o`;# update: output 91; 98;
#print `r.mapcalc 'wiltpt_mc_B = if(wp_mc_b>2.0,1.0,0.9)/100.0' --o`;# update: output 91; 98;


#  strms_30m was created assuming a stream cells is defined by an upslope contributing area of 50 30x30m cells

print `r.mapcalc 'soil_depth_A = if(strms_30m>0,1.0,soil_depth_a)' --o`;  #output 91, 92, 93, 94, 95
#print `r.mapcalc 'soil_depth_A = if(strms_30m>0,1.0,soil_depth_a+20)' --o`;  #Using the soil_depth_a & soil_depth_b below --MZ
print `r.mapcalc 'soil_depth_B = if(strms_30m>0,19.0,soil_depth_b)' --o`;
#print `r.mapcalc 'soil_depth_B = if(strms_30m>0,19.0,soil_depth_b+200)' --o`; # add 300cm to enhao's soil b layer; 96 (+200, update +350)
#print `r.mapcalc 'soil_depth_B = if(strms_30m>0,19.0,soil_depth_b)' --o`; # output 91, 92, 93, 94, 95


#  increase the B horizon soil depth of the deepest soils by 125 cm to make a total soil depth of 3 m
#print `r.mapcalc soil_depth_B='if(soil_depth_B==122.0,soil_depth_B+125.0,soil_depth_B)'`;
#  increase the B horizon soil depth of the deepest soils by 225 cm to make a total soil depth of 4 m
#  increase the B horizon soil depth of the deepest soils by 325 cm to make a total soil depth of 5 m
#print `r.mapcalc soil_depth_B='if(soil_depth_B==122.0,soil_depth_B+325.0,soil_depth_B)'`;

#print `r.mapcalc 'soil_depth_B = if(soil_depth_B==122.0,soil_depth_B+325.0,soil_depth_B)' --o`;   # decrease from 450 to 225 to test if the soil depth  (4m in total) can increase the peak flow 12/26/2016 MZ; Update: commented on 20190424 becasue a new soil depth map from Enhao was added
# print `r.mapcalc 'soil_depth_B = if(soil_depth_B==122.0,soil_depth_B+450.0,soil_depth_B)' --o`;   #original value 12/26/2016 MZ total soil depth of 6.25m

print `r.mapcalc 'sat_mc_A = if(sat_mc_a<58.0,0.65,0.58)' --o`; # sat_mc_A and sat_mc_B are trial values 20190503 MZ; update: output 91;
print `r.mapcalc 'sat_mc_B = if(sat_mc_b>11.0,0.30,0.20)' --o`; # sat_mc_A and sat_mc_B are trial values 20190503 MZ; update: output 91;
#print `r.mapcalc 'sat_mc_A = sat_mc_a/100.0' --o`; # original 51.8%, 58%
#print `r.mapcalc 'sat_mc_B = sat_mc_b/100.0' --o`; # original 13.4%, 10%
print `r.mapcalc 'fieldcap_mc_A = if(fc_mc_a>38.0,0.30,0.15)' --o`; # fieldcap_mc_A and B are trial values 20190503 MZ
print `r.mapcalc 'fieldcap_mc_B = if(fc_mc_b>8.0,0.12,0.072)' --o`; # fieldcap_mc_A and B are trial values 20190503 MZ; update: output 91; output 92 ~
#print `r.mapcalc 'fieldcap_mc_A = if(fc_mc_a>18,17.3,18.8)/100.0' --o`;# update: output 91;
#print `r.mapcalc 'fieldcap_mc_A = fc_mc_a/100.0' --o`; # original 38%, 17.6%
#print `r.mapcalc 'fieldcap_mc_B = fc_mc_b/100.0' --o`; # original 3.7%, 1.6%
print `r.mapcalc 'soil_depth = soil_depth_A+soil_depth_B' --o`;
print `r.mapcalc 'sat_amt = sat_mc_A*soil_depth_A+sat_mc_B*soil_depth_B' --o`;
print `r.mapcalc 'fieldcap_amt_A = soil_depth_A*fieldcap_mc_A-soil_depth_A*residual_mc_A' --o`;
print `r.mapcalc 'fieldcap_amt_B = soil_depth_B*fieldcap_mc_B-soil_depth_B*residual_mc_B' --o`;
print `r.mapcalc 'fieldcap_amt = fieldcap_amt_A+fieldcap_amt_B' --o`;
print `r.mapcalc 'wiltpt_amt = wiltpt_mc_A*soil_depth_A+wiltpt_mc_B*soil_depth_B' --o`;

# ETreduction.mc = Moisture content at which actual ET becomes limited by soil moisture
#				(assumed to be 80% of field capacity) (m^3/m^3)
print `r.mapcalc 'ETreduction_mc = fieldcap_amt*0.8/soil_depth' --o`;


# ------------------------------------------------------------------------------
#    HYDRAULIC CONDUCTIVITIES (cm/day)
#     Ksat_matix =  The lateral matrix Ksat (ignoring macropores) of the soil layer (cm/day)
#     Ksat_mpore =  The lateral Ksat of the soil layer when macropores are active (cm/day)
#     Ksubsurface = The vertical Ksat through the hydraulic restrictive layer (cm/day)
#     Ksubsurface is assumed higher near the buildings because the soil layer appears artificial
#print `r.mapcalc 'Ksat_matrix_A = matrix_ks_a-64.2' --o`; # Output 91;
print `r.mapcalc 'Ksat_matrix_A = matrix_ks_a' --o`;
#print `r.mapcalc 'Ksat_matrix_B = matrix_ks_b-217.19' --o`; # Output 91
print `r.mapcalc 'Ksat_matrix_B = matrix_ks_b-150' --o`; # change to '-50' to make the lateral flow faster 20171220 MZ, and 20171222 MZ;before 95 = matrix_ks_B-150; output 95 (=mattrix_ks_b)
print `r.mapcalc 'Ksat_mpore_A = Ksat_matrix_A*10.0' --o`;
print `r.mapcalc 'Ksat_mpore_B = Ksat_matrix_B/2.0' --o`; # MZ 20200605
#print `r.mapcalc 'Ksubsurface = Ksat_matrix_B/1000.0' --o`;
print `r.mapcalc 'Ksubsurface = matrix_ks_a/500.0' --o`;

#    Kfc is the conductivity at field capacity calculated from Ksat and
#    moisture content relationships from Bresler's formula

print `r.mapcalc 'Kfc_A = Ksat_matrix_A*exp((-13.0/sat_mc_A)*(sat_mc_A-fieldcap_amt_A/soil_depth_A))' --o`;
print `r.mapcalc 'Kfc_B = Ksat_matrix_B*exp((-13.0/sat_mc_B)*(sat_mc_B-fieldcap_amt_B/soil_depth_B))' --o`;

# --------------------------------------------------------------------------------
#    FLOW DIRECTION PROGRAM
#    A "flow unit" is an elevation difference of one unit for an adjacent
#    cell.  An elevation difference of one unit for a diagonal cell would
#    be 1/(square root of 2)=0.707  flow units.


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

#print `r.mapcalc watershed='if(isnull(watershed),0.0,watershed)'

#print `rm -r /home/brooks/gisdata/mica_creek/mica_30m/.tmp`;
# print `rm -r /Users/apple/GrassWorkSpace/Mica_creek/.tmp`;
