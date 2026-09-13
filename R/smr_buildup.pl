#!/usr/bin/perl
# Usage: perl smr_buildup.pl <weather_file> <output_file>
#   weather_file : path to space-delimited weather CSV (from wxData.R)
#   output_file  : path for mass-balance CSV output

use strict;
use warnings;

my $weather_file = $ARGV[0] or die "Usage: perl smr_buildup.pl <weather_file> <output_file>\n";
my $output_file  = $ARGV[1] or die "Usage: perl smr_buildup.pl <weather_file> <output_file>\n";
my $output_mode   = $ENV{SMR_OUTPUT_MODE} // 'full';
die "SMR_OUTPUT_MODE must be full, daily, or final-only (got '$output_mode')\n"
    unless $output_mode =~ /\A(?:full|daily|final-only)\z/;
my ($first_date, $last_date, $days_processed);
open(my $WSHED_TEMPLATE, '<', 'wshed_res_properties.ini') or die "Can't open wshed_res_properties.ini: $!";
my @watershed_records = grep { /\S/ } <$WSHED_TEMPLATE>;
close($WSHED_TEMPLATE);
my $single_watershed = @watershed_records == 1;

print `g.remove -f type=raster name=MASK`;
print `g.region rast=watershed`;
print `r.mapcalc "expression=watershed = if(watershed==1.0,1,1)" --overwrite`;
print `r.to.vect input=watershed output=watershed_extent type=area --overwrite --quiet`;
print `g.region vector=watershed_extent align=watershed`;
print `r.mask raster=watershed --overwrite --quiet`;

# Initialize variables
my $gridsize       = 30.0;  # m
my $time_step      = 24.0;  # hrs
my $temp_time_step = 6.0;   # hrs

# lapse rate variables (Bull Run site: ~900m elevation; Pullman station: 784m)
my $high_site_el               = 900;
my $low_site_el                = 784;
my $mean_annual_onsite_precip  = 130.0;   # cm (Bull Run ~130 cm/yr)
my $mean_annual_sea_level_precip = 15.54; # cm (from original calibration)
my $precip_LR                  = 0.09195; # cm / m
my $temp_LR                    = -3.20;   # C / km

# Land cover maps
print `r.mapcalc "expression=max_canopy_storage_amt = if(landuse==6.0,0.147,if(landuse==5.0,0.15,if(landuse==4.0,0.15,if(landuse==3.0,0.31,if(landuse==2.0,0.1,if(landuse==1.0,0.0,0.0))))))" --overwrite`;
print `r.mapcalc "expression=canopy_cover = if(landuse==6.0,0.50,if(landuse==5.0,0.60,if(landuse==4.0,0.60,if(landuse==3.0,0.80,if(landuse==2.0,0,if(landuse==1.0,0,0.0))))))" --overwrite`;
print `r.mapcalc "expression=root_zone = if(landuse==6.0,soil_depth,if(landuse==5.0,soil_depth_A,if(landuse==4.0,soil_depth,if(landuse==3.0,soil_depth,if(landuse==2.0,0.0,if(landuse==1.0,0.0,0.0))))))" --overwrite`;
print `r.mapcalc "expression=wiltpt_amt = if(landuse==6.0,wiltpt_mc_A*soil_depth_A+wiltpt_mc_B*soil_depth_B,if(landuse==5.0,wiltpt_mc_A*soil_depth_A,if(landuse==4.0,wiltpt_mc_A*soil_depth_A+wiltpt_mc_B*soil_depth_B,if(landuse==3.0,wiltpt_mc_A*soil_depth_A+wiltpt_mc_B*soil_depth_B,if(landuse==2.0,0.0,if(landuse==1.0,0.0,0.0))))))" --overwrite`;
print `r.mapcalc "expression=ETreduction_mc = if(landuse==6.0,fieldcap_amt*0.8/soil_depth,if(landuse==5.0,fieldcap_amt_A*0.8/soil_depth_A,if(landuse==4.0,fieldcap_amt*0.8/soil_depth,if(landuse==3.0,fieldcap_amt*0.7/soil_depth,if(landuse==2.0,0.0,if(landuse==1.0,0.0,0.0))))))" --overwrite`;
print `r.mapcalc "expression=tmax_rain = if(landuse==6.0,1.8,if(landuse==5.0,-0.5,if(landuse==4.0,1.8,if(landuse==3.0,3.1,if(landuse==2.0,-0.5,if(landuse==1.0,-0.5,0))))))" --overwrite`;
print `r.mapcalc "expression=tmin_snow = if(landuse==6.0,-1.0,if(landuse==5.0,-3.0,if(landuse==4.0,-1.0,if(landuse==3.0,0.9,if(landuse==2.0,-3.0,if(landuse==1.0,-3.0,0.0))))))" --overwrite`;

# Initial conditions
print `r.mapcalc "expression=storage_amt = sat_amt*0.4" --overwrite`;
print `r.mapcalc "expression=storage_amt_A = (wiltpt_mc_A+(fieldcap_mc_A-wiltpt_mc_A)*0.2)*soil_depth_A" --overwrite`;
print `r.mapcalc "expression=canopy_storage_amt = 0.0" --overwrite`;
print `r.mapcalc "expression=swe = 0.0" --overwrite`;
print `r.mapcalc "expression=snow.age = 1.0" --overwrite`;
print `r.mapcalc "expression=swe.yesterday = 0.0" --overwrite`;
print `r.mapcalc "expression=albedo = 0.2" --overwrite`;
print `r.mapcalc "expression=liquid.water = 0.0" --overwrite`;
print `r.mapcalc "expression=ice.content = 0.0" --overwrite`;
print `r.mapcalc "expression=tsnow_surf = -2.0" --overwrite`;
print `r.mapcalc "expression=u.surface = 24.0*tsnow_surf*(2.1*1000.0*min(swe/100.0,0.02)+1000.0*2.1*max(0.0,0.02-swe/100.0))" --overwrite`;
print `r.mapcalc "expression=tsnow.pack = -0.5" --overwrite`;
print `r.mapcalc "expression=u.total = 24.0*tsnow.pack*(2.1*1000.0*swe/100.0+2.1*1000.0*0.4)" --overwrite`;
print `r.mapcalc "expression=mass_balance_total = 0.0" --overwrite`;
print `r.mapcalc "expression=canopy_storage_amt_pre = 0.0" --overwrite`;
print `r.mapcalc "expression=storage_amt_pre = storage_amt" --overwrite`;
print `r.mapcalc "expression=max_storage_A_horizon = 0.0" --overwrite`;
print `r.mapcalc "expression=max_storage_B_horizon = 0.0" --overwrite`;
print `r.mapcalc "expression=max_storage = 0.0" --overwrite`;
print `r.mapcalc "expression=max_runoff = 0.0" --overwrite`;

# Compute slope/aspect inside GRASS
system('r.slope.aspect elevation=el slope=slope aspect=aspect --overwrite --quiet');

# Read watershed/reservoir properties
my %area_;
my %res_vol_;
my %res_coeff_;
my %base_flow_;

open(my $WATERSHED_INIT, '<', 'wshed_res_properties.ini') or die "Cannot open wshed_res_properties.ini: $!";
while (my $line = <$WATERSHED_INIT>) {
    chomp $line;
    my ($wshed_id, $area_cells, $res_vol, $res_coeff) = split(' ', $line);
    if ($wshed_id > 0) {
        $area_{$wshed_id}      = $area_cells * $gridsize * $gridsize;
        $res_vol_{$wshed_id}   = $res_vol;
        $res_coeff_{$wshed_id} = $res_coeff;
        $base_flow_{$wshed_id} = 0.0;
    }
}
close($WATERSHED_INIT);

# Annual accumulator maps — initialize for years in the run
# (will be re-used across years; just initialise once)
my %year_init_done;

# Weather loop
open(my $WEATHER, '<', $weather_file) or die "Cannot open weather file '$weather_file': $!";

while (<$WEATHER>) {
    chomp;
    my ($date,$year,$month,$day,$doy,$tmax,$tmin,$tavg,$tdew,$precip,$pet,
        $hour_1,$hour_6,$hour_12,$hour_18,$l_turb,$cloud,
        $cc_water,$cc_urban,$cc_forest,$cc_shrub,$cc_grass,$cc_row_crop,
        $rh_water,$rh_urban,$rh_forest,$rh_shrub,$rh_grass,$rh_row_crop,
        $output) = split(' ', $_);
    $first_date //= $date;
    $last_date = $date;
    $days_processed++;

    # init annual accumulator maps on first occurrence of each year
    unless ($year_init_done{$year}) {
        print `r.mapcalc "expression=runoff_$year = 0.0" --overwrite`;
        print `r.mapcalc "expression=runoff_feb_$year = 0.0" --overwrite`;
        print `r.mapcalc "expression=Psat_$year = 0.0" --overwrite`;
        $year_init_done{$year} = 1;
    }

    my @hourly_tmp_array = ($hour_1, $hour_6, $hour_12, $hour_18);

    print "\n|----- SNOW ACCUMULATION AND MELT MODEL -----|\n";

    # Precipitation lapse and partitioning
    print `r.mapcalc "expression=precip = (el*$precip_LR+$mean_annual_sea_level_precip) / ($mean_annual_onsite_precip) * $precip" --overwrite`;
    print `r.mapcalc "expression=tavg = $tavg + ($temp_LR/1000) * (el - $low_site_el)" --overwrite`;
    print `r.mapcalc "expression=tdew = $tdew + ($temp_LR/1000) * (el - $low_site_el - el)" --overwrite`;
    print `r.mapcalc "expression=rain = if(tavg>tmax_rain,precip,if(tavg<tmin_snow,0.0,(tavg-tmin_snow)/(tmax_rain-tmin_snow)*precip))" --overwrite`;
    print `r.mapcalc "expression=snow = precip-rain" --overwrite`;

    # PET lapse
    print `r.mapcalc "expression=pet_data = (1.8281*$pet-0.2486)" --overwrite`;
    print `r.mapcalc "expression=pet = (pet_data - $pet)*(el - $low_site_el)/($low_site_el - $high_site_el)+pet_data" --overwrite`;

    # Canopy interception
    print `r.mapcalc "expression=canopy_storage_amt = canopy_storage_amt + rain" --overwrite`;
    print `r.mapcalc "expression=throughfall = if(canopy_storage_amt>max_canopy_storage_amt,canopy_storage_amt-max_canopy_storage_amt,0.0)" --overwrite`;
    print `r.mapcalc "expression=canopy_storage_amt = canopy_storage_amt-throughfall" --overwrite`;
    print `r.mapcalc "expression=canopy_evap = if(rain>0.0,min(canopy_storage_amt,pet),min(canopy_storage_amt,pet))" --overwrite`;
    print `r.mapcalc "expression=pet = max(0.0,pet-canopy_evap)" --overwrite`;
    print `r.mapcalc "expression=canopy_storage_amt = canopy_storage_amt-canopy_evap" --overwrite`;

    # SAM — energy balance
    print `r.mapcalc "expression=rh = if(swe>0,$rh_urban/(1.0-(canopy_cover/3.0)),$rh_row_crop/(1.0-(canopy_cover/3.0)))" --overwrite`;
    print `r.mapcalc "expression=snow.age = if(snow>0.0 && throughfall==0.0,1.0,snow.age+1.0)" --overwrite`;
    print `r.mapcalc "expression=albedo = if(swe.yesterday+snow>0.0,min(0.95,0.7383*snow.age^(-0.1908)),0.2)" --overwrite`;

    print `r.mapcalc "expression=t_diff = if($cloud==1.0,0.1,0.75*(1-$cloud)*(1-exp(-0.6*$cloud/((0.75-0.4)*(1-$cloud))))) * watershed" --overwrite`;
    print `r.mapcalc "expression=coefdh = t_diff/0.10" --overwrite`;
    print `r.mapcalc "expression=coefbh = max(0.0,(1-$cloud)-t_diff/0.75)" --overwrite`;
    print `r.sun elevation=el aspect=aspect slope=slope linke_value=$l_turb albedo=albedo coeff_bh=coefbh coeff_dh=coefdh beam_rad=beam_rad diff_rad=diff_rad refl_rad=refl_rad day=$doy nprocs=4 --overwrite --quiet`;
    print `r.mapcalc "expression=q.srad = (1.0-canopy_cover)*(1.0-albedo)*(beam_rad+diff_rad+refl_rad) * 3600.0/1000.0" --overwrite`;
    print `r.mapcalc "expression=q.lw = 5.67*(10.0^-8.0)*3600.0*24.0/1000.0*(max(0.72,((9.2*(10^-6.0)*(tavg+273.15)^2.0)*(1.0-0.84*$cloud)+0.84*$cloud)*(1-canopy_cover)+0.92*canopy_cover)*(tavg+273.15)^4.0-0.98*(tsnow_surf+273.15)^4.0)" --overwrite`;

    print `r.mapcalc "expression=vap.d.snow_sat = exp((16.78*tsnow_surf-116.9)/(tsnow_surf+237.3))/(273.15+tsnow_surf)/0.4615" --overwrite`;
    print `r.mapcalc "expression=vap.d.snow = if(tsnow_surf<0.0,vap.d.snow_sat*exp(335.0*tsnow_surf/(0.4615*(273.15+tsnow_surf)^2.0)),vap.d.snow_sat)" --overwrite`;
    print `r.mapcalc "expression=vap.d.air = exp((16.78*tdew-116.9)/(tdew+237.3))/(tdew+273.15)/0.4615" --overwrite`;

    print `r.mapcalc "expression=q.latent = (2500.0+335.0)*(vap.d.air-vap.d.snow)/(rh/(3600.0*24.0))" --overwrite`;
    print `r.mapcalc "expression=condens = if((swe.yesterday+snow)>0.0,q.latent/((2500.0+335.0)*1000.0)*100.0,0.0)" --overwrite`;
    print `r.mapcalc "expression=q.sensible = 1.0*1.29*(tavg-tsnow_surf)/(rh/(3600.0*24.0))" --overwrite`;
    print `r.mapcalc "expression=q.rain.ground = 4.2*1000.0*(throughfall/100.0)*(max(0.0,tavg)-tsnow_surf)+2.1*1000.0*(snow/100.0)*(min(0,tavg)-tsnow_surf)+0.0" --overwrite`;
    print `r.mapcalc "expression=q.total = q.srad+q.lw+q.sensible+q.latent+q.rain.ground" --overwrite`;

    print `r.mapcalc "expression=refreeze = min((throughfall+liquid.water),max(0.0,-1.0*100.0*(u.total+q.total)/(335.0*1000.0)))-min((snow+ice.content+condens),max(0.0,100.0*(u.total+q.total)/(335.0*1000.0)))" --overwrite`;
    print `r.mapcalc "expression=u.total = u.total+q.total+refreeze/100.0*335.0*1000.0" --overwrite`;
    print `r.mapcalc "expression=ice.content = if(u.total>0.000001,0.0,ice.content+snow+refreeze+condens)" --overwrite`;
    print `r.mapcalc "expression=liquid.water = if(u.total>0.000001||u.total<-0.000001,0.0,min(0.03/(1.0-0.03)*ice.content,liquid.water+throughfall-refreeze))" --overwrite`;
    print `r.mapcalc "expression=swe = liquid.water+ice.content" --overwrite`;
    print `r.mapcalc "expression=snowmelt = (swe.yesterday - swe) + snow+condens" --overwrite`;
    print `r.mapcalc "expression=swe.y.save = swe.yesterday" --overwrite`;
    print `r.mapcalc "expression=swe.yesterday = swe" --overwrite`;
    print `r.mapcalc "expression=u.surface = if(u.total==0.0,0.0,if(swe>0.0,min(0.0,u.surface+(q.total+refreeze/100.0*335.0*1000.0)*(min(swe/100.0,0.02)/0.05+max(0.02-swe/100.0,0.0)/0.3)),u.surface+(q.total+refreeze/100.0*335.0*1000.0)*(min(swe/100.0,0.02)/0.05+max(0.02-swe/100.0,0.0)/0.3)))" --overwrite`;
    print `r.mapcalc "expression=tsnow_surf = u.surface/24.0/(2.1*1000.0*min(swe/100.0,0.02)+1000.0*2.1*max(0.0,0.02-swe/100.0))" --overwrite`;

    # Road runoff
    print `r.mapcalc "expression=road_runoff = if(landuse==2.0,max(0.0,snowmelt+throughfall-1.5)/$gridsize,0.0)" --overwrite`;
    print `r.mapcalc "expression=water_input = snowmelt+throughfall" --overwrite`;
    print `r.mapcalc "expression=snowmelt = if(water_input>0.0,snowmelt*(water_input-road_runoff)/water_input,snowmelt)" --overwrite`;
    print `r.mapcalc "expression=throughfall = if(water_input>0.0,throughfall*(water_input-road_runoff)/water_input,throughfall)" --overwrite`;

    # Daily accumulators
    print `r.mapcalc "expression=temp_sum = 0.0" --overwrite`;
    print `r.mapcalc "expression=actualET_daily_flow = 0.0" --overwrite`;
    print `r.mapcalc "expression=perc_daily_flow = 0.0" --overwrite`;
    print `r.mapcalc "expression=runoff_daily_flow = 0.0" --overwrite`;

    # 6-hour sub-step loop
    foreach my $hrly_tmp (@hourly_tmp_array) {
        print `r.mapcalc "expression=temp = $hrly_tmp-($temp_LR/1000)*($low_site_el - el)" --overwrite`;
        print `r.mapcalc "expression=temp_sum = max(0.0,temp)+temp_sum" --overwrite`;
    }

    foreach my $hrly_tmp (@hourly_tmp_array) {
        print `r.mapcalc "expression=temp = $hrly_tmp-($temp_LR/1000)*($low_site_el - el)" --overwrite`;

        # Water distribution to soil layers
        print `r.mapcalc "expression=input = snowmelt*max(temp,0.0)/if(temp_sum==0.0,1.0,temp_sum)+throughfall*$temp_time_step/24.0" --overwrite`;
        print `r.mapcalc "expression=storage_amt = storage_amt+input" --overwrite`;
        print `r.mapcalc "expression=storage_amt_A_tmp = if(storage_amt<fieldcap_amt,min(input+storage_amt_A,fieldcap_amt_A),storage_amt_A)" --overwrite`;
        print `r.mapcalc "expression=storage_amt_B = if(storage_amt<fieldcap_amt,storage_amt-storage_amt_A_tmp,if((storage_amt-fieldcap_amt_A)<sat_mc_B*soil_depth_B,storage_amt-fieldcap_amt_A,sat_mc_B*soil_depth_B))" --overwrite`;
        print `r.mapcalc "expression=storage_amt_A = if(storage_amt<fieldcap_amt,storage_amt_A_tmp,if((storage_amt-storage_amt_B)<sat_mc_A*soil_depth_A,storage_amt-storage_amt_B,sat_mc_A*soil_depth_A))" --overwrite`;

        # Lateral flow (Darcy + Bresler)
        print `r.mapcalc "expression=effK_A = if(storage_amt_A<fieldcap_amt_A,Ksat_matrix_A*exp((-13.0/sat_mc_A)*(sat_mc_A-storage_amt_A/soil_depth_A)),if(storage_amt_A>=sat_mc_A*soil_depth_A,Ksat_mpore_A,(Ksat_mpore_A-Kfc_A)*(storage_amt_A/soil_depth_A-fieldcap_amt_A/soil_depth_A)/(sat_mc_A-fieldcap_amt_A/soil_depth_A)+Kfc_A))*$temp_time_step/24.0" --overwrite`;
        print `r.mapcalc "expression=effK_B = if(storage_amt_B<fieldcap_amt_B,Ksat_matrix_B*exp((-13.0/sat_mc_B)*(sat_mc_B-storage_amt_B/soil_depth_B)),if(storage_amt_B>=sat_mc_B*soil_depth_B,Ksat_mpore_B,(Ksat_mpore_B-Kfc_B)*(storage_amt_B/soil_depth_B-fieldcap_amt_B/soil_depth_B)/(sat_mc_B-fieldcap_amt_B/soil_depth_B)+Kfc_B))*$temp_time_step/24.0" --overwrite`;

        print `r.mapcalc "expression=lateral_flow_A = min(effK_A*slope/100.0*soil_depth_A/$gridsize/100.0,storage_amt_A)" --overwrite`;
        print `r.mapcalc "expression=lateral_flow_A = if(flowunits==0.0,0.0,lateral_flow_A)" --overwrite`;
        print `r.mapcalc "expression=lateral_flow_B = min(effK_B*slope/100.0*soil_depth_B/$gridsize/100.0,storage_amt_B)" --overwrite`;
        print `r.mapcalc "expression=lateral_flow_B = if(flowunits==0.0,0.0,lateral_flow_B)" --overwrite`;
        print `r.mapcalc "expression=lateral_flow = lateral_flow_A + lateral_flow_B" --overwrite`;

        print `r.mapcalc "expression=storage_amt_A = storage_amt_A - lateral_flow_A  + (if(isnull(lateral_flow_A[-1,0]),0.0,lateral_flow_A[-1,0])*north + if(isnull(lateral_flow_A[-1,1]),0.0,lateral_flow_A[-1,1])*northeast+if(isnull(lateral_flow_A[0,1]),0.0,lateral_flow_A[0,1])*east+if(isnull(lateral_flow_A[1,1]),0.0,lateral_flow_A[1,1])*southeast+if(isnull(lateral_flow_A[1,0]),0.0,lateral_flow_A[1,0])*south+if(isnull(lateral_flow_A[1,-1]),0.0,lateral_flow_A[1,-1])*southwest + if(isnull(lateral_flow_A[0,-1]),0.0,lateral_flow_A[0,-1])*west +if(isnull(lateral_flow_A[-1,-1]),0.0,lateral_flow_A[-1,-1])*northwest)" --overwrite`;
        print `r.mapcalc "expression=storage_amt_B = storage_amt_B - lateral_flow_B  + (if(isnull(lateral_flow_B[-1,0]),0.0,lateral_flow_B[-1,0])*north + if(isnull(lateral_flow_B[-1,1]),0.0,lateral_flow_B[-1,1])*northeast+if(isnull(lateral_flow_B[0,1]),0.0,lateral_flow_B[0,1])*east+if(isnull(lateral_flow_B[1,1]),0.0,lateral_flow_B[1,1])*southeast+if(isnull(lateral_flow_B[1,0]),0.0,lateral_flow_B[1,0])*south+if(isnull(lateral_flow_B[1,-1]),0.0,lateral_flow_B[1,-1])*southwest + if(isnull(lateral_flow_B[0,-1]),0.0,lateral_flow_B[0,-1])*west +if(isnull(lateral_flow_B[-1,-1]),0.0,lateral_flow_B[-1,-1])*northwest)" --overwrite`;
        print `r.mapcalc "expression=storage_amt = storage_amt - lateral_flow  + (if(isnull(lateral_flow[-1,0]),0.0,lateral_flow[-1,0])*north + if(isnull(lateral_flow[-1,1]),0.0,lateral_flow[-1,1])*northeast+if(isnull(lateral_flow[0,1]),0.0,lateral_flow[0,1])*east+if(isnull(lateral_flow[1,1]),0.0,lateral_flow[1,1])*southeast+if(isnull(lateral_flow[1,0]),0.0,lateral_flow[1,0])*south+if(isnull(lateral_flow[1,-1]),0.0,lateral_flow[1,-1])*southwest + if(isnull(lateral_flow[0,-1]),0.0,lateral_flow[0,-1])*west +if(isnull(lateral_flow[-1,-1]),0.0,lateral_flow[-1,-1])*northwest)" --overwrite`;

        # Re-distribute between A and B
        print `r.mapcalc "expression=storage_amt_A_tmp = if(storage_amt<fieldcap_amt,min(storage_amt_A,fieldcap_amt_A),storage_amt_A)" --overwrite`;
        print `r.mapcalc "expression=storage_amt_B = if(storage_amt<fieldcap_amt,storage_amt-storage_amt_A_tmp,if((storage_amt-fieldcap_amt_A)<sat_mc_B*soil_depth_B,storage_amt-fieldcap_amt_A,sat_mc_B*soil_depth_B))" --overwrite`;
        print `r.mapcalc "expression=storage_amt_A = if(storage_amt<fieldcap_amt,storage_amt_A_tmp,if((storage_amt-storage_amt_B)<sat_mc_A*soil_depth_A,storage_amt-storage_amt_B,sat_mc_A*soil_depth_A))" --overwrite`;

        # ET (Thornthwaite-Mather)
        print `r.mapcalc "expression=ET_coeff = if(landuse==6.0,$cc_row_crop,if(landuse==5.0,$cc_grass,if(landuse==4.0,$cc_shrub,if(landuse==3.0,$cc_forest,if(landuse==2.0,$cc_urban,if(landuse==1.0,$cc_water,$cc_water))))))" --overwrite`;
        print `r.mapcalc "expression=root_storage_amt = if(landuse==6.0,storage_amt,if(landuse==5.0,storage_amt,if(landuse==4.0,storage_amt,if(landuse==3.0,storage_amt_A,if(landuse==2.0,storage_amt,if(landuse==1.0,storage_amt,storage_amt))))))" --overwrite`;
        print `r.mapcalc "expression=actualET_flow = if(root_storage_amt>(root_zone*ETreduction_mc),min(max(0.0,(root_storage_amt-wiltpt_amt)),pet*ET_coeff*$temp_time_step/24.0),if(root_storage_amt>wiltpt_amt,min(max(0.0,(root_storage_amt-wiltpt_amt)),max(0.0,pet*((root_storage_amt/root_zone-wiltpt_amt/root_zone)/(ETreduction_mc-wiltpt_amt/root_zone))*ET_coeff*$temp_time_step/24.0)),0.0))" --overwrite`;
        print `r.mapcalc "expression=actualET_daily_flow = actualET_daily_flow+actualET_flow" --overwrite`;
        print `r.mapcalc "expression=storage_amt = storage_amt-actualET_flow" --overwrite`;
        print `r.mapcalc "expression=storage_amt_A = if(landuse==6.0,if(storage_amt<fieldcap_amt,fieldcap_amt_A*storage_amt/fieldcap_amt,fieldcap_amt_A),if(landuse==5.0,if(storage_amt<fieldcap_amt,fieldcap_amt_A*storage_amt/fieldcap_amt,fieldcap_amt_A),if(landuse==4.0,if(storage_amt<fieldcap_amt,fieldcap_amt_A*storage_amt/fieldcap_amt,fieldcap_amt_A),if(landuse==3.0,storage_amt_A-actualET_flow,if(landuse==2.0,if(storage_amt<fieldcap_amt,fieldcap_amt_A*storage_amt/fieldcap_amt,fieldcap_amt_A),if(landuse==1.0,if(storage_amt<fieldcap_amt,fieldcap_amt_A*storage_amt/fieldcap_amt,fieldcap_amt_A),0))))))" --overwrite`;
        print `r.mapcalc "expression=storage_amt_B = storage_amt-storage_amt_A" --overwrite`;

        # Percolation
        print `r.mapcalc "expression=perc = min(min(max(storage_amt-fieldcap_amt,0.0),storage_amt_B),Ksubsurface*$temp_time_step/24.0)" --overwrite`;
        print `r.mapcalc "expression=storage_amt = storage_amt - perc" --overwrite`;
        print `r.mapcalc "expression=storage_amt_B = storage_amt_B - perc" --overwrite`;
        print `r.mapcalc "expression=storage_amt_A_tmp = if(storage_amt<fieldcap_amt,min(storage_amt_A,fieldcap_amt_A),storage_amt_A)" --overwrite`;
        print `r.mapcalc "expression=storage_amt_B = if(storage_amt<fieldcap_amt,storage_amt-storage_amt_A_tmp,if((storage_amt-fieldcap_amt_A)<sat_mc_B*soil_depth_B,storage_amt-fieldcap_amt_A,sat_mc_B*soil_depth_B))" --overwrite`;
        print `r.mapcalc "expression=storage_amt_A = if(storage_amt<fieldcap_amt,storage_amt_A_tmp,if((storage_amt-storage_amt_B)<sat_mc_A*soil_depth_A,storage_amt-storage_amt_B,sat_mc_A*soil_depth_A))" --overwrite`;
        print `r.mapcalc "expression=perc_daily_flow = perc_daily_flow + perc" --overwrite`;

        # Saturation-excess runoff
        print `r.mapcalc "expression=runoff_flow = if(storage_amt>sat_amt,storage_amt-sat_amt,0.0)" --overwrite`;
        print `r.mapcalc "expression=storage_amt = storage_amt-runoff_flow" --overwrite`;
        print `r.mapcalc "expression=runoff_daily_flow = runoff_daily_flow+runoff_flow" --overwrite`;
        print `r.mapcalc "expression=saturation = storage_amt/sat_amt*100.0" --overwrite`;
        print `r.mapcalc "expression=runoff_$year = runoff_flow+runoff_$year" --overwrite`;
    }

    print `r.mapcalc "expression=runoff_daily_flow = runoff_daily_flow+road_runoff" --overwrite`;
    print `r.mapcalc "expression=Psat_$year = Psat_$year + if(saturation>=100,1,0)" --overwrite`;

    # February output maps
    if ($month == 2) {
        `r.mapcalc "expression=runoff_feb_$year = runoff_feb_$year + runoff_daily_flow" --overwrite`;
    }
    if ($doy == 61) {
        print `r.mapcalc "expression=average_feb_runoff_$year = runoff_feb_$year / 29" --overwrite`;
        print `r.out.ascii input=average_feb_runoff_$year output=feb_avg_runoff_$year`;
    }

    # Daily reporting is intentionally separate from state evolution.  The
    # final-only path preserves final maps while avoiding all zonal aggregation.
    if ($output_mode ne 'final-only') {
    for my $watershed_record (@watershed_records) {
        chomp $watershed_record;
        my ($wshed_id, $area_cells, $res_vol, $res_coeff) = split(' ', $watershed_record);

        my %out;
        my @report_vars = $output_mode eq 'full'
            ? qw(runoff_daily_flow precip rain actualET_daily_flow canopy_evap
                 snowmelt storage_amt throughfall canopy_storage_amt perc_daily_flow
                 swe condens snow q.srad q.latent q.sensible q.lw q.rain.ground q.total)
            : qw(runoff_daily_flow precip actualET_daily_flow storage_amt perc_daily_flow);
        for my $var (@report_vars) {
            if ($single_watershed) {
                # The active watershed MASK already limits the map.  Direct
                # r.univar avoids a temporary zonal raster and second process.
                my $summary = `r.univar -g map=$var`;
                my ($sum) = $summary =~ /^sum=([^\r\n]+)/m;
                die "r.univar returned no sum for $var\n" unless defined $sum;
                $out{$var} = $sum * 1;
            } else {
                my $safe = $var; $safe =~ s/\./_/g;
                print `r.stats.zonal base=watershed cover=$var output=tmp_$safe method=sum --overwrite --quiet`;
                my $val = `r.stats -A -n -N input=tmp_$safe`;
                chomp $val;
                $out{$var} = $val * 1;
            }
        }

        my $runoff_cms = ($out{runoff_daily_flow} / 100) * ($gridsize * $gridsize) / ($time_step * 3600.0);
        my $perc_cms   = $out{perc_daily_flow}    / 100.0 / ($time_step * 3600.0);

        $res_vol_{$wshed_id}  += $perc_cms - $base_flow_{$wshed_id};
        $base_flow_{$wshed_id} = $res_coeff_{$wshed_id} * $res_vol_{$wshed_id};
        my $Q = $base_flow_{$wshed_id} + $runoff_cms;

        open(my $OUT, '>>', $output_file) or die "Cannot open output file '$output_file': $!";
        if ($output_mode eq 'full') {
            print $OUT join(" ",
                $wshed_id, $date, $year,
                $out{runoff_daily_flow}, $out{precip}, $out{rain},
                $out{actualET_daily_flow}, $out{canopy_evap}, $out{snowmelt},
                $out{storage_amt}, $out{throughfall}, $out{canopy_storage_amt},
                $out{perc_daily_flow}, $Q,
                $out{swe}, $out{condens}, $out{snow},
                $base_flow_{$wshed_id},
                $out{'q.srad'}, $out{'q.latent'}, $out{'q.sensible'}, $out{'q.lw'},
                $out{'q.rain.ground'}, $out{'q.total'}
            ) . "\n";
        } else {
            print $OUT join(" ", $wshed_id, $date, $year,
                $out{runoff_daily_flow}, $out{precip}, $out{actualET_daily_flow},
                $out{storage_amt}, $out{perc_daily_flow}, $Q) . "\n";
        }
        close($OUT);
    }
    }

    # These maps are daily intermediates: reporting above has consumed every
    # one, and none is a next-day state.  Removing them bounds mapset growth
    # without changing the ordered state transitions retained below.
    my @daily_ephemeral_maps = qw(
        precip rain snow pet_data pet temp temp_sum input storage_amt_A_tmp
        effK_A effK_B lateral_flow_A lateral_flow_B lateral_flow
        actualET_flow actualET_daily_flow perc perc_daily_flow runoff_flow runoff_daily_flow
        canopy_evap throughfall water_input road_runoff
        rh t_diff coefdh coefbh beam_rad diff_rad refl_rad
        q.srad q.lw q.latent q.sensible q.rain.ground q.total
        vap.d.snow_sat vap.d.snow vap.d.air condens snowmelt refreeze
    );
    print `g.remove -f type=raster name=@{[join(',', @daily_ephemeral_maps)]} --quiet`;
}
close($WEATHER);

if ($output_mode eq 'final-only') {
    open(my $SUMMARY, '>', "$output_file.summary.txt")
        or die "Cannot write final-only summary '$output_file.summary.txt': $!";
    print $SUMMARY "output_mode=final-only\n";
    print $SUMMARY "days_processed=$days_processed\n";
    print $SUMMARY "first_date=$first_date\n";
    print $SUMMARY "last_date=$last_date\n";
    print $SUMMARY "final_state_maps=storage_amt,saturation\n";
    close($SUMMARY);
}
