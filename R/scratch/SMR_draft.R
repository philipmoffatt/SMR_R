# This is a rough draft the SMR as R functions with the long term to build an 
# R SMR package, allows for the seamless initialization, running, and reporting of a base model
# of a base model for any small watershed.  
#
# This is not is not in draft for use now. 
#
# author Ames Fowler andapdted from cripts by Erin Brooks and co. 
# fowler53@msu.edu 

DEM = rast("C:/Users/amesf/Dropbox/SMR_R/processed_data/ASC/MFC/dem_breached.asc")

cli_scale = (DEM/543.728)^.5  #dummy -make file 
# plot(DEM)

precip = 0#1.5# # dummy wrap up in Json parameter file. 
tavg = .6#C     # dummy wrap up in Json parameter file.


# Functions 1 --------------------
## setup function 1.1X
 ## write out the SMR_data_setup in functions for 1 getting data, 
 ## 2 processing DEMs, 3 processing soisl 

### 1.10 call asc into environment
load_asc_data <- function(Asc_dir = "C:/Users/amesf/Dropbox/SMR_R/processed_data/ASC/MFC/"){ # dummy make dynamic 
  ASC_files <- list.files(Asc_dir, full.names = T)
  temp_stack <- rast(ASC_files)
  return(temp_stack)
}
### 1.11 setup storage-----------
init_bucket_stack <- function(temp_stack){
  
  ## conceptual model 
  # ref_stack = spatial parameters and fluxes 
  # bucket_stack = Max Storage 
  # storage_stack = current water storage (depth) in each cell down through tops soil
  #
  # canopy snow
  # canopy water
  # snow
  # surface detention 
  # soil a ## dummy add dynamic layers and Piezo depth. Ask Matt for root model. 
  # soil b 
 
  bucket_stack <-  10 ##need a snow and water if we are to model it all here
  names(bucket_stack) = "canopy_storage"
  bucket_stack$snow = bucket_stack$canopy_storage*0 ## watch non max snow 
  bucket_stack$surface = temp_stack$NLCD_simple^0*5 #mm
  bucket_stack$soil_A = temp_stack$depth_A*temp_stack$wsatiated.r_A/100*.98
  bucket_stack$soil_B = temp_stack$depth_B*temp_stack$wsatiated.r_B/100*.98
  bucket_stack$soil_A_FC = temp_stack$depth_A*temp_stack$wthirdbar.r_A/100
  bucket_stack$soil_A_WP = temp_stack$depth_A*temp_stack$wfifteenbar.r_A/100
  bucket_stack$soil_B_FC = temp_stack$depth_B*temp_stack$wthirdbar.r_B/100
  bucket_stack$soil_B_WP = temp_stack$depth_B*temp_stack$wfifteenbar.r_B/100
  return(bucket_stack)
}

### 1.12 setup storage-----------
init_storage_stack <- function(bucket_stack, init_moisture){
  storage_stack <- ref_stack*init_moisture  ## cheap initial moisture 
                                            ## fine with spinnup 
  return(storage_stack)}

### 1.13 setup storage -----------
init_ref_stack <- function(Asc_dir){
  ASC_files <- list.files(Asc_dir, full.names = T)
  temp_stack <- rast(ASC_files)
  
  # Flow rates
  ref_stack <-  temp_stack$ksat.r_A
  names(ref_stack) = "soil_A_ksat"
  ref_stack$soil_A_ksat_macro = temp_stack$ksat.r_A*10 #needs citation 
  ref_stack$soil_A_ksat_dp = temp_stack$ksat.r_A/10 #needs citation
  ref_stack$soil_B_ksat = temp_stack$ksat.r_B
  ref_stack$soil_B_ksat_macro = temp_stack$ksat.r_B*10 #needs citation 
  ref_stack$soil_B_ksat_dp = temp_stack$ksat.r_B/10 #needs citation
  #add unit test to log file 
  
  #add calc for percent or do every thing in depth terms?? not now. 
  
  ##Add KFC - do all static calcs once. 
  #    Kfc is the conductivity at field capacity calculated from Ksat and
  #    moisture content relationships from Bresler's formula
  
  ref_stack$Kfc_A = ref_stack$soil_A_ksat*
    exp((-13.0/(ref_stack$soil_A/temp_stack$depth_A))*
          ((ref_stack$soil_A-ref_stack$soil_A_FC)/temp_stack$depth_A))
  
  ref_stack$Kfc_B = ref_stack$soil_B_ksat*
    exp((-13.0/(ref_stack$soil_B/temp_stack$depth_B))*
          ((ref_stack$soil_B-ref_stack$soil_B_FC)/temp_stack$depth_B))
  
  #    FLOW DIRECTION PROGRAM
  #    A "flow unit" is an elevation difference of one unit for an adjacent
  #    cell.  An elevation difference of one unit for a diagonal cell would
  #    be 1/(square root of 2)=0.707  flow units.
  el = temp_stack$dem_breached

  # Calculate the differences with the 8-neighbors cells 
  
  ## check direction is not flipped ####
  diff_n <- el - (terra::shift(el,0,-30) %>% terra::extend(.,el) %>% terra::crop(.,el))
  diff_ne <- (el - (terra::shift(el, -30, -30) %>% terra::extend(.,el) %>% terra::crop(.,el)))* 0.707
  diff_e <- el - (terra::shift(el, -30, 0) %>% terra::extend(.,el) %>% terra::crop(.,el))
  diff_se <- (el - (terra::shift(el, -30, 30) %>% terra::extend(.,el) %>% terra::crop(.,el)))* 0.707
  diff_s <- el - (terra::shift(el, 0,30) %>% terra::extend(.,el) %>% terra::crop(.,el))
  diff_sw <- (el - ((terra::shift(el, 30, 30) %>% terra::extend(.,el) %>% terra::crop(.,el))))* 0.707
  diff_w <- el - (terra::shift(el, 30, 0) %>% terra::extend(.,el) %>% terra::crop(.,el))
  diff_nw <- (el - (terra::shift(el, 30, -30) %>% terra::extend(.,el) %>% terra::crop(.,el)))* 0.707
  
  
  # Take the maximum with 0.0
  # 
  # diff_ne[diff_ne < 0] <- 0
  # diff_e[diff_e < 0] <- 0
  # diff_se[diff_se < 0] <- 0
  # diff_s[diff_s < 0] <- 0
  # diff_sw[diff_sw < 0] <- 0
  # diff_w[diff_w < 0] <- 0
  # diff_nw[diff_nw < 0] <- 0
  # 
  slope_delta <- rast(list(diff_n, diff_ne, diff_e, diff_se, diff_s, diff_sw, diff_w, diff_nw))
  names(slope_delta) <- c("diff_n", "diff_ne", "diff_e", "diff_se", "diff_s", "diff_sw", "diff_w", "diff_nw")
  slope_delta[slope_delta < 0] <- 0
  
  
  # Sum the differences
  flowunits <- sum(slope_delta)
  
  #    The following maps are the percent of the cell's neighbors in a direction
  #    that will flow to the cell.  For example, "north", is the percent of the
  #    lateral flow out of the cell to the nor:th that will flow to the cell.
  
  percent_flow = slope_delta/flowunits
  names(percent_flow) <- c("n", "ne", "e", "se", "s", "sw", "w", "nw")
  
  ref_stack = c(ref_stack, percent_flow)
  return(ref_stack)
  }
  
### 1.21 Run precipitation -----------
run_precip <- function(precip, tavg, precip_scale =cli_scale,
                       el=DEM, tmax_rain=2, tmin_snow = -1, 
                       storage_stack = storage_stack){
  #temp
  tavg = tavg-0.00531*(el-599) ## dummy fix and make dynamic
  

  #prcip partition
    #what ET model are we using what is the assumption for days with p
    # recipitation/high relative humidity sublimation??
  
  if(precip>0){
    precip = cli_scale*(precip)
    
    if(tmax_rain == tmin_snow){
      if(tavg>=tmax_rain){
        rain_percent = 1
      }else{
        rain_percent = (tavg-(tmin_snow))/(tmax_rain-(tmin_snow)) 
        rain_percent[rain_percent<=0] = 0
        rain_percent[rain_percent>=1] = 1
      }
      
      et = 0 
      rain = precip_scale*precip
    }
  }else{ ## we may not want to assume low/no et with precipt
      et = ((tavg>0)*tavg*1.2) ### dummy add several precipitation models. 
      rain = 0*precip_scale
  }
  
  #snow melt model 
  if(sum(storage_stack$snow)>0){
    
    melt = snow_melt*(snow_melt>storage_stack$snow)
    rain = rain +melt
  }
  
  atm_flux = rast(list("rain" = rain,
                       "snow" = precip-rain, 
                       "et" = et))
  # length(names(temp))
  # names(temp) = c("rain", "snow", "et")
  return(temp) 
}



### 1.22 Run vertical redistribution -----------

run_vert_flow(){
  
  ## canopy model
  
  ## infiltration model 
  
  ## redistribution 

}
### 1.23 Run lateral flow  -----------

run_precip <- function(,, piezo = F){  ## ********************************* FIX here *******************************
    #    Lateral flow is the amount of moisture that LEAVES each cell.  North,
    #    northeast, etc., are fraction of neighboring cells' outgoing lateral
    #    flow that flows towards the cell
    
    if(storage_amt<fieldcap_amt){
        effK = Ksat_matrix*exp((-13.0/sat_mc)*(sat_mc-storage_amt/soil_depth))
      }elseif(storage_amt>=sat_mc*soil_depth){
        effK = Ksat_mpore,(Ksat_mpore-Kfc)*
          (storage_amt/soil_depth-fieldcap_amt/soil_depth)/
          (sat_mc-fieldcap_amt/soil_depth)
      }else{
        effK = Kfc
      }
        
      
    print `r.mapcalc 'lateral_flow = min(effK*slope/100.0*soil_depth/$gridsize/100.0,storage_amt)
    print `r.mapcalc 'lateral_flow = if(flowunits==0.0,0.0,lateral_flow)
  
  
  if(piezo == F){
   #    Lateral flow is based on Darcy's Law, with gradient
    #    equal to land slope, and direction maps (north, northeast, etc)
    #    calculated from the elevation map in the program smr.setup. 
    
    
  }else(piezo == T){
    # Transmissivity for a single exponential Ksat profile
    # the scaled piezo is calculated to prevent over-prediction of Ksat in shallow soils
    # Ksat and transmissivity is in cm/day and cm^2/day in smr.setup program
    
    #print `r.mapcalc scaled_piezo='if(soil_depth>65.0,piezo,piezo/soil_depth*65.0)'`; # cm
    #print `r.mapcalc transmissivity='(Ks_s_mp/Ks_exp*(exp(Ks_exp*D_mp)-exp(Ks_exp*(D_mp-max(0.0,piezo-D_no_mp))))+min(D_no_mp,piezo)*Ks_d_mp)*$time_step/24.0'`; # cm^2/hr
    #print `r.mapcalc ksat_avg='if(scaled_piezo==0.0,0.0,transmissivity/scaled_piezo)'`; # cm/hr
    #print `echo "lateral_flow=min(storage_amt,ksat_avg*piezo/($gridsize*100.0)*slope/100.0)" | r.mapcalc`; # cm
    #print `echo "lateral_flow=min(storage_amt,transmissivity*slope/100.0/($gridsize*100.0))" | r.mapcalc`; # cm
    #print `r.mapcalc lateral_flow='if(flowunits==0.0,0.0,lateral_flow)'`; # cm
  }
  storage_amt = storage_amt - lateral_flow  + 
    (if(isnull(lateral_flow[-1,0]),0.0,lateral_flow[-1,0])*north + 
    if(isnull(lateral_flow[-1,1]),0.0,lateral_flow[-1,1])*northeast+
    if(isnull(lateral_flow[0,1]),0.0,lateral_flow[0,1])*east+
    if(isnull(lateral_flow[1,1]),0.0,lateral_flow[1,1])*southeast+
    if(isnull(lateral_flow[1,0]),0.0,lateral_flow[1,0])*south+
    if(isnull(lateral_flow[1,-1]),0.0,lateral_flow[1,-1])*southwest + 
    if(isnull(lateral_flow[0,-1]),0.0,lateral_flow[0,-1])*west +
    if(isnull(lateral_flow[-1,-1]),0.0,lateral_flow[-1,-1])*northwest)
  lateral_flow_$wateryear = lateral_flow_$wateryear+lateral_flow
  
}
   

### 1.24 Run percolation -----------

### 1.25 Run runoff -----------

### 

### 2 set up SMR

routing_stack <- init_routing_stack(Asc_dir)
temp <- init_routing_stack(Asc_dir)
build_stacks <- 

weatherdata <- read_csv("C:/Users/amesf/Dropbox/ames ascii/prw_2yr_weather.txt.csv",col_names =  F)
names(weatherdata) = c("wateryear", "date","year","tmin","tmax","tavg", "precip","pet_sntjohn", "cc_forest","cc_partial","cc_wheat","cc_grass","cc_bare","cc_pea", "output")

i=1

tic()

stoagestack$depth_A*1
toc()
