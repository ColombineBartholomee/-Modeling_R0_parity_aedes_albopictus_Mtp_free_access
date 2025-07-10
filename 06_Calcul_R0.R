# This script calculates and simulates theoretical R0 values 
# for a DENV based on parameter inputs and for different human exposure scenario


##################################'Packages used
rm(list=ls())
library(tidyverse)  # ‘2.0.0’
library(dplyr)# ‘1.1.2’"
library(lubridate) #  ‘1.9.3’
library(mgcv)  #‘1.9.1’
library(purrr)  #‘1.0.1’
library(boot) # ‘1.3.31’

#################################'Function used
source("analysis/functions_R0.R")

##################################'Dataset
df<-read.csv("data/df_model.csv")%>%
  dplyr::mutate(NB_TOT=NB_ALBO_PP+NB_ALBO_F_NP, DATE_COLLECTE=parse_date(DATE_COLLECTE,"%Y-%m-%d"), MONTH=as.integer(month(DATE_COLLECTE)), 
                PROP_PP=ifelse(NB_TOT>0, NB_ALBO_PP/NB_TOT, 0), AREA=as.factor(AREA), ZONE=as.factor(ZONE))%>%
  relocate(NB_TOT, PROP_PP,.before = NB_ALBO_F_ND)%>%
  relocate(DATE_COLLECTE, MONTH,.after = DATE_POSE)%>%
  # filter(NB_ALBO_F_PP_NP>0)%>%
  mutate(across(
    contains("lsm_c_area_mn_LCG"), 
    ~ . / 0.0001
  ), JOUR=as.factor(JOUR),GL_20=BUFFER_20, GL_50=BUFFER_50, GL_100=BUFFER_100, GL_250=BUFFER_250)%>%
  mutate(MONTH=as.integer(MOIS))%>%
  arrange(MONTH)%>%
  dplyr::select(ID_COLLECTE, ID_PIEGE,num_session,LATITUDE,LONGITUDE,  DATE_POSE, DATE_COLLECTE,  DATE_HEURE_COLLECTE,MOIS,MONTH, AREA, ZONE,  NB_ALBO_F,PROP_PP, NB_ALBO_F_NP,  NB_ALBO_PP, NB_TOT, NB_ALBO_F_DIS ,
                TMEAN_collection)%>%
  na.omit()

## data of exposure
#exposure<-read.csv('data/aedes_human_exposure.csv')%>%
#  dplyr::select(HOUR, ID_PIEGE,HUMAN_EXPOSURE , HUMAN_EXPOSURE_ACT,ACTIVITY_TOTAL , ACTIVITY_RELATIVE       )

##################################'Calculate R0 and all intermediate computations for each observation
df_output <- df %>%
  dplyr::mutate(
    TMN = TMEAN_collection, # temperature
    ma = NB_ALBO_F_DIS , # daily abundance in a trap, proxy for m.a, the mosquito:human ratio or no of daily bites per human
    parity_rate = PROP_PP,
    #  HUMAN_EXPOSURE=1,
    # ACTIVITY=1,
    expo=1,
    b = map2_dbl("DENV", TMN, compute_b),
    c = map_dbl("DENV", compute_c), 
    r = map_dbl("DENV", compute_r),
    N = map2_dbl("DENV",TMN, compute_N),
    phi=compute_phi(),## human feeding indice, which depends also of the time of exposure of each host
    gc = map_dbl(TMN, compute_gc),
    a = phi/gc, ## tran et al. 2005 : a = host trophic preference/ gonotrophic cyle ; and according to Garret-Jones (1969) : human biting habit (a) = frequency of human biting (1/gc) * proportion of human blood meal
    p = map2_dbl(parity_rate, gc, compute_p),
    capacity = pmap_dbl(list(ma, a, p,expo, N), ~compute_capacity(..1, ..2, ..3, ..4,..5)),
    R0_basic = pmap_dbl(list(capacity, r, b, c), ~compute_R0(..1, ..2, ..3, ..4)),
    R0_poletti = pmap(list(a, b, c, r, N, p, expo, ma), compute_R0_poletti)
  ) %>%
  unnest_wider(R0_poletti, names_sep = "_") 

##################################'Calculate low and high values of R0 with simulation of parity rate and of ma using gams, 

### Calcul of the parity rate for 1000 simulations
df_parity<-compute_parity_rate_ic(df_output,  method = "area_mois",nsim=1000)%>%
  select(ID_COLLECTE, sim_id, parity_rate_sim)

### Calcul of ma for 1000 simulations
df_ma<-compute_ma_counts(df_output, method = "area_mois", nsim=1000)%>%
  dplyr::select(ID_COLLECTE, sim_id, ma_sim)

### Putting all the parameters together and to simulate different human exposure scenarios + calculate R0 and R0 poletti with raw data and different human exposure scenario

# Selection of hours of exposition 
#hours_list <- unique(exposure$HOUR)

#  exposure values
expo_values <- c(0.10, 0.20, 0.30)
df_sim <- map_dfr(expo_values, function(expo_val) {
  
  df_tmp <- df_parity%>%
    left_join(df_ma)%>%
    left_join(df_output)%>%
    mutate(
      expo = expo_val,  # adding the value of exposure
      capacity = pmap_dbl(list(ma, a, p,expo, N), ~compute_capacity(..1, ..2, ..3, ..4,..5)),
      R0_basic = pmap_dbl(list(capacity, r, b, c), ~compute_R0(..1, ..2, ..3, ..4)),
      R0_poletti =pmap(list(a, b, c, r, N, p, expo, ma), compute_R0_poletti),
      R0_poletti_r0 = map_dbl(R0_poletti, "r0"),
      R0_poletti_r0VH = map_dbl(R0_poletti, "r0VH"),
      R0_poletti_r0HV = map_dbl(R0_poletti, "r0HV")
    ) %>%
    select(-R0_poletti)%>%
    dplyr::mutate(
      p_sim = map2_dbl(parity_rate_sim, gc, compute_p),
      capacity_sim = pmap_dbl(list(ma_sim, a, p_sim, expo, N), 
                              ~compute_capacity(..1, ..2, ..3, ..4, ..5)),
      R0_basic_sim = pmap_dbl(list(capacity_sim, r, b, c), 
                              ~compute_R0(..1, ..2, ..3, ..4)),
      R0_poletti_sim =pmap(list(a, b, c, r, N, p_sim, expo, ma_sim), compute_R0_poletti),
      R0_poletti_sim_r0 = map_dbl(R0_poletti_sim, "r0"),
      R0_poletti_sim_r0VH = map_dbl(R0_poletti_sim, "r0VH"),
      R0_poletti_sim_r0HV = map_dbl(R0_poletti_sim, "r0HV")
    ) %>%
    select(-R0_poletti_sim)
  
  return(df_tmp)
})


### To group by ID_COLLECTE, and to calculate the Ic low and high
df_fin<-df_sim%>%
  group_by(ID_COLLECTE, ZONE, AREA, MOIS, ma, parity_rate, TMN,b, c, r, N, gc,phi, a, p,expo, capacity, R0_basic, R0_poletti_r0, R0_poletti_r0VH,R0_poletti_r0HV)%>%
  summarise(# Compute average values
    capacity = mean(capacity, na.rm = TRUE),
    R0_basic_sim_mean = mean(R0_basic_sim, na.rm = TRUE),
    R0_poletti_sim_mean = mean(R0_poletti_sim_r0, na.rm = TRUE),
    
    # Confidence intervals using empirical quantiles (2.5% and 97.5%)
    R0_basic_low_quantile = quantile(R0_basic_sim, 0.025, na.rm = TRUE),
    R0_basic_high_quantile = quantile(R0_basic_sim, 0.975, na.rm = TRUE),
    R0_poletti_low_quantile = quantile(R0_poletti_sim_r0, 0.025, na.rm = TRUE),
    R0_poletti_high_quantile = quantile(R0_poletti_sim_r0, 0.975, na.rm = TRUE)) 

#write.csv(df_sim, "data/df_sim.csv") ## to save simulations
write.csv(df_fin, "data/df_fin.csv") ## to save simulations averaged

### To calculate a R0 with a daily mean exposure
df_fin<-read.csv("data/df_fin.csv")%>%
  group_by(ID_COLLECTE, ZONE, AREA, MOIS, ma, parity_rate, TMN,b, c, r, N, gc,phi, a, p)%>%
  summarise(# Compute average values
    expo=mean(expo, na.rm=T), 
    capacity=mean(capacity, na.rm=T), 
    R0_basic=mean(R0_basic, na.rm=T), 
    R0_poletti_r0=mean(R0_poletti_r0, na.rm=T),
    R0_poletti_r0VH=mean(R0_poletti_r0VH, na.rm=T),
    R0_poletti_r0HV=mean(R0_poletti_r0HV, na.rm=T),
    capacity = mean(capacity, na.rm = TRUE),
    R0_basic_sim_mean = mean(R0_basic_sim_mean, na.rm = TRUE),
    R0_poletti_sim_mean = mean(R0_poletti_sim_mean, na.rm = TRUE),
    # Confidence intervals using empirical quantiles (2.5% and 97.5%)
    R0_basic_low_quantile = mean(R0_basic_low_quantile,  na.rm = TRUE),
    R0_basic_high_quantile = mean(R0_basic_high_quantile,  na.rm = TRUE),
    R0_poletti_low_quantile = mean(R0_poletti_low_quantile,  na.rm = TRUE),
    R0_poletti_high_quantile = mean(R0_poletti_high_quantile,  na.rm = TRUE)) 

write.csv(df_fin, "data/df_fin_averaged.csv")
