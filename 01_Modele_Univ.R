###############################'Univariate modeling

# This script runs univariate GLMMs to test the effect of different explanatory variables
# on mosquito parity rate

################# Open packages
rm(list=ls())
library(tidyverse)# ‘2.0.0’
library(glmmTMB) # ‘1.9.9’
library(purrr)# ‘1.0.1’
library(furrr)# ‘0.3.1’
library(correlation)# 0.8.6’
library(performance)# ‘0.13.0’
library(minqa)# ‘1.2.8’
tolerance = 1e-12

################# Open dataset containing the dependant and independent variables

df<-read.csv("data/df_model.csv")%>%
  dplyr::filter(NB_ALBO_F_DIS>0)%>%
  dplyr::mutate(NB_TOT=NB_ALBO_PP+NB_ALBO_F_NP, PROP_PP=NB_ALBO_PP/NB_TOT, num_session=as.factor(num_session))%>% ## Creation of the variable proportion of females parous
  relocate(BUFFER_20, BUFFER_50, BUFFER_100, BUFFER_250, .after = ZONE )%>%
  mutate(across(
    contains("lsm_c_area_mn_LCG"), 
    ~ . / 0.0001))%>% ## Putting in m2 the surfaces which are in ha
  relocate(DATE_COLLECTE, .after=ID_COLLECTE)%>%
  relocate(NB_TOT, PROP_PP,.after=NB_ALBO_GRAVIDE)%>%
  relocate(MOIS, JOUR, .before=ZONE)%>%
 dplyr::filter(NB_TOT>0)%>%
  na.omit()


predictors <- setdiff(colnames(df), c("X","X.1" , "ID_COLLECTE","NB_ALBO_DISS","NB_ALBO_ND","ID_PIEGE", "DATE_COLLECTE","HEURE_COLLECTE","DATE_HEURE_COLLECTE","DATE_POSE", "DATE_HEURE_POSE","TYPE_PIEGE" ,             
                                         "LATITUDE","LONGITUDE","AREA","num_session","MOIS","JOUR", "NB_ALBO_F","NB_ALBO_F_DIS","NB_ALBO_F_ND","NB_ALBO_F_NP","NB_ALBO_PP","NB_ALBO_GRAVIDE",         
                                       "NB_TOT","PROP_PP"))


#############################
######### univariate modeling
#############################

###### GLMM #####
df_glmm<-df%>%
  dplyr::select(PROP_PP, NB_TOT,ID_PIEGE,AREA, num_session, MOIS, NB_ALBO_F,predictors) %>% 
  mutate( PROP_PP = as.character(PROP_PP), NB_TOT=as.character(NB_TOT),  NB_ALBO_F=as.character(NB_ALBO_F), MOIS=as.character(MOIS)) %>%
  mutate_if(is.numeric, ~scale(., center = TRUE, scale = FALSE)) %>%  # we center  to be able to compare the magnitudes (centering also helps with allowing easier interpretation of variable coefficients associated with different magnitudes, e.g. when one regressor is measured on a very small order, while another on a very large order.  )
  mutate(PROP_PP = as.numeric(PROP_PP), NB_TOT=as.integer(NB_TOT), NB_ALBO_F=as.integer(NB_ALBO_F), MOIS=as.integer(MOIS))%>%
  arrange(MOIS)%>%
  na.omit()


mod_nul<-glmmTMB(PROP_PP ~ 1 , ## null model 
                    data = df_glmm, family = binomial(link = "logit"), weight=NB_TOT)

mod_nul_al<-glmmTMB(PROP_PP ~ 1 +(1|num_session), ## null model with random effects
                 data = df_glmm, family = binomial(link = "logit"), weight=NB_TOT)


func <- function(x){ ## function which calculates univariate modelling 
  ret <- glmmTMB(as.formula(paste0("PROP_PP ~ ",x,"+ (1|num_session)")), data = df_glmm, family = binomial(link = "logit"), weight=NB_TOT)
  return(ret)
}

possible_a <- possibly(func, otherwise = NA_real_)

glmms_univs <- future_map(colnames(df_glmm[7:ncol(df_glmm)]), possible_a) ## Apply function func to all explicatives variables of the glmmm

glmm_to_rm <- NULL
glmm_to_rm2 <- NULL

for(i in 1:length(glmms_univs)){ ## to delete empty models
  ifelse(is.na(glmms_univs[[i]]), glmm_to_rm <- c(glmm_to_rm,i),
         ifelse(is.na(summary(glmms_univs[[i]])$AICtab[1]),c(glmm_to_rm2,i),glmms_univs[[i]]))
}

glmm_to_rm <- c(glmm_to_rm,glmm_to_rm2)
glmms_univs <- glmms_univs[-glmm_to_rm]

func2 <- function(x){ ## Fucntion which cerates a dataframe with the effects of each variable
  ret <- broom.mixed::tidy(x, conf.int = TRUE,  exponentiate = TRUE)
  ret$r2 <- performance::r2_nakagawa(x, tolerance = tolerance, null_model = mod_nul_al)$R2_marginal
  return(ret)
}

possible_b <- possibly(func2, otherwise = NULL)
glmm_univ_prop <- future_map(glmms_univs, possible_b) ## apply func2 to all the models

glmm_univ_prop_2<-do.call(rbind.data.frame, glmm_univ_prop) %>%
  filter(effect == "fixed" & term!="(Intercept)")

write.csv(glmm_univ_prop_2,"data/glmm_prop_parous_traps.csv", row.names = F)

