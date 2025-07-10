###############################'Graphci representation of multivariate GLMM

#######################'Packages
library(tidyverse) # ‘2.0.0’
library(broom.mixed) # ‘0.2.9.5’
library(patchwork) # ‘1.3.0’
library(glmmTMB)  # ‘1.9.9’
library(ggeffects) # ggeffects
rm(list=ls())

###########'GLMM
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

centers <- sapply(df %>% select_if(is.numeric), function(x) attr(scale(x, center = TRUE, scale = FALSE), 'scaled:center')) ## allows to have the value of variables before centering
df_glmm<-df%>%
  mutate( PROP_PP = as.character(PROP_PP), NB_TOT=as.character(NB_TOT),  NB_ALBO_F=as.character(NB_ALBO_F), MOIS=as.character(MOIS),LATITUDE = as.character(LATITUDE),LONGITUDE = as.character(LONGITUDE)) %>%
  mutate_if(is.numeric, ~scale(., center = TRUE, scale = FALSE)) %>%  # we center  to be able to compare the magnitudes (centering also helps with allowing easier interpretation of variable coefficients associated with different magnitudes, e.g. when one regressor is measured on a very small order, while another on a very large order.  )
  mutate(PROP_PP = as.numeric(PROP_PP), NB_TOT=as.integer(NB_TOT), NB_ALBO_F=as.integer(NB_ALBO_F), MOIS=as.integer(MOIS),
         LATITUDE = as.numeric(LATITUDE),LONGITUDE = as.numeric(LONGITUDE))%>%
  arrange(MOIS)%>%
  na.omit()


####'Full model 
mod_select<-glmmTMB( PROP_PP ~ lsm_c_area_mn_LCG_100_13 + lsm_c_area_mn_LCG_250_12 +      lsm_c_te_LCG_20_10 + RHMAX_collection + (1 | num_session),
                     data = df_glmm,weights = NB_TOT, family = binomial(link = "logit"))

glmm_coeffs_pvalues <- broom.mixed::tidy(mod_select, conf.int = TRUE, exponentiate = TRUE)

#### Plotting the partial dependence plots for the retained independent variables
predictors_pdp <- glmm_coeffs_pvalues %>% 
  filter(p.value < 0.05 &effect == "fixed" & term != "estimate" & term != "(Intercept)")
  
predictors_pdp <- unique(predictors_pdp$term)

pdps_plots = list()
i<-1

for(i in 1:length(predictors_pdp)){
  df<-df_glmm
  var<-predictors_pdp[i]
  eff<-ggpredict(mod_select, terms = paste0(var, "[all]"), bias_correction = T)
  df[,var]<-df[,var]+centers[var]
  df$x<-df[,var]
  eff$x<-eff$x+centers[var]
  eff$GLMM<-eff$predicted
  p<-ggplot() +
   # geom_line(color = "blue", alpha=0.5) +
    geom_line(data = eff, aes(x = x, y = GLMM, color = "GLMM Prediction"), size = 1)+
    geom_ribbon(data = eff, aes(x = x, ymin = conf.low, ymax = conf.high, color="95% CI"), alpha = 0.2)  +
    theme_minimal()+
    ylab("Probability of parity predicted")+
    xlab(var)+
    ylim(c(0,1))+
    geom_point(data = df %>% count(x, PROP_PP), 
               aes(x = x, y = PROP_PP, size = n, color = "Observations"), alpha = 0.8) +
    
    # Scales adjusted
    scale_color_manual(name = "Legend", 
                       values = c("GLMM Prediction" = "blue", "Observations" = "grey50", "95% CI"="grey80")) +
   # scale_fill_manual(name = "Legend", values = c("95% CI" = "grey80")) +
    scale_size_continuous(name = "Observations", limits = c(0, 10),  breaks = seq(0,10, by = 2)) +
  
    theme_minimal() +
    labs(y = "Proportion of parous predicted", x = var) +
    ylim(0, 1)
  pdps_plots[[i]] <- p
}

# plot all the PDPs in one single figure
glmm_pdps <- patchwork::wrap_plots(pdps_plots) + plot_layout(guides = "collect")


ggsave(filename = "data/plots/prop_parous/plot_PDPs_glmm_parous.pdf",plot =glmm_pdps, device = "pdf", width = 7, height = 8)

