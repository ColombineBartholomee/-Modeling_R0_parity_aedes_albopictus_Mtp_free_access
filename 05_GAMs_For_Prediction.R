################################################################################
# Description:
# This script models:
#   - the proportion of parous females,
#   - and the number of females caught,
# using Generalized Additive Models (GAMs), stratified by ZONE and MONTH or AREA and MONTH.
# 
# Only the retained GAMs for both responses are presented.
################################################################################

##########################
# Load Required Packages #
##########################
library(ggplot2) #  ‘3.5.1’
library(tidyverse) # ‘2.0.0’
library(dplyr)  # ‘1.1.2’"
library(mgcv) #‘1.9.1’

####################
# Load and Prepare #
####################
df<-read.csv("data/df_model.csv")%>%
  dplyr::filter(NB_ALBO_F_DIS>0)%>%
  dplyr::mutate(NB_TOT=NB_ALBO_PP+NB_ALBO_F_NP, PROP_PP=ifelse(NB_TOT>0, NB_ALBO_PP/NB_TOT, 0), num_session=as.factor(num_session))%>% ## Creation of the variable proportion of females parous
  relocate(BUFFER_20, BUFFER_50, BUFFER_100, BUFFER_250, .after = ZONE )%>%
  mutate(across(
    contains("lsm_c_area_mn_LCG"), 
    ~ . / 0.0001))%>% ## Putting in m2 the surfaces which are in ha
  relocate(DATE_COLLECTE, .after=ID_COLLECTE)%>%
  relocate(NB_TOT, PROP_PP,.after=NB_ALBO_GRAVIDE)%>%
  relocate(MOIS, JOUR, .before=ZONE)%>%
  na.omit()
df$ZONE<-as.factor(df$ZONE)
df$AREA<-as.factor(df$AREA)

#########################################
# Function: Compare a list of GAM models
#########################################

compare_gams <- function(gam_list) {
  results <- data.frame(
    Model = character(),
    AIC = numeric(),
    BIC=numeric(),
    Explained_Variance = numeric(),
    R2 = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(gam_list)) {
    model <- gam_list[[i]]
    model_summary <- summary(model)
    
    results[i, ] <- list(
      Model = paste0("GAM_", i),
      AIC = AIC(model),
      BIC=BIC(model),
      Explained_Variance = round(model_summary$dev.expl * 100, 2),  # En pourcentage
      R2 = round(model_summary$r.sq, 4)
    )
  }
  
  return(results)
}


#######################################
# 1. Modelling the proportion of parous females
#######################################

### Null Model
M0<-gam(PROP_PP ~  1,
        data = df,
        family = binomial,
        weights = NB_TOT)

### --- A. By ZONE and MONTH ---

# Step 1: Compare spline types
gam_pp_cr <- gam(PROP_PP ~ ZONE + s(MOIS, k=5, bs='cr'), 
                 data = df,family = binomial, weights = NB_TOT,method='REML')
gam_pp_cs <- gam(PROP_PP ~ ZONE + s(MOIS, k=5, bs='cs') ,
                 data = df,
                 family = binomial,weights = NB_TOT, method='REML')
gam_pp_tp <- gam(PROP_PP ~ ZONE + s(MOIS, k=5, bs='tp') ,
                 data = df,
                 family = binomial,weights = NB_TOT, method='REML')
comp_gam<-compare_gams(list(gam_pp_cr,gam_pp_cs,gam_pp_tp,  M0))  # 'tp' chosen

# Step 2: Interaction test
gam_pp_tp_int <- gam(PROP_PP ~ ZONE + s(MOIS, k=5, bs='tp')+s(MOIS, by=ZONE, k=5), 
                 data = df,family = binomial, weights = NB_TOT,method='REML')
comp_gam<-compare_gams(list(gam_pp_tp_int,gam_pp_tp,  M0)) # no interaction retained

# Step 3: Random effect
gam_pp_al<-gam(PROP_PP ~ ZONE +   s(MOIS, k=5, bs='tp') + s(num_session, bs="re"),
               data = df,
               family = binomial,
               weights = NB_TOT,  select = TRUE, method='REML')
comp_gam<-compare_gams(list(gam_pp_al,gam_pp_tp,gam_pp_tp_int, M0)) # no random effect retained

# Final Model Summary and Diagnostics
sum_pp<-summary(gam_pp_tp)
#R-sq.(adj) =  0.144   Deviance explained = 16.6%
plot(gam_pp_tp, all.terms = TRUE, pages = 1, shift = coef(gam_pp_tp)[1])
vis.gam(gam_pp_tp, theta = 120, color = "heat")
par(mfrow = c(2,2))
gam.check(gam_pp_tp) ## ok

# Observed vs Predicted
df$pred <- predict(gam_pp_tp, type = "response")
ggplot(df, aes(x = pred, y = PROP_PP)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # y = x
  labs(
    title = "Observation vs Prediction (GAM)",
    x = "Prediction",
    y = "Observation"
  ) +
  annotate("text", x = 0.1, y = 0.95, label = paste0("R² = ", round(sum_pp$r.sq, 3)), size = 5, hjust = 0) +
  theme_minimal()

# Monthly Effect Prediction (by ZONE)
mois_grid <- expand.grid(# Creation of a dataframe
  MOIS = seq(5, 11, length.out = 300),
  ZONE = unique(df$ZONE),
  NB_TOT = 1
)
pred <- predict(gam_pp_tp, newdata = mois_grid, type = "link", se.fit = TRUE)
mois_grid$fit_link <- pred$fit
mois_grid$se_link <- pred$se.fit
mois_grid$fit <- gam_pp_tp$family$linkinv(mois_grid$fit_link)
mois_grid$lower <- gam_pp_tp$family$linkinv(mois_grid$fit_link - 2 * mois_grid$se_link)
mois_grid$upper <- gam_pp_tp$family$linkinv(mois_grid$fit_link + 2 * mois_grid$se_link)
ggplot() +
  geom_point(data = df, aes(x = MOIS, y = PROP_PP), alpha = 0.5) +
  geom_ribbon(data = mois_grid, aes(x = MOIS, ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
  geom_line(data = mois_grid, aes(x = MOIS, y = fit), color = "blue", size = 1.2) +
  labs(
    title = "Effet du mois (GAM) avec IC sur l'échelle de réponse",
    y = "Proportion prédite / observée"
  ) +
  theme_minimal()+
  facet_wrap(.~ZONE)


### --- B. By AREA and MONTH ---

# Spline type selection
gam_pp_tp_area <- gam(PROP_PP ~ AREA + s(MOIS, k=5, bs='tp'), data = df, family = binomial, weights = NB_TOT, method='REML')
gam_pp_cr_area <- gam(PROP_PP ~ AREA + s(MOIS, k=5, bs='cr'), data = df, family = binomial, weights = NB_TOT, method='REML')
gam_pp_cs_area <- gam(PROP_PP ~ AREA + s(MOIS, k=5, bs='cs'), data = df, family = binomial, weights = NB_TOT, method='REML')
gam_pp_tp_area_al <- gam(PROP_PP ~ AREA + s(MOIS, k=5, bs='tp') + s(num_session, bs="re"), data = df, family = binomial, weights = NB_TOT, method='REML')

compare_gams(list(gam_pp_cr_area, gam_pp_tp_area, gam_pp_cs_area, M0))

# Model summary and diagnostics
summary(gam_pp_tp_area)
plot(gam_pp_tp_area, all.terms = TRUE, pages = 1, shift = coef(gam_pp_tp)[1])
vis.gam(gam_pp_tp_area, theta = 120, color = "heat")
par(mfrow = c(2, 2)); gam.check(gam_pp_tp_area)

# Observed vs Predicted
df$pred <- predict(gam_pp_tp_area, type = "response")
ggplot(df, aes(x = pred, y = PROP_PP)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Observed vs Predicted (AREA Model)", x = "Prediction", y = "Observation") +
  annotate("text", x = 0.1, y = 0.95, label = paste0("R² = ", round(summary(gam_pp_tp_area)$r.sq, 3)), size = 5, hjust = 0) +
  theme_minimal()

# Prediction with confidence interval by AREA
mois_grid <- expand.grid(MOIS = seq(5, 11, length.out = 300), AREA = unique(df$AREA), NB_TOT = 1)
pred <- predict(gam_pp_tp_area, newdata = mois_grid, type = "link", se.fit = TRUE)
mois_grid$fit <- gam_pp_tp_area$family$linkinv(pred$fit)
mois_grid$lower <- gam_pp_tp_area$family$linkinv(pred$fit - 2 * pred$se.fit)
mois_grid$upper <- gam_pp_tp_area$family$linkinv(pred$fit + 2 * pred$se.fit)

ggplot() +
  geom_point(data = df, aes(x = MOIS, y = PROP_PP), alpha = 0.5) +
  geom_ribbon(data = mois_grid, aes(x = MOIS, ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
  geom_line(data = mois_grid, aes(x = MOIS, y = fit), color = "blue", size = 1.2) +
  labs(title = "Monthly Effect (with CI)", y = "Predicted / Observed Parturity") +
  facet_wrap(. ~ AREA) +
  theme_minimal()

plot_pred_part<-ggplot() +
  # Points gris : observations
  geom_point(data = df, aes(x = MOIS, y = PROP_PP, color = "Observations"), alpha = 0.5) +
  
  # Bande bleue : intervalle de confiance
  geom_ribbon(data = mois_grid, aes(x = MOIS, ymin = lower, ymax = upper, fill = "95% CI"), alpha = 0.2) +
  
  # Ligne bleue : prédictions
  geom_line(data = mois_grid, aes(x = MOIS, y = fit, color = "Prediction"), size = 1.2) +
  
  labs(
    title = "GAM-Predicted Parturity Rate by Area and Month",
    y = "Parturity Rate",
    x = "Month",
    color = "Legend",
    fill = "Legend"
  ) +
  
  scale_color_manual(values = c("Observations" = "grey50", "Prediction" = "blue")) +
  scale_fill_manual(values = c("95% CI" = "blue")) +
  
  facet_wrap(. ~ AREA) +
  theme_minimal() +
  theme(legend.position = "top")+
  theme(
    margin(t = 1, r = 2, b = 1, l = 1, unit = "cm"),
    legend.position = "top",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle=45, hjust = 1),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )

ggsave("article/prop_area_gam.svg", plot_pred_part,
       width = 18, height = 17, units = "cm", dpi = 300, limitsize = F)




### --- B. By AREA and MONTH ---

# Spline type selection
gam_pp_tp_area <- gam(PROP_PP ~ AREA + s(MOIS, k=5, bs='tp'), data = df, family = binomial, weights = NB_TOT, method='REML')
gam_pp_cr_area <- gam(PROP_PP ~ AREA + s(MOIS, k=5, bs='cr'), data = df, family = binomial, weights = NB_TOT, method='REML')
gam_pp_cs_area <- gam(PROP_PP ~ AREA + s(MOIS, k=5, bs='cs'), data = df, family = binomial, weights = NB_TOT, method='REML')
gam_pp_tp_area_al <- gam(PROP_PP ~ AREA + s(MOIS, k=5, bs='tp') + s(num_session, bs="re"), data = df, family = binomial, weights = NB_TOT, method='REML')

compare_gams(list(gam_pp_tp_area, M0))

# Model summary and diagnostics
summary(gam_pp_tp_area)
plot(gam_pp_tp_area, all.terms = TRUE, pages = 1, shift = coef(gam_pp_tp)[1])
vis.gam(gam_pp_tp_area, theta = 120, color = "heat")
par(mfrow = c(2, 2)); gam.check(gam_pp_tp_area)

# Observed vs Predicted
df$pred <- predict(gam_pp_tp_area, type = "response")
ggplot(df, aes(x = pred, y = PROP_PP)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Observed vs Predicted (AREA Model)", x = "Prediction", y = "Observation") +
  annotate("text", x = 0.1, y = 0.95, label = paste0("R² = ", round(summary(gam_pp_tp_area)$r.sq, 3)), size = 5, hjust = 0) +
  theme_minimal()

# Prediction with confidence interval by AREA
mois_grid <- expand.grid(MOIS = seq(5, 11, length.out = 300), AREA = unique(df$AREA), NB_TOT = 1)
pred <- predict(gam_pp_tp_area, newdata = mois_grid, type = "link", se.fit = TRUE)
mois_grid$fit <- gam_pp_tp_area$family$linkinv(pred$fit)
mois_grid$lower <- gam_pp_tp_area$family$linkinv(pred$fit - 2 * pred$se.fit)
mois_grid$upper <- gam_pp_tp_area$family$linkinv(pred$fit + 2 * pred$se.fit)
pred_prop<-ggplot() +
   geom_point(data = df, aes(x = MOIS, y = PROP_PP), alpha = 0.5) +
     geom_ribbon(data = mois_grid, aes(x = MOIS, ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
     geom_line(data = mois_grid, aes(x = MOIS, y = fit), color = "blue", size = 1.2) +
     labs(title = "Monthly Effect (with CI)", y = "Predicted / Observed Proportion") +
    facet_wrap(. ~ AREA) +
     theme_minimal()
prop_pp<-ggplot() +
  geom_point(data = df, aes(x = MOIS, y = PROP_PP), alpha = 0.5) +
  geom_smooth(data = df, aes(x = MOIS, y = PROP_PP), color = "blue", size = 1.2,se = FALSE) +
  labs(title = "Monthly Effect (with CI)", y = "Observed Proportion") +
  ylim(c(0,1))+
  facet_wrap(. ~ AREA) +
  theme_minimal()

prop_area <- arrangeGrob(pred_prop, prop_pp, ncol = 2)
ggsave("plots/prop_area.pdf", prop_area,
       width = 16.5, height = 11.7, units = "in")

################################################################################
# 2. Modelling the number of females caught (NB_ALBO_F_DIS)
################################################################################

### Null Model
M0_dis <- gam(NB_ALBO_F_DIS ~ 1, data = df, family = nb())

### --- A. By ZONE and MONTH ---

# Step 1: Compare spline types
gam_dis_cr <- gam(NB_ALBO_F_DIS ~ ZONE + s(MOIS, k=5, bs='cr'), data = df, family = nb(), method='REML')
gam_dis_cs <- gam(NB_ALBO_F_DIS ~ ZONE + s(MOIS, k=5, bs='cs'), data = df, family = nb(), method='REML')
gam_dis_tp <- gam(NB_ALBO_F_DIS ~ ZONE + s(MOIS, k=5, bs='tp'), data = df, family = nb(), method='REML')

compare_gams(list(gam_dis_cr, gam_dis_cs, gam_dis_tp, M0_dis)) # Select 'tp'

# Step 2: Interaction test
gam_dis_tp_int <- gam(NB_ALBO_F_DIS ~ ZONE + s(MOIS, k=5, bs='tp') + s(MOIS, by=ZONE, k=5), data = df, family = nb, method='REML')
compare_gams(list(gam_dis_tp_int, gam_dis_tp, M0_dis)) # Keep the one with interaction if improved

# Step 3: Random effect (session)
gam_dis_tp_re <- gam(NB_ALBO_F_DIS ~ ZONE + s(MOIS, k=5, bs='tp') + s(num_session, bs='re'), data = df, family = nb, method='REML')
compare_gams(list(gam_dis_tp_re, gam_dis_tp_int, gam_dis_tp, M0_dis)) # Keep best

# Final Model Summary & Diagnostics
summary(gam_dis_tp_int)
plot(gam_dis_tp_int, all.terms = TRUE, pages = 1)
vis.gam(gam_dis_tp_int, theta = 120, color = "heat")
par(mfrow = c(2, 2)); gam.check(gam_dis_tp_int)

# Observed vs Predicted
df$pred_dis <- predict(gam_dis_tp_int, type = "response")
ggplot(df, aes(x = pred_dis, y = NB_ALBO_F_DIS)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Observed vs Predicted Catches", x = "Prediction", y = "Observation") +
  annotate("text", x = max(df$pred_dis)*0.1, y = max(df$NB_ALBO_F_DIS)*0.9,
           label = paste0("R² = ", round(summary(gam_dis_tp_int)$r.sq, 3)), size = 5) +
  theme_minimal()

# Monthly prediction by ZONE
mois_grid <- expand.grid(MOIS = seq(5, 11, length.out = 300), ZONE = unique(df$ZONE))
pred <- predict(gam_dis_tp_int, newdata = mois_grid, type = "link", se.fit = TRUE)
mois_grid$fit <- exp(pred$fit)
mois_grid$lower <- exp(pred$fit - 2 * pred$se.fit)
mois_grid$upper <- exp(pred$fit + 2 * pred$se.fit)

ggplot() +
  geom_point(data = df, aes(x = MOIS, y = NB_ALBO_F_DIS), alpha = 0.3) +
  geom_ribbon(data = mois_grid, aes(x = MOIS, ymin = lower, ymax = upper), fill = "orange", alpha = 0.2) +
  geom_line(data = mois_grid, aes(x = MOIS, y = fit), color = "orange", size = 1.2) +
  labs(title = "Monthly Effect (GAM Prediction)", y = "Predicted / Observed Females Caught") +
  facet_wrap(. ~ ZONE) +
  theme_minimal()

### --- B. By AREA and MONTH ---

# Step 1: Compare spline types
gam_dis_tp_area <- gam(NB_ALBO_F_DIS ~ AREA + s(MOIS, k=5, bs='tp'), data = df, family = nb, method='REML')
gam_dis_cr_area <- gam(NB_ALBO_F_DIS ~ AREA + s(MOIS, k=5, bs='cr'), data = df, family = nb, method='REML')
gam_dis_cs_area <- gam(NB_ALBO_F_DIS ~ AREA + s(MOIS, k=5, bs='cs'), data = df, family = nb, method='REML')
gam_dis_tp_area_re <- gam(NB_ALBO_F_DIS ~ AREA + s(MOIS, k=5, bs='tp') + s(num_session, bs="re"), data = df, family = nb, method='REML')

compare_gams(list(gam_dis_tp_area_re,  M0_dis))

# Final model summary & diagnostics
summary(gam_dis_tp_area)
plot(gam_dis_tp_area, all.terms = TRUE, pages = 1)
vis.gam(gam_dis_tp_area, theta = 120, color = "heat")
par(mfrow = c(2, 2)); gam.check(gam_dis_tp_area)

# Observed vs Predicted
df$pred_dis_area <- predict(gam_dis_tp_area, type = "response")
ggplot(df, aes(x = pred_dis_area, y = NB_ALBO_F_DIS)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Observed vs Predicted (AREA Model)", x = "Prediction", y = "Observation") +
  annotate("text", x = max(df$pred_dis_area)*0.1, y = max(df$NB_ALBO_F_DIS)*0.9,
           label = paste0("R² = ", round(summary(gam_dis_tp_area)$r.sq, 3)), size = 5) +
  theme_minimal()

# Monthly prediction by AREA
mois_grid <- expand.grid(MOIS = seq(5, 11, length.out = 300), AREA = unique(df$AREA))
pred <- predict(gam_dis_tp_area, newdata = mois_grid, type = "link", se.fit = TRUE)
mois_grid$fit <- exp(pred$fit)
mois_grid$lower <- exp(pred$fit - 2 * pred$se.fit)
mois_grid$upper <- exp(pred$fit + 2 * pred$se.fit)

ma_pred<-ggplot() +
  geom_point(data = df, aes(x = MOIS, y = NB_ALBO_F_DIS), alpha = 0.3) +
  geom_ribbon(data = mois_grid, aes(x = MOIS, ymin = lower, ymax = upper), fill = "grey", alpha = 0.2) +
  geom_line(data = mois_grid, aes(x = MOIS, y = fit), color = "blue", size = 1.2) +
  labs(title = "Monthly Effect (by AREA)", y = "Predicted / Observed Females Caught") +
  facet_wrap(. ~ AREA) +
  theme_minimal()

ma_obs<-ggplot() +
  geom_smooth(data = df, aes(x = MOIS, y = NB_ALBO_F_DIS), fill = "blue", size = 1.2, se=F) +
  labs(title = "Monthly Effect (by AREA)", y = "Observed Females Caught") +
  facet_wrap(. ~ AREA) +
  theme_minimal()

ma_area <- arrangeGrob(ma_pred, ma_obs, ncol = 2)
ggsave("plots/ma_area.pdf", ma_area,
       width = 16.5, height = 11.7, units = "in")


plot_pred_ma<-ggplot() +
  # Points gris : observations
  geom_point(data = df, aes(x = MOIS, y = NB_ALBO_F_DIS, color = "Observations"), alpha = 0.5) +
  
  # Bande bleue : intervalle de confiance
  geom_ribbon(data = mois_grid, aes(x = MOIS, ymin = lower, ymax = upper, fill = "95% CI"), alpha = 0.2) +
  
  # Ligne bleue : prédictions
  geom_line(data = mois_grid, aes(x = MOIS, y = fit, color = "Prediction"), size = 1.2) +
  
  labs(
    title = "GAM-Predicted Female abundance by Area and Month",
    y = "Female abundance",
    x = "Month",
    color = "Legend",
    fill = "Legend"
  ) +
  
  scale_color_manual(values = c("Observations" = "grey50", "Prediction" = "blue")) +
  scale_fill_manual(values = c("95% CI" = "blue")) +
  
  facet_wrap(. ~ AREA) +
  theme_minimal() +
  theme(legend.position = "top")+
  theme(
    margin(t = 1, r = 2, b = 1, l = 1, unit = "cm"),
    legend.position = "top",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle=45, hjust = 1),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )

ggsave("article/ma_area_gam.svg", plot_pred_ma,
       width = 18, height = 17, units = "cm", dpi = 300, limitsize = F)


gam_area <- arrangeGrob(plot_pred_part, plot_pred_ma, ncol = 2)
ggsave("article/gam_area.svg", plot_pred_ma,
       width = 18, height = 17, units = "cm", dpi = 300, limitsize = F)
