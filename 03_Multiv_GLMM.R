###############################'To build multivariate GLMM to explain parity rate

#######'Packages
rm(list=ls())
library(ggplot2) #  ‘3.5.1’
library(correlation) # 0.8.6’
library(glmmTMB)  # ‘1.9.9’
library(performance)# ‘0.13.0’
library(DHARMa) # '0.4.6 '
library(MuMIn) #  ‘1.48.4’
library(spdep) # ‘1.3.5’
library(sf) # ‘1.0.16’
library(tidyverse) # ‘2.0.0’
library(dplyr) # ‘1.1.2’"


#######'Opening data

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

df<-cbind(df, model.matrix(~ZONE - 1, data = df)) # Allow to transform a categorial variable in a boolean variable 

#######'Graphic representation

predictors<-c("NB_ALBO_F",
              "RHMAX_collection", 
              "BUFFER_20",
              "lsm_c_np_LCG_20_13", "lsm_c_te_LCG_100_12", "lsm_c_pland_LCG_20_12", "lsm_c_area_mn_LCG_250_12",
              "lsm_c_te_LCG_100_13", "lsm_c_pland_LCG_100_13","lsm_c_area_mn_LCG_100_13",
              "lsm_c_te_LCG_20_11", "lsm_c_pland_LCG_20_11","lsm_c_area_mn_LCG_20_11",
              "lsm_c_te_LCG_20_10", "lsm_c_pland_LCG_20_10","lsm_c_area_mn_LCG_20_10",
              "POP_50_sum", "ZONEResidential")



# Plot the bivariate relationship between proportion of parous and each selected predictor 

p_prop <- df %>% 
  dplyr::select(PROP_PP,predictors) %>%
  pivot_longer(-PROP_PP) %>%
 ggplot(aes(y = PROP_PP, x = value)) +
  geom_point() + 
  ylim(c(0,1)) +
  geom_smooth(method = "lm", formula = y ~ x) +
  # geom_smooth(method='lm') +
  facet_wrap(.~name, scales = "free") + 
  theme_bw() + 
  ggtitle("Prop parous ~Variables")

#######'Correlation between variables
predictors<-c("NB_ALBO_F",
              "RHMAX_collection", 
              "BUFFER_20",
              "lsm_c_np_LCG_20_13", "lsm_c_te_LCG_100_12", "lsm_c_pland_LCG_20_12", "lsm_c_area_mn_LCG_250_12",
              "lsm_c_te_LCG_100_13", "lsm_c_pland_LCG_100_13","lsm_c_area_mn_LCG_100_13",
              "lsm_c_te_LCG_20_11", "lsm_c_pland_LCG_20_11","lsm_c_area_mn_LCG_20_11",
              "lsm_c_te_LCG_20_10", "lsm_c_pland_LCG_20_10","lsm_c_area_mn_LCG_20_10",
              "POP_50_sum",  "ZONEResidential")
m <- cor(df[,predictors], method = "pearson", use = "pairwise.complete.obs")
index <- which(abs(m) > .7 & abs(m) < 1,arr.ind = T) 
df_cor <- subset(as.data.frame(index) , row <= col)
p <- cbind.data.frame(stock1 = rownames(m)[df_cor[,1]], stock2 = colnames(m)[df_cor[,2]])


#######'VIF

# Function for calculating VIFs for a model based on a specific target
calculate_vif<- function(data, target_column) {
  # Extract names predictors
  predictors <- setdiff(colnames(data), target_column)
  
  # Calculate VIF for each predictor
  vif_values <- sapply(predictors, function(var) {
    formula <- as.formula(paste(target_column, "~ ."))
    lm_fit <- lm(formula, data = data[, c(target_column, predictors), drop = FALSE])
    r_squared <- summary(lm(lm_fit$model[[var]] ~ ., data = lm_fit$model[, -which(names(lm_fit$model) == var)]))$r.squared
    vif_value <- 1 / (1 - r_squared)
    return(vif_value)
  })
  
  return(vif_values)
}

# Function to reduce VIF-based variables for a given target
reduce_vif <- function(data, target_column, vif_threshold = 10) {
  working_data <- data
  predictors <- setdiff(colnames(data), target_column)
  
  while (TRUE) {
    # Calculate VIF
    vif_values <- calculate_vif(working_data, target_column)
    cat("Current VIF values:\n")
    print(vif_values)
    
    # Check that all VIFs are below the threshold
    if (all(vif_values < vif_threshold)) break
    
    # Identify the variable with the highest VIF
    variable_to_remove <- names(which.max(vif_values))
    cat("Removing variable with high VIF:", variable_to_remove, "\n")
    
    # Delete variables
    predictors <- setdiff(predictors, variable_to_remove)
    working_data <- working_data[, c(target_column, predictors), drop = FALSE]
  }
  
  return(working_data)
}

var_selected<-reduce_vif(df[,c(predictors, "PROP_PP")],target_column="PROP_PP", vif_threshold = 3)
c<-colnames(var_selected)

############################################################'GLMM final 

df_glmm<-df%>%
  dplyr::select(PROP_PP, NB_TOT,ID_PIEGE,AREA, num_session, MOIS, NB_ALBO_F,ZONE,LATITUDE, LONGITUDE, predictors) %>% 
  mutate( PROP_PP = as.character(PROP_PP), NB_TOT=as.character(NB_TOT),  NB_ALBO_F=as.character(NB_ALBO_F), MOIS=as.character(MOIS),LATITUDE = as.character(LATITUDE),LONGITUDE = as.character(LONGITUDE)) %>%
  mutate_if(is.numeric, ~scale(., center = TRUE, scale = FALSE)) %>%  # we center  to be able to compare the magnitudes (centering also helps with allowing easier interpretation of variable coefficients associated with different magnitudes, e.g. when one regressor is measured on a very small order, while another on a very large order.  )
  mutate(PROP_PP = as.numeric(PROP_PP), NB_TOT=as.integer(NB_TOT), NB_ALBO_F=as.integer(NB_ALBO_F), MOIS=as.integer(MOIS),
         LATITUDE = as.numeric(LATITUDE),LONGITUDE = as.numeric(LONGITUDE))%>%
  arrange(MOIS)%>%
  na.omit()

####'Full model witjout variables selection
null <- glmmTMB(PROP_PP ~ 1 ,
                                      data = df_glmm, family = binomial(link = "logit"), weights = NB_TOT) #" nul model
null_model <- glmmTMB(PROP_PP ~ 1+(1|num_session) ,
                data = df_glmm, family = binomial(link = "logit"), weights = NB_TOT) ## nul model with random effects

mod_fin<-glmmTMB(PROP_PP ~  RHMAX_collection +   BUFFER_20+ lsm_c_te_LCG_20_10+lsm_c_area_mn_LCG_250_12+lsm_c_area_mn_LCG_100_13+lsm_c_area_mn_LCG_20_10+POP_50_sum+
                                       NB_ALBO_F +  (1|num_session),
                                           data = df_glmm, family = binomial(link = "logit"),  weights = NB_TOT) ## full mdoel 

AIC(null, null_model, mod_fin)
#  df      AIC
#null        1 357.7970
#null_model  2 354.9140
#mod_fin    10 343.7697


####'Automated selection based on aic via dredge
options(na.action = "na.fail")  
dredged <- dredge(mod_fin, trace = TRUE, rank = "AIC")
best_model <- get.models(dredged, 1)[[1]]

mod_fin_select<-glmmTMB( PROP_PP ~ lsm_c_area_mn_LCG_100_13 + lsm_c_area_mn_LCG_250_12 +      lsm_c_te_LCG_20_10 + RHMAX_collection + (1 | num_session),
                       data = df_glmm,weights = NB_TOT, family = binomial(link = "logit"))
                        
summary(mod_fin_select)
table_mod<-broom.mixed::tidy(mod_fin_select, conf.int = TRUE,  exponentiate = TRUE)
performance::r2_nakagawa(mod_fin_select)
DHARMa::simulateResiduals(mod_fin_select, plot=T) ## OK
summary(mod_select_fin)
AIC(null,null_model, mod_fin, mod_fin_select)
r.squaredGLMM(mod_fin_select)
#df      AIC
#null            1 357.7970
#null_model      2 354.9140
#mod_fin        10 343.7697
#mod_fin_select  6 335.9415
performance::check_model(mod_fin_select)

####'Spatial auto correlation 
## Distance matrice
coords <- st_as_sf(df_glmm, coords = c("LONGITUDE", "LATITUDE"), crs = 4326)
k_neigh <- knearneigh(st_coordinates(coords), k = 5, longlat = TRUE)
nb <- knn2nb(k_neigh)
# Weight matrices
listw <- nb2listw(nb, style = "W", zero.policy = TRUE)
df_glmm$residus <- residuals(mod_fin_select, type = "pearson")
moran_result <- moran.test(df_glmm$residus, listw, zero.policy = TRUE)
print(moran_result) ## no spatial auto correlation 

#######'Temporal auto correlation
res<-DHARMa::simulateResiduals(mod_fin_select)
res_agg <- recalculateResiduals(res, group = df$MOIS)
testTemporalAutocorrelation(res_agg, time = unique(df$MOIS)) ## no temporal autocorrelation

#######'% variance explained
vars<-insight::get_variance(mod_fin_select)
total_var <- vars$var.fixed + vars$var.random + vars$var.residual

pct_fixed  <- vars$var.fixed   / total_var * 100
pct_random <- vars$var.random  / total_var * 100
pct_resid  <- vars$var.residual / total_var * 100

cat("Fixed effects :", round(pct_fixed, 1), "%\n")
cat("Random effects :", round(pct_random, 1), "%\n")
cat("Residual :", round(pct_resid, 1), "%\n")


######################'Cross validation on AREA


cv_col <- "AREA"
indices_cv <- CAST::CreateSpacetimeFolds(df_glmm, spacevar = cv_col, 
                                         k = length(unique(unlist(df_glmm[, cv_col]))))
df_glmm$predicted <- NA

for (i in seq_along(indices_cv$index)) {
  

  train_indices <- indices_cv$index[[i]]
  test_indices <- indices_cv$indexOut[[i]]
  
  train_data <- df_glmm[train_indices, ]
  test_data <- df_glmm[test_indices, ]
  
  glmm_model<-glmmTMB( PROP_PP ~ lsm_c_area_mn_LCG_100_13 + lsm_c_area_mn_LCG_250_12 +      lsm_c_te_LCG_20_10 + RHMAX_collection + (1 | num_session),
                           data = df_glmm,weights = NB_TOT, family = binomial(link = "logit"))
  test_data$predicted <- predict(glmm_model, newdata = test_data, type = "response")
  df_glmm$predicted[test_indices] <- test_data$predicted
}


df_glmm$num_session<-as.numeric(df_glmm$num_session)

## Model evaluation plots 
plot_eval_parity_model <- df_glmm %>%
  dplyr::group_by(AREA,num_session) %>%   
  dplyr::summarise(pred = mean(predicted), obs = mean(PROP_PP)) %>%
  as_tibble() %>%
  pivot_longer(c('pred','obs')) %>%
  dplyr::mutate(name = ifelse(name=="pred","Predicted","Observed")) %>%
  ggplot(aes(x=num_session, y = value, color = name)) +
  geom_point() + 
  geom_line() + 
  facet_wrap(~AREA, scales = "free") + 
  theme_bw() + 
  scale_colour_manual(values=c("#009E73","#E69F00"),na.translate = F) + 
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9)) +
  ylim(c(0,1))+
  xlab("Sampling session") +
  ylab("Parturity rate") + 
  labs(color='Parturity rate of Ae. albopictus') + 
  theme(legend.position="bottom") + 
  ggtitle('Parturity rate model: observed vs. predicted values')+theme(
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


ggsave(filename = "article/prop_parous_plot_prediction_vs_observtion.svg",plot =plot_eval_parity_model, width = 18, height = 17, units = "cm", dpi = 300, limitsize = F)
