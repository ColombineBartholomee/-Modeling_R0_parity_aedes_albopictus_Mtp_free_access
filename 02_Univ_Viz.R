###############################'Graphic representation of univariate modelling 

rm(list=ls())
library(tidyverse) # ‘2.0.0’
library(purrr) # ‘1.0.1’
library(patchwork) # ‘1.3.0’
library(landscapemetrics) # ‘2.1.4’

source("analysis/functions.R") ## Function to represent graphics of univariate modelling

########################
#### GLMM univariate pour abondance totale
########################

# open results of GLMM
#glmm_univ_parity <- read.csv("data/glmm_parity_traps.csv", stringsAsFactors = F)
glmm_univ_parity <- read.csv("data/glmm_prop_parous_traps.csv", stringsAsFactors = F)
glmm_univ_parity$indicator <- 'parity'

#### Micro-climatic data representation
univ_glmm_microclim <- glmm_univ_parity %>%
  filter(p.value<0.2)%>%
  filter(grepl("TMEAN|TMIN|TMAX|RHMEAN|RHMIN|RHMAX",term)) %>%
  filter(grepl("collection|prec",term)) %>%
  mutate(type = sub('\\_.*', '', term)) %>%
  mutate(type = gsub("RFSUM|RFMAX|RFMIN","Rainfall",type)) %>%
  mutate(type = gsub("TMEAN|TMIN|TMAX|TAMP","Temperature",type)) %>%
  mutate(type = gsub("RHMEAN|RHMIN|RHMAX","Humidity",type)) %>%
  mutate(type=gsub("Patmean|Patmin|Patmax|Patm_diff_prev_day|Patm_daily_range", "Patmos", type))%>%
  mutate(type=ifelse(type=="Patm", "Patmos", type))%>%
  mutate(buffer = case_when(grepl("collection",term) ~ "during coll.",
                            grepl("24h_48h",term) ~ "24h and 48h before coll.",
                            grepl("24h",term) ~ "24h",
                            grepl("48h",term) ~ "48h")) %>%
  mutate(term = gsub("collection","during collection",term)) %>%
  mutate(term = gsub("24hprec|24h_prec","24h preceding collection",term)) %>%
  mutate(term = gsub("24h_48h_","24h and 48h preceding collection",term)) %>%
  mutate(term = gsub("48hprec|48h_prec","48h preceding collection",term)) %>%
  mutate(label = case_when(grepl("MAX|max",term) ~ "Maximum",
                           grepl("MIN|min",term) ~ "Minimum",
                           grepl("MEAN|SUM|mean",term) ~ "Average")) %>%
  mutate(term = gsub("RFSUM","Cumulative/average",term)) %>%
  mutate(term = gsub("RFMAX","Maximum",term)) %>%
  mutate(term = gsub("RFMIN","Minimum",term)) %>%
  mutate(term = gsub("TMEAN","Cumulative/average",term)) %>%
  mutate(term = gsub("TMIN","Minimum",term)) %>%
  mutate(term = gsub("TMAX","Maximum",term)) %>%
  mutate(term = gsub("TAMP","Amplitude",term)) %>%
  mutate(term = gsub("RHMEAN","Cumulative/average",term)) %>%
  mutate(term = gsub("RHMIN","Minimum",term)) %>%
  mutate(term = gsub("RHMAX","Maximum",term)) %>%
  mutate(term = gsub("_"," ",term)) %>%
  mutate(correlation = ifelse(p.value<=0.2,estimate,NA)) %>%
  mutate(r2 = ifelse(p.value<=0.2,r2,NA)) %>%
  mutate(indicator = forcats::fct_relevel(indicator, c("presence","abundance"))) %>%
  mutate(label = forcats::fct_relevel(label, c("Minimum","Maximum","Average"))) %>%
  mutate(buffer = forcats::fct_relevel(buffer, c("during coll.","24h","48h")))


p_microclim <- fun_plot_tile_univ_spatial_r2(univ_glmm_microclim, metric_name = "glmm", indicator = univ_glmm_microclim$term, lc_source = "Microclimatic conditions", type = "", xlabel = "time before collection")

ggsave(filename = "data/plots/prop_parous/microclim_glmm.pdf",plot = p_microclim, device = "pdf", width = 11, height = 8)


#### Spatial data - landscape metrics
landcover_grouped_veget_data_dict <- read.csv("data/landcover_grouped_veget_data_dic.csv", stringsAsFactors = F) %>% mutate(lc_source = "LCG")
data_dic <- read.csv("data/data_dictionary.csv", stringsAsFactors = F) 

lsm_data_dic <- landcover_grouped_veget_data_dict %>%
  mutate(class = as.character(class))

univ_glmm_lsm <- glmm_univ_parity %>%
  filter(grepl("lsm|POP",term)) %>%
  dplyr::mutate(function_name = sub("\\_L.*", "", term), lc_source = word(term, -3, sep = "_"), buffer = word(term, -2, sep = "_"), class = word(term, -1, sep = "_")) %>%
  dplyr::mutate(function_name=ifelse(grepl("sd", function_name),"sd",
                                      ifelse(grepl("sum", function_name),"sum",function_name)))%>%
  dplyr:: mutate(buffer = forcats::fct_relevel(buffer, c("0","20","50","100","250"))) %>%
  left_join(list_lsm()) %>%
  left_join(lsm_data_dic) %>%
  rename(correlation = estimate) %>%
  dplyr::mutate(label = ifelse(level=="landscape","landscape",label)) %>%
  dplyr::mutate(label = paste0(label, " - ", metric," - ",name)) %>%
  #nest(-c(lc_source,type))
  nest(-lc_source)

plots_univ_glmm_spatial <- univ_glmm_lsm %>%
  #mutate(univ_spatial = pmap(list(data,lc_source,type), ~fun_plot_tile_univ_spatial(correlation_df = ..1, metric_name = "glmm", indicator = ..1$indicator, lc_source = ..2, type = ..3))) %>%
  mutate(univ_spatial = pmap(list(data,lc_source), ~fun_plot_tile_univ_spatial_r2(correlation_df = ..1, metric_name = "glmm", indicator = ..1$indicator, lc_source = ..2, type = "", xlabel = "buffer radius around the collection site (meters)"))) %>%
  dplyr::select(-data)
pmap(list(plots_univ_glmm_spatial$lc_source,plots_univ_glmm_spatial$univ_spatial),~ggsave(filename = paste0("data/plots/prop_parous/",..1,".pdf"),plot = ..2, device = "pdf", width = 7,height = 8))

####Breeding sites

fil <- glmm_univ_parity %>%
  filter(grepl("BUFFER",term)) %>%
  mutate(r2_modif=dplyr::case_when(estimate >=1~r2, 
                                   estimate<1~r2*(-1)),buffer = as.factor(word(term, -1, sep = "_"))) %>%
  filter(p.value<0.25)%>%
  mutate(p.value2 = case_when(
    p.value <= 0.001 ~ "***",
    p.value > 0.001 & p.value <= 0.01  ~  "**",
    p.value > 0.01 & p.value <= 0.05 ~ " *",
    p.value > 0.05 ~ ""
  ))


p <- ggplot(fil, aes(buffer, term)) + 
  geom_tile(aes(fill = r2_modif), color = "white") + 
  facet_grid(~indicator, scales="free_y", space="free_y") +
  #facet_grid(.~indicator, scales="free_y", space="free_y") +
  xlab("") + 
  ylab("") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10, face = "italic"),
        axis.text.y = element_text(size = 8)
  ) +
  #geom_text(aes(label = ifelse(is.na(correlation), "",paste(round(correlation,2), p.value2))), size = 3)
  geom_text(aes(label = ifelse(is.na(estimate), "", p.value2)), size = 3) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", limit = c(-0.1,0.10),midpoint = 0, space = "Lab", name = "r2 ajuste", na.value="grey")


ggsave(filename = "data/plots/prop_parous/breeding_sites.pdf",plot = p, device = "pdf")


#### Sampling environments


fil <- glmm_univ_parity %>%
  filter(grepl("ZONE",term)) %>%
  mutate(r2_modif=dplyr::case_when(estimate >=1~r2, 
                                   estimate<1~r2*(-1)),buffer = as.factor(word(term, -1, sep = "_"))) %>%
  filter(p.value<0.2)%>%
  mutate(p.value2 = case_when(
    p.value <= 0.001 ~ "***",
    p.value > 0.001 & p.value <= 0.01  ~  "**",
    p.value > 0.01 & p.value <= 0.05 ~ " *",
    p.value > 0.05 ~ ""
  ))

p <- ggplot(fil, aes(buffer, term)) + 
  geom_tile(aes(fill = r2_modif), color = "white") + 
  facet_grid(~indicator, scales="free_y", space="free_y") +
  #facet_grid(.~indicator, scales="free_y", space="free_y") +
  xlab("") + 
  ylab("") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10, face = "italic"),
        axis.text.y = element_text(size = 8)
  ) +
  #geom_text(aes(label = ifelse(is.na(correlation), "",paste(round(correlation,2), p.value2))), size = 3)
  geom_text(aes(label = ifelse(is.na(estimate), "", p.value2)), size = 3) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", limit = c(-0.1,0.10),midpoint = 0, space = "Lab", name = "r2 ajuste", na.value="grey")


ggsave(filename = "data/plots/prop_parous/zone.pdf",plot = p, device = "pdf")

