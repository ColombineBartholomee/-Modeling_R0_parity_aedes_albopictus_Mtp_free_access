#############################################
# Title: Plotting R0 Estimates by Exposure, Zone, and Month
# Description: This script visualizes theoretical and simulated R0 values
#              (including confidence intervals) based on environmental 
#              and human exposure scenarios.
#############################################

### Load required packages
rm(list = ls())  # Clear environment

library(ggplot2) #  ‘3.5.1’
library(tidyverse) # ‘2.0.0’
library(dplyr) # ‘1.1.2’"
library(patchwork)# ‘1.3.0’
library(gridExtra) #‘2.3’

### Load dataset
df_fin <- read.csv("data/df_fin.csv")

#############################################
# SECTION 1: Raw Boxplot Visualizations
#############################################

# 1.1 Theoretical R0 (basic)
p_basic <- df_fin %>%
  ggplot(aes(x = as.factor(MOIS), y = R0_basic, fill = ZONE)) +
  geom_boxplot() +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  facet_grid(expo ~ ZONE, scales = "free_y", 
             labeller = labeller(expo = function(x) paste("Expo:", x))) +
  xlab("Month") +
  ggtitle("Theoretical R0 Basic by Month and Environment")

# 1.2 Theoretical R0 (Poletti method)
p_poletti <- df_fin %>%
  ggplot(aes(x = as.factor(MOIS), y = R0_poletti_r0, fill = ZONE)) +
  geom_boxplot() +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  facet_grid(expo ~ ZONE, scales = "free_y",
             labeller = labeller(expo = function(x) paste("Expo:", x))) +
  xlab("Month") +
  ggtitle("Theoretical R0 (Poletti) by Month and Environment")


# 1.3 Simulated R0 (basic)
p_basic_sim_mean <- df_fin %>%
  ggplot(aes(x = as.factor(MOIS), y = R0_basic_sim_mean, fill = ZONE)) +
  geom_boxplot() +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  facet_grid(expo ~ ZONE, scales = "free_y",
             labeller = labeller(expo = function(x) paste("Expo:", x))) +
  xlab("Month") +
  ggtitle("Simulated Average R0 Basic by Month and Environment")

# 1.4 Simulated R0 (Poletti method)
p_poletti_sim_mean <- df_fin %>%
  ggplot(aes(x = as.factor(MOIS), y = R0_poletti_sim_mean, fill = ZONE)) +
  geom_boxplot() +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  facet_grid(expo ~ ZONE, scales = "free_y",
             labeller = labeller(expo = function(x) paste("Expo:", x))) +
  xlab("Month") +
  ggtitle("Simulated Average R0 (Poletti) by Month and Environment")

#############################################
# SECTION 2: Summarized Data with Means and Confidence Intervals by zone and month and exposure
#############################################

# Grouping by zone, month, and exposure to compute means
df_zone <- df_fin %>%
  group_by(ZONE, MOIS, expo) %>%
  summarise(
    ma = mean(ma, na.rm = TRUE),
    parity_rate = mean(parity_rate, na.rm = TRUE),
    R0_basic = mean(R0_basic, na.rm = TRUE),
    p = mean(p, na.rm = TRUE),
    R0_basic_sim_mean = mean(R0_basic_sim_mean, na.rm = TRUE),
    R0_basic_low_quantile = mean(R0_basic_low_quantile, na.rm = TRUE),
    R0_basic_high_quantile = mean(R0_basic_high_quantile, na.rm = TRUE),
    R0_poletti_r0 = mean(R0_poletti_r0, na.rm = TRUE),
    R0_poletti_sim_mean = mean(R0_poletti_sim_mean, na.rm = TRUE),
    R0_poletti_low_quantile = mean(R0_poletti_low_quantile, na.rm = TRUE),
    R0_poletti_high_quantile = mean(R0_poletti_high_quantile, na.rm = TRUE)) %>%
  mutate(MOIS = as.numeric(as.character(MOIS)))  # Ensure numeric month

#############################################
# SECTION 3: Plot with Quantile Confidence Intervals by zone and month and exposure
#############################################

# 3.1 Theoretical vs Simulated R0 with quantile CI
p_basic_ZONE_theorical_simulated <- ggplot(df_zone, aes(x = MOIS)) +
  geom_smooth(aes(y = R0_basic, color = ZONE, linetype = "Theoretical"), method = "loess", se = FALSE, linewidth = 1) +
  geom_smooth(aes(y = R0_basic_sim_mean, color = ZONE, linetype = "Simulated"), method = "loess", se = FALSE, linewidth = 1) +
  geom_smooth(aes(y = R0_basic_low_quantile, color = ZONE, linetype = "CI Lower"), method = "loess", se = FALSE, linewidth = 0.8) +
  geom_smooth(aes(y = R0_basic_high_quantile, color = ZONE, linetype = "CI Upper"), method = "loess", se = FALSE, linewidth = 0.8) +
  geom_hline(aes(yintercept = log1p(1), linetype = "R0 = 1"), color = "black", size = 1) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200)) +
  labs(
    title = "Theoretical and Simulated R0 Basic with 95% CI",
    subtitle = "Theoretical = Solid, Simulated = Dotted, CI = Dashed",
    x = "Month", y = "R0 value (log1p scale)",
    color = "Zone", linetype = "Line type"
  ) +
  scale_linetype_manual(values = c("Theoretical" = "solid", "Simulated" = "dotted", "CI Lower" = "dashed", "CI Upper" = "dashed", "R0 = 1" = "solid")) +
  theme_minimal() +
  theme(legend.position = "top") +
  facet_grid(expo ~ ZONE, scales = "free_y", labeller = labeller(expo = function(x) paste("Expo:", x)))

# 3.2 Theoretical R0 only
p_basic_ZONE_IC_R0_theorical <- ggplot(df_zone, aes(x = MOIS)) +
  geom_smooth(aes(y = R0_basic, color = ZONE, linetype = "Theoretical"), method = "loess", se = FALSE, linewidth = 1) +
  geom_hline(aes(yintercept = log1p(1), linetype = "R0 = 1"), color = "black", size = 1) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200)) +
  labs(title = "Theoretical R0 Basic", x = "Month", y = "R0 value (log1p scale)", color = "Zone", linetype = "Line type") +
  scale_linetype_manual(values = c("Theoretical" = "solid", "R0 = 1" = "solid")) +
  theme_minimal() +
  theme(legend.position = "top") +
  facet_grid(expo ~ ZONE, scales = "free_y", labeller = labeller(expo = function(x) paste("Expo:", x)))

# Combine and save plots
R0_ZONE_basic_quantile <- arrangeGrob(p_basic_ZONE_theorical_simulated, p_basic_ZONE_IC_R0_theorical, ncol = 2)
plot(R0_ZONE_basic_quantile)
ggsave("plots/R0_ZONE_basic_quantile_DENV.pdf", R0_ZONE_basic_quantile, width = 16.5, height = 11.7, units = "in")


#############################################
# SECTION 4: Summarized Data with Means and Confidence Intervals by area and month and exposure
#############################################

df_AREA<-df_fin%>%
  group_by(AREA, ZONE, MOIS, expo)%>%
  summarise(ma=mean(ma, na.rm=T), 
            parity_rate=mean(parity_rate, na.rm=T),
            R0_basic=mean(R0_basic, na.rm=T),
            p=mean(p, na.rm=T),
            R0_basic_sim_mean=mean(R0_basic_sim_mean, na.rm=T),
            R0_basic_low_quantile =mean(R0_basic_low_quantile  , na.rm=T),
            R0_basic_high_quantile= mean(R0_basic_high_quantile, na.rm=T), 
            R0_poletti_r0 = mean(R0_poletti_r0 , na.rm=T),
            R0_poletti_sim_mean= mean(R0_poletti_sim_mean, na.rm=T),
            R0_poletti_low_quantile= mean(R0_poletti_low_quantile , na.rm=T),
            R0_poletti_high_quantile= mean(R0_poletti_high_quantile, na.rm=T))


#############################################
# SECTION 5: Raw Boxplot Visualizations
#############################################

# 6.1 Theoretical R0 (basic)
p_basic <- df_fin %>%
  ggplot(aes(x = as.factor(MOIS), y = R0_basic, fill = AREA)) +
  geom_boxplot() +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  facet_grid(expo ~ AREA,
             scales = "free_y",
             labeller = labeller(expo = function(x) paste("Expo:", x))) +
  xlab("MOIS") +
  ggtitle("Theorical R0 basic according to month, sampling areas and human exposure")

# 6.2 Theoretical R0 (POletti)
p_poletti<-df_fin %>%
  ggplot(aes(x = as.factor(MOIS), y = R0_poletti_r0, fill = AREA)) +
  geom_boxplot() +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  facet_grid(expo ~ AREA,
             scales = "free_y",
             labeller = labeller(expo = function(x) paste("Expo:", x))) +
  xlab("MOIS") +
  ggtitle("Theorical R0 poletti according to month, sampling areas and human exposure")


# 6.3 Simulated  R0 (basic)
p_basic_sim_mean<-df_fin %>%
  ggplot(aes(x = as.factor(MOIS), y= R0_basic_sim_mean, fill = AREA)) +
  geom_boxplot() +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  facet_grid(expo ~ AREA,
             scales = "free_y",
             labeller = labeller(expo = function(x) paste("Expo:", x))) +
  xlab("MOIS") +
  ggtitle("Simulated R0 basic according to month, sampling areas and human exposure")

# 6.4 Simulated  R0 (POletti)
p_poletti_sim_mean<-df_fin%>%
  ggplot(aes(x = as.factor(MOIS), y= R0_poletti_sim_mean, fill = AREA))  +
  geom_boxplot() +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  facet_grid(expo ~ AREA,
             scales = "free_y",
             labeller = labeller(expo = function(x) paste("Expo:", x))) +
  xlab("MOIS") +
  ggtitle("Simulated  R0 poletti according to month, sampling areas and human exposure")

df_AREA <- df_AREA %>%
  mutate(MOIS = as.numeric(as.character(MOIS)))  


#############################################
# SECTION 6: Theorical and Simulated R0 Poletti - With Confidence Intervals (Quantiles) by area, month and exposure
#############################################

p_poletti_AREA_IC <- ggplot(data = df_AREA, aes(x = MOIS)) +
  geom_smooth(aes(y = R0_poletti_r0, color = ZONE, linetype = "Theorical"),
              method = "loess", se = FALSE, linewidth = 1) +
  geom_smooth(aes(y = R0_poletti_sim_mean, color = ZONE, linetype = "Simulated"),
              method = "loess", se = FALSE, linewidth = 1) +
  geom_smooth(aes(y = R0_poletti_low_quantile, color = ZONE, linetype = "IC (Lower)"),
              method = "loess", se = FALSE, linewidth = 0.8) +
  geom_smooth(aes(y = R0_poletti_high_quantile, color = ZONE, linetype = "IC (Upper)"),
              method = "loess", se = FALSE, linewidth = 0.8) +
  geom_hline(aes(yintercept = log1p(1), linetype = "R0 = 1"), 
             color = "black", size = 1, show.legend = TRUE) +
  scale_y_continuous(
    trans = "log1p",
    breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200),
    labels = c(0, 1, 2, 5, 10, 20, 50, 100, 200)
  ) +
  labs(
    title = "Theorical and Simulated R0 Poletti with 95% CIs (Quantiles)",
    subtitle = "Solid = Theorical, Dotted = Simulated, Dashed = IC, R0 = 1",
    x = "Month", y = "R0 Poletti (log1p scale)",
    color = "AREAS", linetype = "Lines"
  ) +
  scale_linetype_manual(
    values = c(
      "Theorical" = "solid",
      "Simulated" = "dotted",
      "IC (Lower)" = "dashed",
      "IC (Upper)" = "dashed",
      "R0 = 1" = "solid"
    )
  ) +
  theme_minimal() +
  theme(legend.position = "top") +
  facet_grid(expo ~ AREA, scales = "free_y", 
             labeller = labeller(expo = function(x) paste("Expo:", x)))



#############################################
# SECTION 7: Combined Poletti Plots and Export - With Confidence Intervals (Quantiles) by area, month and exposure
#############################################

ggsave("plots/R0_AREA_poletti_DENV.pdf", p_poletti_AREA_IC,
       width = 16.5, height = 11.7, units = "in")

