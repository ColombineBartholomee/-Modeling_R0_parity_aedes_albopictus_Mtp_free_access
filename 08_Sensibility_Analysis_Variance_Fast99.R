###################################### This script performs sensitivity analysis using FAST99 functions
###################################### for different parameters and calculates model variance 
###################################### to assess which model is the most robust 
############################### Regarding the calculation of R0.

######### Packages used
rm(list=ls())
library(ggplot2)  #  ‘3.5.1’
library(tidyverse) # ‘2.0.0’
library(sensitivity) #  ‘1.30.1’

######### Functions used
source("analysis/functions_R0.R")

######### R0 analysis for the DENV model with the following parameters depending on temperature

# Base parameter ranges (without b)
base_param_ranges <- data.frame(
  min = c(10, 0.14, 2, 0.1, 1, 0.1,  0.682),  # T, c, r, parity, ma, human exposure,aedes exposure,  human trophic preference
  max = c(45, 0.39, 7, 0.90, 65, 1, 1)
)
rownames(base_param_ranges) <- c("T", "c", "r", "parity", "ma", "expo", "phi")

# -------------------------------------
# Define combinations
# -------------------------------------
combinations <- expand.grid(
  method_N = c("Caminade2016", "Benkimoun2021"),
  method_b = c("Caminade2016", "Benkimoun2021", "Verga-Rua2013", "Metelman2019"),
  method_R0 = c("R0_basic", "R0_poletti"),
  stringsAsFactors = FALSE
)

# -------------------------------------
# Store results
# -------------------------------------
all_results <- list()

# -------------------------------------
# Run Fast 99 analysis for each combination
# -------------------------------------
for (i in 1:nrow(combinations)) {
  combo <- combinations[i, ]
  cat("\n➡️  Combo", i, "/", nrow(combinations), ":", paste(combo, collapse = " / "), "\n")
  res <- run_sensitivity_fast99(
    method_N = combo$method_N,
    method_b = combo$method_b,
    method_R0 = combo$method_R0,
    n = 1000
  )
  all_results[[i]] <- list(
    methods = combo,
    summary = res$summary,
    fast = res$fast_res
  )
}
# Combine  results
summary_df <- do.call(rbind, lapply(all_results, function(res) {
  data.frame(
    method_N = res$methods$method_N,
    method_b = res$methods$method_b,
    method_R0 = res$methods$method_R0,
    param = names(res$fast$X),
    S_main = res$fast$D1 / res$fast$V,
    S_total_normalized = res$fast$Dt / res$fast$V
  )
}))


# Combine results for variance
variance_summary <- do.call(rbind, lapply(all_results, function(res) {
  data.frame(
    method_N = res$methods$method_N,
    method_b = res$methods$method_b,
    method_R0 = res$methods$method_R0,
    variance_R0 = res$fast$V
  )
}))


# -------------------------------------
# Plotting results
# -------------------------------------

# Combine  results by methods used and linked with the description of each method
summary_comparaison <- summary_df %>%
  group_by(method_N, method_b, method_R0, param) %>%
  summarise(
    mean_S_main = mean(S_main, na.rm = TRUE),
    mean_S_total = mean(S_total_normalized, na.rm = TRUE)
  )

method_b_desc <- tibble::tibble( ## Description of method used for b
  method_b = c("Verga-Rua2013", "Metelman2019", "Caminade2016", "Benkimoun2021"),
  desc_b = c(
    "Random draw (uniform 0.56–0.67)",
    "Constant (b = 0.5)",
    "Random draw (uniform 0.1–0.75)",
    "Function of T (Benkimoun 2021)"
  )
)

method_N_desc <- tibble::tibble( ## Description of method used for N
  method_N = c("Caminade2016", "Benkimoun2021"),
  desc_N = c(
    "Function of T (Caminade 2016)",
    "Function of T (Benkimoun 2021)"
  )
)

summary_comparaison <- summary_comparaison %>%
  left_join(method_b_desc, by = "method_b") %>%
  left_join(method_N_desc, by = "method_N") %>%
  mutate(
    desc_b_N = paste("N:", desc_N,"\nb:", desc_b)
  )


# # Plot the graph with the indices for each parameter, by methods

p_fast99<-ggplot(summary_comparaison, aes(x = param, y = mean_S_main, fill = param)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~method_R0  + desc_b_N) +
  labs(
    x = "Parameter", 
    y = "Main sensitivity (S_main)", 
    fill = "N method",
    title = "Sensitivity according to methods for b, N and R0"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("plots/fast99.pdf", p_fast99,
       width = 16.5, height = 11.7, units = "in")
# -------------------------------------
# Calcul of variance
# -------------------------------------

## Selection of combination with 5 lower variances
top_5_variances <- variance_summary %>%
  arrange(variance_R0) %>%
  head(5)

## Looking to combination with lower variance
plot_variance<-ggplot(top_5_variances, aes(x = reorder(paste(method_N, method_b, method_R0, sep = " | "), variance_R0), y = variance_R0, fill=method_N)) +
  geom_col(fill = "#69b3a2") +
  coord_flip()+
  labs(
    title = "Top 10 combinations with the lowest variance",
    x = "Method (N | b | R0)",
    y = "Total variance (V)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 10)
  ) +
  facet_wrap(~method_R0)

