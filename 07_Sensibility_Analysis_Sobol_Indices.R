# R Script: Sensitivity Analysis for a Vector-Borne Pathogen Transmission Model
# This script performs a Sobol sensitivity analysis to evaluate the impact of different parameters 
# on the reproduction and transmission of a pathogen transmitted by a vector. 
# It includes the calculation of vector competence, basic reproduction rate (R0), 
# and sensitivity analysis for these parameters.

######### Load necessary libraries
library(tidyverse) # ‘2.0.0’
library(sensitivity) #  ‘1.30.1’
library(ggplot2) #  ‘3.5.1’
rm(list=ls())
######### Functions used
source("analysis/functions_R0.R")

# -------------------------------------
# Base parameter ranges 
# -------------------------------------
base_param_ranges <- list(
  ma = c(1, 65),
  T = c(10, 45),
  parity = c(0.1, 0.90),
  c = c(0.14, 0.39),
  r = c(2, 7), 
  phi=c(0.684, 1),
  expo=c(0.1, 1),
  a=c(0.078, 0.3), 
  gc=c(2.90, 11),
  b = c(0.1, 10.8)
)

# -------------------------------------
# Define combinations of methods
# -------------------------------------
combinations <- expand.grid(
  method_N = c("Caminade2016", "Benkimoun2021"),
  method_b = c("Caminade2016", "Benkimoun2021", "Verga-Rua2013", "Metelman2019"),
  method_R0 = c("R0_basic", "R0_poletti"),
  stringsAsFactors = FALSE
)

# -------------------------------------
# Run Sobol analysis for each combination
# -------------------------------------
all_sobol_results <- list()

for (i in 1:nrow(combinations)) {
  combo <- combinations[i, ]
  cat("\n➡️  Combo", i, "/", nrow(combinations), ":", paste(combo, collapse = " / "), "\n")
  
  # Run Sobol sensitivity analysis for this combination
  sobol_res <- run_sensitivity_for_combo_sobol(
    method_N = combo$method_N,
    method_b = combo$method_b,
    method_R0 = combo$method_R0,
    n = 100000
  )
  
  # Store results
  all_sobol_results[[i]] <- list(
    methods = combo,
    sobol = sobol_res
  )
}

# Combine Sobol results
sobol_summary <- bind_rows(lapply(all_sobol_results, function(x) x$sobol))%>%
  mutate(
    first_order = ifelse(first_order < 0, 0, first_order),
    total_order = ifelse(total_order < 0, 0, total_order)
  )

# -------------------------------------
# Plotting results
# -------------------------------------
# Calculate the total number of rows in the resulting dataframe
n_sobol_rows <- nrow(sobol_summary)

# Add method descriptions to the results
sobol_summary <- sobol_summary %>%
  mutate(
    # Repeat the method values for each row of Sobol results
    method_N = rep(sapply(all_sobol_results, function(x) x$methods$method_N), each = n_sobol_rows / length(all_sobol_results)),
    method_b = rep(sapply(all_sobol_results, function(x) x$methods$method_b), each = n_sobol_rows / length(all_sobol_results)),
    method_R0 = rep(sapply(all_sobol_results, function(x) x$methods$method_R0), each = n_sobol_rows / length(all_sobol_results))
  )

# Display the final summary of results
print(sobol_summary)

# Descriptions for 'method_b'
method_b_desc <- tibble::tibble(
  method_b = c("Verga-Rua2013", "Metelman2019", "Caminade2016", "Benkimoun2021"),
  desc_b = c(
    "Random draw (uniform 0.56–0.67)",
    "Constant (b = 0.5)",
    "Random draw (uniform 0.1–0.75)",
    "Function of T (Benkimoun 2021)"
  )
)

# Descriptions for 'method_N'
method_N_desc <- tibble::tibble(
  method_N = c("Caminade2016", "Benkimoun2021"),
  desc_N = c(
    "Function of T (Caminade 2016)",
    "Function of T (Benkimoun 2021)"
  )
)

# Join descriptions to the summary dataframe
sobol_summary <- sobol_summary %>%
  left_join(method_b_desc, by = "method_b") %>%
  left_join(method_N_desc, by = "method_N") %>%
  mutate(
    # Create a combined description for 'desc_b_N'
    desc_b_N = paste("N:", desc_N, "\nb:", desc_b),
    
    # Apply the absolute value to the first and second order indices
    first_order = abs(first_order), 
    total_order = abs(total_order)
  )


# Plot the graph 
plot_sobol<-sobol_summary %>%
  gather(key = "order_type", value = "index_value", first_order, total_order) %>%
  ggplot(aes(x = parameter, y = index_value, fill = order_type)) +
  geom_bar(stat = "identity", position = "dodge") +  # Bar chart side by side
  theme_minimal() + 
  facet_wrap(~method_R0 + desc_b_N) +  # Facets for each combination
  labs(
    title = "Sobol Indices: First and Second Order",
    x = "Parameters",
    y = "Sensitivity Index",
    fill = "Order Type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate X-axis labels
    strip.text = element_text(size = 10)  # Adjust the size of facet labels
  ) +
  scale_y_continuous(
    limits = c(0, 1),  # Set Y-axis limits to (0, 1)
    expand = c(0, 0)    # Avoid extra space at the bottom
  ) 

ggsave("plots/sobol_indices.pdf", plot_sobol,
       width = 16.5, height = 11.7, units = "in")

