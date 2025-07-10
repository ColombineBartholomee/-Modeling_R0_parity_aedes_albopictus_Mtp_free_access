# app.R
library(shiny) # ‘1.9.1’
library(tidyverse) # ‘2.0.0’
library(plotly) # ‘4.10.4’
library(purrr) # ‘1.0.1’
library(mgcv)# ‘1.9.1’
library(boot)# ‘1.3.31’
# ---- Load your functions and data ----
source("analysis/functions_R0.R")

df <- read.csv("data/df_model.csv") %>%
  mutate(NB_TOT = NB_ALBO_PP + NB_ALBO_F_NP,
         DATE_COLLECTE = parse_date(DATE_COLLECTE, "%Y-%m-%d"),
         PROP_PP = ifelse(NB_TOT > 0, NB_ALBO_PP / NB_TOT, 0),
         AREA = as.factor(AREA),
         ZONE = as.factor(ZONE)) %>%
  na.omit()

# ---- Helper Functions ----
compute_r0_from_df <- function(df_input, virus, expo, method_N, method_c, method_b) {
  df<-df_input %>%
    mutate(
      TMN = TMEAN_collection,
      ma = NB_ALBO_F,
      parity_rate = PROP_PP,
      expo = expo,
      b = pmap_dbl(list(virus, TMN), ~compute_b(..1, ..2, method=method_b)),
      c = pmap_dbl(list(virus),~compute_c(..1,  method=method_c)),
      r = map_dbl(virus, compute_r),
      N = pmap_dbl(list(virus, TMN), ~compute_N(..1, ..2, method=method_N)),
      phi = compute_phi(),
      gc = map_dbl(TMN, compute_gc),
      a = phi / gc,
      p = map2_dbl(parity_rate, gc, compute_p),
      capacity = pmap_dbl(list(ma, a, p, expo, N), ~compute_capacity(..1, ..2, ..3, ..4, ..5)),
      R0_basic = pmap_dbl(list(capacity, r, b, c), ~compute_R0(..1, ..2, ..3, ..4)),
      R0_poletti = pmap(list(a, b, c, r, N, p, expo, ma), compute_R0_poletti)
    ) %>%
    unnest_wider(R0_poletti, names_sep = "_")
  return(df)
}

simulate_r0 <- function(df_out, virus, expo, nsim = 1000) {
  df_parity <- compute_parity_rate_ic(df_out, method = "area_mois", nsim = nsim)
  df_ma <- compute_ma_counts(df_out, method = "area_mois", nsim = nsim) %>%
    dplyr::select(ID_COLLECTE, sim_id, ma_sim)
  
  df_base <- df_out %>%
    select(ID_COLLECTE, ZONE, AREA, MOIS, ma, parity_rate, TMN, b, c, r, N, gc, phi, a, p, expo, capacity,
           R0_basic, R0_poletti_r0, R0_poletti_r0VH, R0_poletti_r0HV)
  
  
  df<-df_parity %>%
    left_join(df_ma, by = c("ID_COLLECTE", "sim_id")) %>%
    left_join(df_base, by = "ID_COLLECTE", suffix = c("", "")) %>%
    mutate(
      expo = expo,
      p_sim = map2_dbl(parity_rate_sim, gc, compute_p),
      capacity_sim = pmap_dbl(list(ma_sim, a, p_sim, expo, N), ~compute_capacity(..1, ..2, ..3, ..4, ..5)),
      R0_basic_sim = pmap_dbl(list(capacity_sim, r, b, c), ~compute_R0(..1, ..2, ..3, ..4)),
      R0_poletti_sim = pmap(list(a, b, c, r, N, p_sim, expo, ma_sim), compute_R0_poletti),
      R0_poletti_sim_r0 = map_dbl(R0_poletti_sim, "r0"),
      R0_poletti_sim_r0VH = map_dbl(R0_poletti_sim, "r0VH"),
      R0_poletti_sim_r0HV = map_dbl(R0_poletti_sim, "r0HV")
    )
  return(df)
}

summarise_simulated_r0 <- function(df_sim) {
  df<-df_sim %>%
    group_by(ID_COLLECTE, ZONE, AREA, MOIS, ma, parity_rate, TMN, b, c, r, N, gc, phi, a, p, expo, capacity, R0_basic, R0_poletti_r0, R0_poletti_r0VH, R0_poletti_r0HV) %>%
    summarise(
      ma_sim=mean(ma_sim, na.rm=T), 
      parity_rate_sim=mean(parity_rate_sim, na.rm=T), 
      p_sim=mean(p_sim, na.rm=T),
      capacity = mean(capacity, na.rm = TRUE),
      R0_basic_sim_mean = mean(R0_basic_sim, na.rm = TRUE),
      R0_poletti_sim_mean = mean(R0_poletti_sim_r0, na.rm = TRUE),
      R0_basic_low_quantile = quantile(R0_basic_sim, 0.025, na.rm = TRUE),
      R0_basic_high_quantile = quantile(R0_basic_sim, 0.975, na.rm = TRUE),
      R0_poletti_low_quantile = quantile(R0_poletti_sim_r0, 0.025, na.rm = TRUE),
      R0_poletti_high_quantile = quantile(R0_poletti_sim_r0, 0.975, na.rm = TRUE),
      .groups = 'drop'
    ) 
  return(df)
}

# ---- UI ----
ui <- fluidPage(
  titlePanel("Simulation de R0 - Arbovirus (DENV, ZIKV, CHIKV)"),
  sidebarLayout(
    sidebarPanel(
      selectInput("virus", "Choisir le virus :", choices = c("DENV", "ZIKA", "CHIKV")),
      selectInput("method_N", "Méthode pour N", choices = NULL),
      selectInput("method_c", "Méthode pour c", choices = NULL),
      selectInput("method_b", "Méthode pour b", choices = NULL),
      selectInput("grouping", "Grouper par :", choices = c("ZONE", "AREA"), selected = "ZONE"),
      checkboxGroupInput("zone_area", "Choisir les groupes :", choices = unique(c(levels(df$ZONE), levels(df$AREA))), selected = unique(df$ZONE)),
      sliderInput("expo", "Niveau d'exposition :", min = 0.1, max = 1, step = 0.05, value = 1),
      actionButton("run_model", "Lancer la simulation", icon = icon("play"))
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("R0", plotlyOutput("r0Plot")),
        tabPanel("Détails biologiques",
                 plotlyOutput("plot_ma"),
                 plotlyOutput("plot_p"),
                 plotlyOutput("plot_parity")),
        tabPanel("Détails environnementaux",
                 plotlyOutput("plot_gc"),
                 plotlyOutput("plot_N"),
                 plotlyOutput("plot_temp"))
      )
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {
  methods_list <- list(
    DENV = list(
      method_N = c("Caminade2016", "Benkimoun2021"),
      method_c = c("Metelman2019"),
      method_b = c("Verga-Rua2013", "Benkimoun2021")
    ),
    ZIKA = list(
      method_N = c("Caminade2016","Winokur2020","Tesla2018"),
      method_c = c("Caminade2016","Tesla2018"),
      method_b = c("DiLuca2016")
    ),
    CHIKV = list(
      method_N = c("Christofferson2023"),
      method_c = c("Solimini2018"),
      method_b = c("Verga-Rua2013","method name")
    )
  )
  
  # Mettre à jour les selectInput des méthodes dès le lancement de l'app
  observe({
    virus <- input$virus
    if (is.null(virus) || virus == "") return()
    updateSelectInput(session, "method_N",
                      choices = methods_list[[virus]]$method_N,
                      selected = methods_list[[virus]]$method_N[1])
    updateSelectInput(session, "method_c",
                      choices = methods_list[[virus]]$method_c,
                      selected = methods_list[[virus]]$method_c[1])
    updateSelectInput(session, "method_b",
                      choices = methods_list[[virus]]$method_b,
                      selected = methods_list[[virus]]$method_b[1])
  })
  
  # Afficher les sélections pour test/debug
  output$selected <- renderPrint({
    list(
      virus = input$virus,
      method_N = input$method_N,
      method_c = input$method_c,
      method_b = input$method_b
    )
  })
  
  # Mise à jour dynamique des choix du checkboxGroup selon le grouping choisi
  observeEvent(input$grouping, {
    req(input$grouping)
    choices <- if (input$grouping == "ZONE") levels(df$ZONE) else levels(df$AREA)
    updateCheckboxGroupInput(session, "zone_area", choices = choices, selected = choices)
  })
  
  # Filtrage réactif des données selon les sélections
  filtered_data <- reactive({
    req(input$zone_area, input$grouping)
    
    df %>%
      filter(MOIS %in% 5:11) %>%
      filter(!!rlang::sym(input$grouping) %in% input$zone_area)
  })
  
  r0_results <- eventReactive(input$run_model, {
    req(input$grouping)
    grouping_var <- input$grouping
    group_vars <- c(grouping_var, "MOIS", "expo")
    
    withProgress(message = "Calcul en cours...", value = 0, {
      incProgress(0.2, detail = "Filtrage des données")
      Sys.sleep(0.1)
      
      df_input <- filtered_data()
      virus <- input$virus
      expo_val <- input$expo
      
      incProgress(0.4, detail = "Calcul des valeurs théoriques")
      Sys.sleep(0.1)
      
      df_out <-compute_r0_from_df(
        df_input,
        virus = virus,
        expo = expo_val,
        method_N = input$method_N,
        method_c = input$method_c,
        method_b = input$method_b
      )
      
      incProgress(0.6, detail = "Simulation")
      Sys.sleep(0.1)
      
      df_sim <- simulate_r0(df_out, virus = virus, expo = expo_val)
      
      incProgress(0.8, detail = "Résumé des résultats")
      Sys.sleep(0.1)
      
      df_res <- summarise_simulated_r0(df_sim)
      
      df_agg <- df_res %>%
        group_by(across(all_of(group_vars))) %>%
        summarise(
          gc=mean(gc, na.rm=T),
          N=mean(N, na.rm=T),
          TMN=mean(TMN, na.rm=T),
          ma = mean(ma, na.rm = TRUE),
          parity_rate = mean(parity_rate, na.rm = TRUE),
          R0_basic = mean(R0_basic, na.rm = TRUE),
          p = mean(p, na.rm = TRUE),
          ma_sim=mean(ma_sim, na.rm=T), 
          parity_rate_sim=mean(parity_rate_sim, na.rm=T), 
          p_sim=mean(p_sim, na.rm=T),
          R0_basic_sim_mean = mean(R0_basic_sim_mean, na.rm = TRUE),
          R0_basic_low_quantile = mean(R0_basic_low_quantile, na.rm = TRUE),
          R0_basic_high_quantile = mean(R0_basic_high_quantile, na.rm = TRUE),
          R0_poletti_r0 = mean(R0_poletti_r0, na.rm = TRUE),
          R0_poletti_sim_mean = mean(R0_poletti_sim_mean, na.rm = TRUE),
          R0_poletti_low_quantile = mean(R0_poletti_low_quantile, na.rm = TRUE),
          R0_poletti_high_quantile = mean(R0_poletti_high_quantile, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(MOIS = as.numeric(as.character(MOIS)))
      
      incProgress(1, detail = "Terminé")
      df_agg
    })
  })
  
  output$r0Plot <- renderPlotly({
    req(r0_results())
    df_plot <- r0_results()
    grouping_var <- input$grouping
    ci_type <- input$r0_ci_type
    
    group_sym <- sym(grouping_var)
    # IC quantile
    p <- ggplot(df_plot, aes(x = MOIS)) +
      # R0 simulé
      geom_smooth(aes(y = R0_basic_sim_mean, color = !!group_sym, linetype = "Simulé"),
                  method = "loess", se = FALSE, linewidth = 1) +
      geom_smooth(aes(y = R0_basic_low_quantile, color = !!group_sym, linetype = "IC inférieur"),
                  method = "loess", se = FALSE, linewidth = 0.8) +
      geom_smooth(aes(y = R0_basic_high_quantile, color = !!group_sym, linetype = "IC supérieur"),
                  method = "loess", se = FALSE, linewidth = 0.8) +
      # R0 théorique (non simulé), avec smooth
      geom_smooth(aes(y = R0_basic, color = !!group_sym, linetype = "Théorique"),
                  method = "loess", se = FALSE, linewidth = 1, alpha = 0.7) +
      geom_hline(yintercept = 1, linetype = "solid", color = "black", size = 1) +
      scale_y_continuous(trans = "log1p", breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200)) +
      labs(
        title = paste("R0 Basic simulé et théorique avec IC quantile pour", input$virus),
        subtitle = "Simulé et Théorique = traits pleins, IC = traits pointillés",
        x = "Mois", y = "Valeur de R0",
        color = grouping_var,
        linetype = "Type de ligne"
      ) +
      scale_linetype_manual(values = c(
        "Simulé" = "dashed",
        "Théorique" = "solid",
        "IC inférieur" = "dotted",
        "IC supérieur" = "dotted"
      )) +
      theme_minimal() +
      theme(legend.position = "top") +
      facet_wrap(as.formula(paste("~", grouping_var)), scales = "free_y")
    
    
    ggplotly(p, height = 700, width = NULL)
  })
  
  output$plot_ma <- renderPlotly({
    req(r0_results())
    df_plot <- r0_results()
    grouping_var <- input$grouping
    
    p <- ggplot(df_plot, aes_string(x = "MOIS")) +
      geom_smooth(aes_string(y = "ma_sim", color = grouping_var, linetype = "'Simulé'"),
                  method = "loess", se = FALSE, linewidth = 1) +
      geom_smooth(aes_string(y = "ma", color = grouping_var, linetype = "'Théorique'"),
                  method = "loess", se = FALSE, linewidth = 1) +
      labs(title = "Nombre moyen moustiques femelles (ma)",
           x = "Mois", y = "ma",
           color = grouping_var, linetype = "Type de ligne") +
      scale_linetype_manual(values = c("Théorique" = "solid", "Simulé" = "dashed")) +
      theme_minimal() +
      theme(legend.position = "top") +
      facet_wrap(as.formula(paste("~", grouping_var)), scales = "free_y")+
      ylim(c(0,70))
    
    ggplotly(p)
  })
  
  output$plot_p <- renderPlotly({
    req(r0_results())
    df_plot <- r0_results()
    grouping_var <- input$grouping
    
    p <- ggplot(df_plot, aes_string(x = "MOIS")) +
      geom_smooth(aes_string(y = "p_sim", color = grouping_var, linetype = "'Simulé'"),
                  method = "loess", se = FALSE, linewidth = 1) +
      geom_smooth(aes_string(y = "p", color = grouping_var, linetype = "'Théorique'"),
                  method = "loess", se = FALSE, linewidth = 1) +
      labs(title = "Taux de survie quotidien (p)",
           x = "Mois", y = "p",
           color = grouping_var, linetype = "Type de ligne") +
      scale_linetype_manual(values = c("Théorique" = "solid", "Simulé" = "dashed")) +
      theme_minimal() +
      theme(legend.position = "top") +
      facet_wrap(as.formula(paste("~", grouping_var)), scales = "free_y")+
      ylim(c(0,1))
    
    ggplotly(p)
  })
  
  output$plot_parity <- renderPlotly({
    req(r0_results())
    df_plot <- r0_results()
    grouping_var <- input$grouping
    
    p <- ggplot(df_plot, aes_string(x = "MOIS")) +
      geom_smooth(aes_string(y = "parity_rate_sim", color = grouping_var, linetype = "'Simulé'"),
                  method = "loess", se = FALSE, linewidth = 1) +
      geom_smooth(aes_string(y = "parity_rate", color = grouping_var, linetype = "'Théorique'"),
                  method = "loess", se = FALSE, linewidth = 1) +
      labs(title = "Taux de parité",
           x = "Mois", y = "Taux de parité",
           color = grouping_var, linetype = "Type de ligne") +
      scale_linetype_manual(values = c("Théorique" = "solid", "Simulé" = "dashed")) +
      theme_minimal() +
      theme(legend.position = "top") +
      facet_wrap(as.formula(paste("~", grouping_var)), scales = "free_y")+
      ylim(c(0,1))
    
    ggplotly(p)
  })
  
  output$plot_gc <- renderPlotly({
    req(r0_results())
    df_plot <- r0_results()
    grouping_var <- input$grouping
    
    p <- ggplot(df_plot, aes_string(x = "MOIS")) +
      geom_smooth(aes_string(y = "gc", color = grouping_var),
                  method = "loess", se = FALSE, linewidth = 1) +
      labs(title = "Cycle gonotrophique",
           x = "Mois", y = "Cycle gonotrophique",
           color = grouping_var) +
      theme_minimal() +
      theme(legend.position = "top") +
      facet_wrap(as.formula(paste("~", grouping_var)), scales = "free_y")
    
    ggplotly(p)
  })
  
  output$plot_N <- renderPlotly({
    req(r0_results())
    df_plot <- r0_results()
    grouping_var <- input$grouping
    
    p <- ggplot(df_plot, aes_string(x = "MOIS")) +
      geom_smooth(aes_string(y = "N", color = grouping_var),
                  method = "loess", se = FALSE, linewidth = 1) +
      labs(title = "EIP",
           x = "Mois", y = "EIP",
           color = grouping_var) +
      theme_minimal() +
      theme(legend.position = "top") +
      facet_wrap(as.formula(paste("~", grouping_var)), scales = "free_y")
    
    ggplotly(p)
  })
  
  output$plot_temp <- renderPlotly({
    req(r0_results())
    df_plot <- r0_results()
    grouping_var <- input$grouping
    
    p <- ggplot(df_plot, aes_string(x = "MOIS")) +
      geom_smooth(aes_string(y = "TMN", color = grouping_var),
                  method = "loess", se = FALSE, linewidth = 1) +
      labs(title = "Temperature moyenne ",
           x = "Mois", y = "Temperature",
           color = grouping_var) +
      theme_minimal() +
      theme(legend.position = "top") +
      facet_wrap(as.formula(paste("~", grouping_var)), scales = "free_y")
    
    ggplotly(p)
  })
}

shinyApp(ui = ui, server = server)


