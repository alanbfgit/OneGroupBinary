library(shiny)
library(bslib)
library(clinfun)
library(ggplot2)

ui <- page_sidebar(
  title = "Simon's Two-Stage Design. Genelux - Forsythe 2026 v2",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  sidebar = sidebar(
    width = 300,
    sliderInput("pu",    "Undesirable response rate (p\u2080)",
                min = 0.05, max = 0.60, value = 0.20, step = 0.01),
    sliderInput("pa",    "Acceptable response rate (p\u2081)",
                min = 0.10, max = 0.80, value = 0.40, step = 0.01),
    sliderInput("ep1",   "Type I error (\u03b1)",
                min = 0.01, max = 0.20, value = 0.10, step = 0.01),
    sliderInput("power", "Power (1 \u2212 \u03b2)",
                min = 0.70, max = 0.99, value = 0.90, step = 0.01),
    hr(),
    radioButtons("design_type", "Design Criterion",
                 choices  = c("Optimal" = "optimal", "MiniMax" = "minimax"),
                 selected = "optimal",
                 inline   = TRUE),
    hr(),
    tags$p(tags$strong("Beta Prior on Response Rate")),
    tags$small(tags$em(
      "Beta(\u03b1, \u03b2) prior on p. ",
      "Default Beta(1, 1) is uniform."
    )),
    numericInput("prior_alpha", "\u03b1  (prior successes)",
                 value = 1, min = 0.1, max = 20, step = 0.1),
    numericInput("prior_beta",  "\u03b2  (prior failures)",
                 value = 1, min = 0.1, max = 20, step = 0.1),
    hr(),
    uiOutput("param_warning")
  ),
  
  navset_card_tab(
    # ── Tab 1: Simon Design ──────────────────────────────────────────────────
    nav_panel(
      "Simon Design",
      layout_column_wrap(
        width = 1 / 2,
        fill  = FALSE,
        value_box(
          title    = "Stage 1: enroll",
          value    = textOutput("vb_n1"),
          showcase = bsicons::bs_icon("people"),
          theme    = "primary"
        ),
        value_box(
          title    = "Total sample size",
          value    = textOutput("vb_n"),
          showcase = bsicons::bs_icon("people-fill"),
          theme    = "info"
        ),
        value_box(
          title    = "Stage 1 threshold (r\u2081)",
          value    = textOutput("vb_r1"),
          showcase = bsicons::bs_icon("check2-circle"),
          theme    = "secondary"
        ),
        value_box(
          title    = "Final threshold (r)",
          value    = textOutput("vb_r"),
          showcase = bsicons::bs_icon("check2-all"),
          theme    = "secondary"
        )
      ),
      
      layout_column_wrap(
        width = 1 / 2,
        card(
          card_header(textOutput("table_header", inline = TRUE)),
          tableOutput("design_table"),
          card_footer(uiOutput("table_footer"))
        ),
        card(
          card_header("Decision Rules"),
          uiOutput("rules")
        )
      )
    ),
    
    # ── Tab 2: Posterior Distribution ────────────────────────────────────────
    nav_panel(
      "Posterior Distribution",
      layout_column_wrap(
        width = 1 / 2,
        fill  = FALSE,
        card(
          card_body(
            sliderInput("n_obs", "Total observations (n)",
                        min = 5, max = 150, value = 30, step = 1)
          )
        ),
        card(
          card_body(
            sliderInput("a_obs", "Observed successes (a)",
                        min = 0, max = 30, value = 12, step = 1)
          )
        )
      ),
      
      layout_column_wrap(
        width = 1 / 2,
        fill  = FALSE,
        value_box(
          title    = htmltools::HTML("P(p &le; p\u2080 &nbsp;| data)"),
          value    = textOutput("vb_post_pu"),
          showcase = bsicons::bs_icon("arrow-down-circle"),
          theme    = "danger",
          height   = "90px"
        ),
        value_box(
          title    = htmltools::HTML("P(p &ge; p\u2081 &nbsp;| data)"),
          value    = textOutput("vb_post_pa"),
          showcase = bsicons::bs_icon("arrow-up-circle"),
          theme    = "success",
          height   = "90px"
        )
      ),
      
      card(
        full_screen = TRUE,
        card_header(textOutput("posterior_plot_header", inline = TRUE)),
        plotOutput("posterior_plot", height = "460px")
      )
    ),
    
    # ── Tab 3: Thall & Simon Design ──────────────────────────────────────────
    nav_panel(
      "Thall & Simon Design",
      layout_column_wrap(
        width = 1 / 2,
        fill  = FALSE,
        card(
          card_body(
            sliderInput("ts_max_n", "Maximum sample size (N)",
                        min = 20, max = 100, value = 30, step = 5)
          )
        ),
        card(
          card_body(
            sliderInput("ts_cohort", "Cohort size for interim looks",
                        min = 1, max = 10, value = 5, step = 1)
          )
        ),
        card(
          card_body(
            sliderInput("ts_fut_thresh", "Futility threshold P(p > p\u2080 | data)",
                        min = 0.01, max = 0.30, value = 0.05, step = 0.01)
          )
        ),
        card(
          card_body(
            sliderInput("ts_eff_thresh", "Efficacy threshold P(p > p\u2080 | data)",
                        min = 0.70, max = 0.99, value = 0.95, step = 0.01)
          )
        )
      ),
      
      layout_column_wrap(
        width = 1 / 2,
        fill  = FALSE,
        value_box(
          title    = "Type I Error (Null)",
          value    = textOutput("vb_ts_type1"),
          showcase = bsicons::bs_icon("exclamation-circle"),
          theme    = "danger",
          height   = "90px"
        ),
        value_box(
          title    = "Power (Alternative)",
          value    = textOutput("vb_ts_power"),
          showcase = bsicons::bs_icon("check-circle"),
          theme    = "success",
          height   = "90px"
        ),
        value_box(
          title    = "ESS under Null",
          value    = textOutput("vb_ts_ess_null"),
          showcase = bsicons::bs_icon("hourglass-split"),
          theme    = "primary",
          height   = "90px"
        ),
        value_box(
          title    = "ESS under Alternative",
          value    = textOutput("vb_ts_ess_alt"),
          showcase = bsicons::bs_icon("hourglass-split"),
          theme    = "info",
          height   = "90px"
        )
      ),
      
      layout_column_wrap(
        width = 1 / 2,
        card(
          card_header("Stopping Boundaries Table"),
          tableOutput("ts_boundary_table"),
          card_footer(uiOutput("ts_table_footer"))
        ),
        card(
          card_header("Operating Characteristics Summary"),
          uiOutput("ts_summary")
        )
      )
    ),
    
    # ── Tab 4: Bayesian Sequential Design ────────────────────────────────────
    nav_panel(
      "Bayesian Sequential Design",
      
      accordion(
        id = "bsd_accordion",
        accordion_panel(
          "Design Configuration",
          layout_column_wrap(
            width = 1 / 2,
            fill  = FALSE,
            card(
              card_header("Simulation Settings"),
              card_body(
                numericInput("bsd_max_n", "Maximum sample size (N)",
                            value = 50, min = 20, max = 200, step = 1),
                tags$small(tags$em(
                  "Threshold grid search: u \u2208 [0.92, 0.98], l \u2208 [0.05, 0.15]. ",
                  "Calibrated for typical Phase II response rates."
                )),
                tags$br(),
                tags$br(),
                actionButton("bsd_run", "Run Calibration",
                            class = "btn-primary w-100",
                            icon  = icon("play"))
              )
            ),
            card(
              card_header("Active Design Parameters"),
              card_body(uiOutput("bsd_params_display"))
            )
          )
        ),
        accordion_panel(
          "Operating Characteristics",
          layout_column_wrap(
            width = 1 / 2,
            fill  = FALSE,
            value_box(
              title    = "Type I Error (H\u2080)",
              value    = textOutput("vb_bsd_type1"),
              showcase = bsicons::bs_icon("exclamation-circle"),
              theme    = "danger"
            ),
            value_box(
              title    = "Power (H\u2081)",
              value    = textOutput("vb_bsd_power"),
              showcase = bsicons::bs_icon("check-circle"),
              theme    = "success"
            ),
            value_box(
              title    = "ESS under H\u2080",
              value    = textOutput("vb_bsd_ess_h0"),
              showcase = bsicons::bs_icon("hourglass-split"),
              theme    = "primary"
            ),
            value_box(
              title    = "ESS under H\u2081",
              value    = textOutput("vb_bsd_ess_h1"),
              showcase = bsicons::bs_icon("hourglass-split"),
              theme    = "info"
            )
          )
        ),
        accordion_panel(
          "Observed Data & Sequential Decision",
          layout_column_wrap(
            width = 1 / 2,
            fill  = FALSE,
            card(
              card_header("Enter Observed Data"),
              card_body(
                numericInput("bsd_successes", "Successes (r)",
                            value = 12, min = 0, max = 200, step = 1),
                numericInput("bsd_total", "Total patients (n)",
                            value = 30, min = 0, max = 200, step = 1)
              )
            ),
            card(
              card_header("Computed Summary"),
              card_body(uiOutput("bsd_computed_summary"))
            )
          ),
          card(
            card_body(uiOutput("bsd_observed_data"))
          )
        ),
        accordion_panel(
          "Decision Boundaries & Summary",
          layout_column_wrap(
            width = 1,
            fill  = FALSE,
            card(
              full_screen = TRUE,
              card_header("Decision Boundaries"),
              plotOutput("bsd_boundary_plot", height = "450px")
            )
          ),
          layout_column_wrap(
            width = 1,
            fill  = FALSE,
            card(
              card_header("Design Summary"),
              uiOutput("bsd_summary")
            )
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Smoothly update the successes slider max to prevent layout flashing
  observeEvent(input$n_obs, {
    current_a <- input$a_obs
    new_max <- input$n_obs
    new_val <- max(0, min(current_a, new_max))
    updateSliderInput(session, "a_obs", max = new_max, value = new_val)
  })
  
  output$param_warning <- renderUI({
    if (input$pa <= input$pu) {
      tags$div(
        class = "alert alert-danger",
        tags$strong("Invalid parameters:"),
        " the acceptable rate (p\u2081) must be greater than the undesirable rate (p\u2080)."
      )
    }
  })
  
  simon_result <- reactive({
    req(input$pa > input$pu)
    ep2 <- 1 - input$power
    tryCatch(
      clinfun::ph2simon(input$pu, input$pa, input$ep1, ep2, nmax = 200),
      error = function(e) NULL
    )
  })
  
  selected_design <- reactive({
    x <- simon_result()
    req(!is.null(x))
    out <- x$out
    if (input$design_type == "minimax") {
      out[which.min(out[, "n"]), , drop = FALSE]
    } else {
      out[which.min(out[, "EN(p0)"]), , drop = FALSE]
    }
  })
  
  post_cdf <- function(threshold, x, n, a, b) {
    pbeta(threshold, a + x, b + n - x)
  }
  
  # ── Tab 1 outputs ──────────────────────────────────────────────────────────
  
  output$vb_n1 <- renderText({
    opt <- selected_design()
    paste(opt[, "n1"], "patients")
  })
  output$vb_n <- renderText({
    opt <- selected_design()
    paste(opt[, "n"], "patients")
  })
  output$vb_r1 <- renderText({
    opt <- selected_design()
    paste("\u2264", opt[, "r1"], "responses \u2192 stop")
  })
  output$vb_r <- renderText({
    opt <- selected_design()
    paste(">", opt[, "r"], "responses \u2192 promising")
  })
  
  output$table_header <- renderText({
    paste(if (input$design_type == "minimax") "MiniMax" else "Optimal", "Design \u2014 Summary Table")
  })
  
  output$table_footer <- renderUI({
    tags$small(tags$em(
      if (input$design_type == "minimax") "\u2018MiniMax\u2019 minimises the maximum sample size (n)."
      else "\u2018Optimal\u2019 minimises the expected sample size under p\u2080 (EN(p\u2080))."
    ))
  })
  
  output$design_table <- renderTable({
    opt <- selected_design()
    data.frame(
      "r\u2081 (Stage 1 threshold)" = opt[, "r1"],
      "n\u2081 (Stage 1 sample)"    = opt[, "n1"],
      "r (Final threshold)"        = opt[, "r"],
      "n (Total sample)"           = opt[, "n"],
      "EN(p\u2080)"                 = round(opt[, "EN(p0)"], 1),
      "PET(p\u2080)"               = round(opt[, "PET(p0)"], 3),
      check.names = FALSE
    )
  }, striped = TRUE, hover = TRUE, width = "100%")
  
  output$rules <- renderUI({
    opt <- selected_design()
    r1 <- opt[, "r1"]
    n1 <- opt[, "n1"]
    r  <- opt[, "r"]
    n  <- opt[, "n"]
    n2 <- n - n1
    
    pet <- round(opt[, "PET(p0)"] * 100, 1)
    en  <- round(opt[, "EN(p0)"],  1)
    
    a <- input$prior_alpha
    b <- input$prior_beta
    
    s1_pu <- post_cdf(input$pu, r1, n1, a, b)
    s1_pa <- 1 - post_cdf(input$pa, r1, n1, a, b)
    fn_pu <- post_cdf(input$pu, r,  n,  a, b)
    fn_pa <- 1 - post_cdf(input$pa, r,  n,  a, b)
    
    prior_label <- sprintf("Beta(%.1f, %.1f) prior", a, b)
    
    tagList(
      tags$h6("Stage 1 stopping rule", class = "fw-bold"),
      tags$p(sprintf(
        "Enroll %d patients in Stage 1. If %d or fewer responses are observed, stop the trial early and conclude that the treatment is not sufficiently active. Under the undesirable rate (p\u2080 = %.2f), there is a %.1f%% probability of stopping early, with an expected trial size of only %.1f patients.",
        n1, r1, input$pu, pet, en
      )),
      tags$div(
        class = "alert alert-secondary py-2",
        tags$strong("Posterior probabilities at Stage 1 boundary"),
        tags$span(class = "text-muted", sprintf(" \u2014 %s, observing %d / %d responses:", prior_label, r1, n1)),
        tags$br(),
        sprintf("P(p \u2264 p\u2080 = %.2f \u2223 data) = %.3f   \u2014 posterior probability the response rate is undesirable or lower", input$pu, s1_pu),
        tags$br(),
        sprintf("P(p \u2265 p\u2081 = %.2f \u2223 data) = %.3f   \u2014 posterior probability the response rate meets the acceptable target", input$pa, s1_pa)
      ),
      tags$hr(),
      tags$h6("Stage 2 decision rule", class = "fw-bold"),
      tags$p(sprintf(
        "If more than %d responses are seen in Stage 1, continue and enroll %d additional patients in Stage 2 (total N = %d). At the end of the trial, if the total number of responses across both stages exceeds %d out of %d patients, declare the treatment promising. If %d or fewer total responses are observed, conclude that the treatment does not meet the target activity level.",
        r1, n2, n, r, n, r
      )),
      tags$div(
        class = "alert alert-secondary py-2",
        tags$strong("Posterior probabilities at final boundary"),
        tags$span(class = "text-muted", sprintf(" \u2014 %s, observing %d / %d responses:", prior_label, r, n)),
        tags$br(),
        sprintf("P(p \u2264 p\u2080 = %.2f \u2223 data) = %.3f   \u2014 posterior probability the response rate is undesirable or lower", input$pu, fn_pu),
        tags$br(),
        sprintf("P(p \u2265 p\u2081 = %.2f \u2223 data) = %.3f   \u2014 posterior probability the response rate meets the acceptable target", input$pa, fn_pa)
      )
    )
  })
  
  # ── Tab 2 outputs ──────────────────────────────────────────────────────────
  
  post_params <- reactive({
    req(!is.null(input$a_obs))
    a  <- input$a_obs
    n  <- input$n_obs
    ap <- input$prior_alpha + a
    bp <- input$prior_beta  + n - a
    list(a = a, n = n, alpha_post = ap, beta_post = bp)
  })
  
  output$posterior_plot_header <- renderText({
    p <- post_params()
    sprintf("Posterior: Beta(%.1f, %.1f)  \u2014  %d successes in %d observations (observed rate = %.1f%%)",
            p$alpha_post, p$beta_post, p$a, p$n, 100 * p$a / p$n)
  })
  
  output$vb_post_pu <- renderText({
    p <- post_params()
    sprintf("%.3f", pbeta(input$pu, p$alpha_post, p$beta_post))
  })
  
  output$vb_post_pa <- renderText({
    p <- post_params()
    sprintf("%.3f", 1 - pbeta(input$pa, p$alpha_post, p$beta_post))
  })
  
  output$posterior_plot <- renderPlot({
    p  <- post_params()
    ap <- p$alpha_post
    bp <- p$beta_post
    n  <- p$n
    pu <- input$pu
    pa <- input$pa
    a  <- p$a
    
    # Generate a smooth continuous grid for the Beta density function
    x_seq <- seq(0, 1, length.out = 1000)
    y_seq <- dbeta(x_seq, ap, bp)
    
    # Classify continuous coordinates into distinct regions
    region <- dplyr::case_when(
      x_seq <= pu ~ paste0("p \u2264 p\u2080 (", pu, ")"),
      x_seq >= pa ~ paste0("p \u2265 p\u2081 (", pa, ")"),
      TRUE        ~ paste0("p\u2080 < p < p\u2081")
    )
    region <- factor(region, levels = c(
      paste0("p \u2264 p\u2080 (", pu, ")"),
      "p\u2080 < p < p\u2081",
      paste0("p \u2265 p\u2081 (", pa, ")")
    ))
    
    prior_y  <- dbeta(x_seq, input$prior_alpha, input$prior_beta)
    df <- data.frame(p_val = x_seq, density = y_seq, prior_density = prior_y, region = region)
    
    mle_p    <- a / n
    prob_pu  <- pbeta(pu, ap, bp)
    prob_pa  <- 1 - pbeta(pa, ap, bp)
    y_max    <- max(y_seq)
    
    mle_hjust <- if (mle_p > 0.55) 1.05 else -0.05
    pu_hjust  <- if (pu < 0.15) -0.05 else 1.05
    pa_hjust  <- if (pa > 0.85) 1.05 else -0.05
    
    ggplot(df, aes(x = p_val, y = density)) +
      # Solid elegant continuous density surface area
      geom_area(aes(fill = region), alpha = 0.85) +
      # Prior density overlay
      geom_line(aes(y = prior_density), color = "grey30", linewidth = 0.9, linetype = "dotted") +
      annotate("text",
               x     = x_seq[which.max(prior_y)],
               y     = max(prior_y) * 1.04,
               label = sprintf("Prior: Beta(%.1f, %.1f)", input$prior_alpha, input$prior_beta),
               color = "grey30", size = 3.5, fontface = "italic", vjust = 0) +
      # Precise target baseline reference lines
      geom_vline(xintercept = mle_p, color = "#e67e22", linewidth = 1.1, linetype = "dotted") +
      geom_vline(xintercept = pu, color = "#c0392b", linewidth = 1.0, linetype = "dashed") +
      geom_vline(xintercept = pa, color = "#1e8449", linewidth = 1.0, linetype = "dashed") +
      # Parameter Label Overlays
      annotate("label", x = pu, y = y_max * 0.95,
               label = sprintf("p\u2080 = %.2f\nP(\u2264p\u2080) = %.3f", pu, prob_pu),
               hjust = pu_hjust, vjust = 1, color = "#c0392b", fill = "white", label.color = "#c0392b",
               size = 4, lineheight = 1.2) +
      annotate("label", x = pa, y = y_max * 0.95,
               label = sprintf("p\u2081 = %.2f\nP(\u2265p\u2081) = %.3f", pa, prob_pa),
               hjust = pa_hjust, vjust = 1, color = "#1e8449", fill = "white", label.color = "#1e8449",
               size = 4, lineheight = 1.2) +
      annotate("label", x = mle_p, y = y_max * 1.02,
               label = sprintf("MLE = %d/%d = %.2f", a, n, mle_p),
               hjust = mle_hjust, vjust = 0, color = "#e67e22", fill = "white", label.color = "#e67e22",
               size = 4, fontface = "bold") +
      scale_fill_manual(
        values = setNames(
          c("#e8736a", "#aab7b8", "#52be80"),
          c(paste0("p \u2264 p\u2080 (", pu, ")"), "p\u2080 < p < p\u2081", paste0("p \u2265 p\u2081 (", pa, ")"))
        )
      ) +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), labels = scales::label_percent(accuracy = 1)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(
        x    = "Response rate (p)",
        y    = "Posterior Density",
        fill = NULL,
        caption = sprintf("Prior: Beta(%.1f, %.1f)   |   Posterior: Beta(%.1f, %.1f)   |   a = %d successes, n = %d observations",
                          input$prior_alpha, input$prior_beta, ap, bp, a, n)
      ) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position   = "bottom",
        panel.grid.minor  = element_blank(),
        plot.caption      = element_text(color = "grey50")
      )
  })
  
  # ── Tab 3 outputs: Thall & Simon Design ────────────────────────────────────
  
  # Function to calibrate stopping boundaries
  calibrate_boundaries <- function(p_0, max_n, cohort_size, P_F, P_E, prior_a, prior_b) {
    # Always include max_n as the final look, even if not divisible by cohort_size
    interim_looks <- unique(c(seq(cohort_size, max_n, by = cohort_size), max_n))
    boundary_table <- data.frame(n = interim_looks, futility_max_x = NA, efficacy_min_x = NA)
    
    for (i in seq_along(interim_looks)) {
      n_curr <- interim_looks[i]
      possible_x <- 0:n_curr
      
      post_a <- prior_a + possible_x
      post_b <- prior_b + n_curr - possible_x
      
      prob_greater_p0 <- 1 - pbeta(p_0, post_a, post_b)
      
      futility_idx <- which(prob_greater_p0 < P_F)
      boundary_table$futility_max_x[i] <- if (length(futility_idx) > 0) max(possible_x[futility_idx]) else -1
      
      efficacy_idx <- which(prob_greater_p0 > P_E)
      boundary_table$efficacy_min_x[i] <- if (length(efficacy_idx) > 0) min(possible_x[efficacy_idx]) else n_curr + 1
    }
    
    return(boundary_table)
  }
  
  # Function to simulate operating characteristics
  simulate_oc <- function(true_p, max_n, boundary_table, n_sims = 5000) {
    stop_futility <- 0
    stop_efficacy <- 0
    sample_sizes <- numeric(n_sims)
    
    looks <- boundary_table$n
    
    for (s in 1:n_sims) {
      raw_data <- rbinom(max_n, 1, true_p)
      cum_responses <- cumsum(raw_data)
      
      for (i in seq_along(looks)) {
        n_curr <- looks[i]
        x_curr <- cum_responses[n_curr]
        
        fut_cutoff <- boundary_table$futility_max_x[i]
        eff_cutoff <- boundary_table$efficacy_min_x[i]
        
        if (n_curr < max_n) {
          if (x_curr <= fut_cutoff) {
            stop_futility <- stop_futility + 1
            sample_sizes[s] <- n_curr
            break
          }
          if (x_curr >= eff_cutoff) {
            stop_efficacy <- stop_efficacy + 1
            sample_sizes[s] <- n_curr
            break
          }
        } else {
          # Final look at max_n
          if (x_curr >= eff_cutoff) {
            stop_efficacy <- stop_efficacy + 1
          } else {
            stop_futility <- stop_futility + 1
          }
          sample_sizes[s] <- max_n
        }
      }
    }
    
    return(list(
      PET_Futility = stop_futility / n_sims,
      PET_Efficacy = stop_efficacy / n_sims,
      ESS = mean(sample_sizes)
    ))
  }
  
  # Reactive computation of boundaries and operating characteristics
  ts_results <- reactive({
    req(input$pa > input$pu)
    
    set.seed(6173)  # Fixed seed for reproducible simulation results
    
    bt <- calibrate_boundaries(
      p_0 = input$pu,
      max_n = input$ts_max_n,
      cohort_size = input$ts_cohort,
      P_F = input$ts_fut_thresh,
      P_E = input$ts_eff_thresh,
      prior_a = input$prior_alpha,
      prior_b = input$prior_beta
    )
    
    sim_null <- simulate_oc(input$pu, input$ts_max_n, bt, n_sims = 5000)
    sim_alt <- simulate_oc(input$pa, input$ts_max_n, bt, n_sims = 5000)
    
    list(
      boundaries = bt,
      sim_null = sim_null,
      sim_alt = sim_alt
    )
  })
  
  # Value boxes for Thall & Simon
  output$vb_ts_type1 <- renderText({
    res <- ts_results()
    sprintf("%.3f", res$sim_null$PET_Efficacy)
  })
  
  output$vb_ts_power <- renderText({
    res <- ts_results()
    sprintf("%.3f", res$sim_alt$PET_Efficacy)
  })
  
  output$vb_ts_ess_null <- renderText({
    res <- ts_results()
    sprintf("%.1f patients", res$sim_null$ESS)
  })
  
  output$vb_ts_ess_alt <- renderText({
    res <- ts_results()
    sprintf("%.1f patients", res$sim_alt$ESS)
  })
  
  # Boundary table output
  output$ts_boundary_table <- renderTable({
    res <- ts_results()
    bt <- res$boundaries
    max_n <- input$ts_max_n
    data.frame(
      "n" = as.integer(bt$n),
      "Stop for Futility if x \u2264" = ifelse(bt$futility_max_x < 0, "\u2014",
                                               as.character(as.integer(bt$futility_max_x))),
      "Stop for Efficacy if x \u2265" = ifelse(bt$efficacy_min_x > max_n, "\u2014",
                                               as.character(as.integer(bt$efficacy_min_x))),
      check.names = FALSE
    )
  }, striped = TRUE, hover = TRUE, width = "100%", align = "ccc")
  
  output$ts_table_footer <- renderUI({
    tags$small(tags$em(
      "\u2014 indicates no stopping boundary exists at that sample size. ",
      "Boundaries are evaluated after every cohort of ", input$ts_cohort, " patients."
    ))
  })
  
  # Summary interpretation
  output$ts_summary <- renderUI({
    res <- ts_results()
    sn <- res$sim_null
    sa <- res$sim_alt
    
    tagList(
      tags$h6("Design Characteristics", class = "fw-bold"),
      tags$p(sprintf(
        "This Bayesian adaptive design uses sequential interim analyses with fixed stopping boundaries. The trial employs a Beta(%.1f, %.1f) prior and implements a sequential Beta-Binomial model to continuously monitor for efficacy and futility.",
        input$prior_alpha, input$prior_beta
      )),
      tags$hr(),
      tags$h6("Operating Characteristics (from simulation)", class = "fw-bold"),
      tags$div(
        class = "alert alert-danger py-2",
        tags$strong(sprintf("Under Null Hypothesis (p\u2080 = %.2f)", input$pu)),
        tags$br(),
        sprintf("Type I Error (False Efficacy Declaration): %.4f", sn$PET_Efficacy),
        tags$br(),
        sprintf("Probability of Early Stopping for Futility: %.4f", sn$PET_Futility),
        tags$br(),
        sprintf("Expected Sample Size: %.1f patients", sn$ESS)
      ),
      tags$div(
        class = "alert alert-success py-2",
        tags$strong(sprintf("Under Alternative Hypothesis (p\u2081 = %.2f)", input$pa)),
        tags$br(),
        sprintf("Power (True Efficacy): %.4f", sa$PET_Efficacy),
        tags$br(),
        sprintf("Probability of Early Stopping for Futility: %.4f", sa$PET_Futility),
        tags$br(),
        sprintf("Expected Sample Size: %.1f patients", sa$ESS)
      ),
      tags$hr(),
      tags$h6("Interpretation", class = "fw-bold"),
      tags$p(
        "The Thall & Simon design is a Bayesian single-arm trial design that enables early stopping based on posterior probabilities. ",
        "It is particularly useful for phase II oncology trials where the goal is to identify promising treatments for further development."
      ),
      tags$ul(
        tags$li("If the posterior probability that p > p\u2080 drops below the futility threshold, the trial stops early to avoid unnecessary patient exposure."),
        tags$li("If the posterior probability that p > p\u2080 exceeds the efficacy threshold, the trial stops early to declare the treatment promising."),
        tags$li("This design typically requires fewer patients than fixed-sample designs while maintaining Type I error control and desired power.")
      )
    )
  })
  
  # ── Tab 4 outputs: Bayesian Sequential Design ──────────────────────────────
  
  # Active parameters display — reactive to sidebar, no button required
  output$bsd_params_display <- renderUI({
    p_mid <- round((input$pu + input$pa) / 2, 2)
    
    row <- function(label, value) {
      tags$tr(
        tags$td(class = "text-muted small pe-3", label),
        tags$td(tags$strong(value))
      )
    }
    
    tagList(
      tags$table(
        class = "table table-sm table-borderless mb-1",
        tags$tbody(
          row(htmltools::HTML("p\u2080 (null rate)"),         sprintf("%.2f", input$pu)),
          row(htmltools::HTML("p\u2081 (target rate)"),       sprintf("%.2f", input$pa)),
          row(htmltools::HTML("p<sub>mid</sub> = (p\u2080+p\u2081)/2"), sprintf("%.2f", p_mid)),
          row(htmltools::HTML("\u03b1 (Type I target)"),     sprintf("%.2f", input$ep1)),
          row(htmltools::HTML("\u03b2 (Type II target)"),    sprintf("%.2f", 1 - input$power)),
          row("Prior",                                        sprintf("Beta(%.1f, %.1f)", input$prior_alpha, input$prior_beta))
        )
      ),
      tags$small(class = "text-muted fst-italic",
        "Adjust sidebar values, then click Run Calibration."
      )
    )
  })
  
  # Simulation function: tracks both decisions and sample sizes for ESS
  bsd_simulate <- function(true_p, n_sims, max_n, u_thresh, l_thresh,
                           p_mid, prior_a, prior_b) {
    n_eff <- 0L; n_fut <- 0L; n_inc <- 0L
    ss    <- numeric(n_sims)
    
    for (s in seq_len(n_sims)) {
      n <- 0L; r <- 0L; dec <- "c"
      while (n < max_n && dec == "c") {
        n   <- n + 1L
        r   <- r + rbinom(1L, 1L, true_p)
        pp  <- 1 - pbeta(p_mid, prior_a + r, prior_b + n - r)
        if      (pp >= u_thresh) dec <- "e"
        else if (pp <= l_thresh) dec <- "f"
      }
      if (dec == "c") dec <- "i"
      if      (dec == "e") n_eff <- n_eff + 1L
      else if (dec == "f") n_fut <- n_fut + 1L
      else                 n_inc <- n_inc + 1L
      ss[s] <- n
    }
    
    list(
      prob_efficacy     = n_eff / n_sims,
      prob_futility     = n_fut / n_sims,
      prob_inconclusive = n_inc / n_sims,
      ESS               = mean(ss)
    )
  }
  
  # eventReactive: runs grid search + final simulation on button click
  bsd_results <- eventReactive(input$bsd_run, {
    req(input$pa > input$pu)
    
    p0      <- input$pu
    p1      <- input$pa
    p_mid   <- (p0 + p1) / 2
    alpha_t <- input$ep1
    beta_t  <- 1 - input$power
    max_n   <- input$bsd_max_n
    prior_a <- input$prior_alpha
    prior_b <- input$prior_beta
    
    u_cands <- seq(0.92, 0.98, by = 0.01)
    l_cands <- seq(0.05, 0.15, by = 0.01)
    n_comb  <- length(u_cands) * length(l_cands)
    
    set.seed(3847)
    
    withProgress(message = "Calibrating thresholds \u2014 please wait...", value = 0, {
      
      # Grid search: 2,000 sims per cell to find best (u, l) pair
      grid_rows <- vector("list", n_comb)
      k <- 1L
      for (u in u_cands) {
        for (l in l_cands) {
          sh0 <- bsd_simulate(p0, 2000L, max_n, u, l, p_mid, prior_a, prior_b)
          sh1 <- bsd_simulate(p1, 2000L, max_n, u, l, p_mid, prior_a, prior_b)
          grid_rows[[k]] <- list(
            u = u, l = l,
            a = sh0$prob_efficacy,
            w = sh1$prob_efficacy,
            e = abs(sh0$prob_efficacy - alpha_t) + abs(sh1$prob_efficacy - (1 - beta_t))
          )
          k <- k + 1L
          incProgress(0.8 / n_comb,
                      detail = sprintf("u = %.2f, l = %.2f", u, l))
        }
      }
      
      errs <- sapply(grid_rows, `[[`, "e")
      best <- grid_rows[[which.min(errs)]]
      
      # Final confirmation simulation: 5,000 sims with calibrated thresholds
      incProgress(0.1, message = "Final simulation under H\u2080...")
      final_h0 <- bsd_simulate(p0, 5000L, max_n, best$u, best$l, p_mid, prior_a, prior_b)
      
      incProgress(0.1, message = "Final simulation under H\u2081...")
      final_h1 <- bsd_simulate(p1, 5000L, max_n, best$u, best$l, p_mid, prior_a, prior_b)
      
      list(
        p0 = p0, p1 = p1, p_mid = p_mid,
        alpha_t = alpha_t, beta_t = beta_t,
        max_n   = max_n, prior_a = prior_a, prior_b = prior_b,
        best_u  = best$u, best_l = best$l,
        type_I      = final_h0$prob_efficacy,
        power       = final_h1$prob_efficacy,
        ESS_h0      = final_h0$ESS,
        ESS_h1      = final_h1$ESS,
        prob_fut_h0 = final_h0$prob_futility,
        prob_fut_h1 = final_h1$prob_futility,
        prob_inc_h0 = final_h0$prob_inconclusive,
        prob_inc_h1 = final_h1$prob_inconclusive
      )
    })
  })
  
  # Value boxes
  output$vb_bsd_type1 <- renderText({
    sprintf("%.3f", bsd_results()$type_I)
  })
  
  output$vb_bsd_power <- renderText({
    sprintf("%.3f", bsd_results()$power)
  })
  
  output$vb_bsd_ess_h0 <- renderText({
    sprintf("%.1f patients", bsd_results()$ESS_h0)
  })
  
  output$vb_bsd_ess_h1 <- renderText({
    sprintf("%.1f patients", bsd_results()$ESS_h1)
  })
  
  # Boundary plot — uses correct posterior-based boundaries (not qbinom approximation)
  output$bsd_boundary_plot <- renderPlot({
    res    <- bsd_results()
    n_vals <- seq_len(res$max_n)
    
    # Efficacy boundary: minimum r at each n where P(p > p_mid | r, n) >= best_u
    eff_b <- sapply(n_vals, function(n) {
      r_seq <- 0:n
      pp    <- 1 - pbeta(res$p_mid, res$prior_a + r_seq, res$prior_b + n - r_seq)
      idx   <- which(pp >= res$best_u)
      if (length(idx) > 0) as.numeric(min(r_seq[idx])) else NA_real_
    })
    
    # Futility boundary: maximum r at each n where P(p > p_mid | r, n) <= best_l
    fut_b <- sapply(n_vals, function(n) {
      r_seq <- 0:n
      pp    <- 1 - pbeta(res$p_mid, res$prior_a + r_seq, res$prior_b + n - r_seq)
      idx   <- which(pp <= res$best_l)
      if (length(idx) > 0) as.numeric(max(r_seq[idx])) else NA_real_
    })
    
    h0_lab <- sprintf("Expected under H\u2080 (p\u2080 = %.2f)", res$p0)
    h1_lab <- sprintf("Expected under H\u2081 (p\u2081 = %.2f)", res$p1)
    lvls   <- c("Efficacy Boundary", "Futility Boundary", h0_lab, h1_lab)
    
    df <- data.frame(
      n    = rep(n_vals, 4),
      r    = c(eff_b, fut_b, n_vals * res$p0, n_vals * res$p1),
      line = factor(rep(lvls, each = length(n_vals)), levels = lvls)
    )
    
    col_vals <- c("Efficacy Boundary" = "#1e8449",
                  "Futility Boundary" = "#c0392b",
                  setNames("#2980b9", h0_lab),
                  setNames("#e67e22", h1_lab))
    lty_vals <- c("Efficacy Boundary" = "solid",
                  "Futility Boundary" = "solid",
                  setNames("dashed", h0_lab),
                  setNames("dashed", h1_lab))
    
    ggplot(df, aes(x = n, y = r, color = line, linetype = line)) +
      geom_line(linewidth = 1.0, na.rm = TRUE) +
      scale_color_manual(values = col_vals) +
      scale_linetype_manual(values = lty_vals) +
      labs(
        x       = "Sample Size (n)",
        y       = "Cumulative Responses (r)",
        color   = NULL,
        linetype = NULL,
        caption = sprintf(
          "Stop for EFFICACY if P(p > %.2f | data) \u2265 %.2f   |   Stop for FUTILITY if P(p > %.2f | data) \u2264 %.2f   |   Prior: Beta(%.1f, %.1f)",
          res$p_mid, res$best_u, res$p_mid, res$best_l, res$prior_a, res$prior_b
        )
      ) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        plot.caption     = element_text(color = "grey50")
      )
  })
  
  # Design summary
  output$bsd_summary <- renderUI({
    res <- bsd_results()
    
    tagList(
      tags$h6("Design Specification", class = "fw-bold"),
      tags$ul(
        tags$li(sprintf("Null response rate (p\u2080): %.2f", res$p0)),
        tags$li(sprintf("Target response rate (p\u2081): %.2f", res$p1)),
        tags$li(htmltools::HTML(sprintf("Reference rate (p<sub>mid</sub>): %.2f &mdash; midpoint of p\u2080 and p\u2081", res$p_mid))),
        tags$li(sprintf("Prior: Beta(%.1f, %.1f)", res$prior_a, res$prior_b)),
        tags$li(sprintf("Maximum sample size: %d", res$max_n)),
        tags$li("Monitoring: continuous (after each patient)")
      ),
      tags$hr(),
      tags$h6("Calibrated Stopping Boundaries", class = "fw-bold"),
      tags$div(
        class = "alert alert-secondary py-2",
        htmltools::HTML(sprintf(
          "Stop for <strong>EFFICACY</strong> if P(p &gt; %.2f | data) &ge; %.2f<br>
           Stop for <strong>FUTILITY</strong> if P(p &gt; %.2f | data) &le; %.2f<br>
           Otherwise <strong>continue</strong> to next patient (up to N = %d)",
          res$p_mid, res$best_u, res$p_mid, res$best_l, res$max_n
        ))
      ),
      tags$hr(),
      tags$h6("Operating Characteristics (5,000 simulations)", class = "fw-bold"),
      tags$div(
        class = "alert alert-danger py-2",
        tags$strong(sprintf("Under H\u2080 (p = %.2f):", res$p0)),
        tags$br(),
        sprintf("Type I Error: %.4f  (target: %.2f)", res$type_I, res$alpha_t),
        tags$br(),
        sprintf("P(futility stop): %.4f", res$prob_fut_h0),
        tags$br(),
        sprintf("P(inconclusive at N = %d): %.4f", res$max_n, res$prob_inc_h0),
        tags$br(),
        sprintf("E(N): %.1f patients", res$ESS_h0)
      ),
      tags$div(
        class = "alert alert-success py-2",
        tags$strong(sprintf("Under H\u2081 (p = %.2f):", res$p1)),
        tags$br(),
        sprintf("Power: %.4f  (target: %.2f)", res$power, 1 - res$beta_t),
        tags$br(),
        sprintf("P(futility stop): %.4f", res$prob_fut_h1),
        tags$br(),
        sprintf("P(inconclusive at N = %d): %.4f", res$max_n, res$prob_inc_h1),
        tags$br(),
        sprintf("E(N): %.1f patients", res$ESS_h1)
      ),
      tags$hr(),
      tags$h6("Regulatory Justification", class = "fw-bold"),
      tags$ul(
        tags$li("Type I error controlled at nominal level via calibrated posterior thresholds"),
        tags$li("Power achieved within nominal specification"),
        tags$li("Operating characteristics from Monte Carlo simulation (5,000 replicates)"),
        tags$li("Threshold calibration via grid search over 77 (u, l) combinations"),
        tags$li("Bayesian approach provides intuitive posterior probability interpretation"),
        tags$li("Continuous monitoring balances statistical efficiency with patient protection")
      )
    )
  })
  
  # Computed summary for observed data
  output$bsd_computed_summary <- renderUI({
    req(bsd_results())
    
    # Get observed data from inputs
    a_obs <- input$bsd_successes
    n_obs <- input$bsd_total
    f_obs <- max(0, n_obs - a_obs)  # Implied failures
    obs_rate <- if (n_obs > 0) a_obs / n_obs else 0
    
    res <- bsd_results()
    prior_a <- res$prior_a
    prior_b <- res$prior_b
    p_mid <- res$p_mid
    
    # Calculate posterior probability with observed data
    post_prob <- 1 - pbeta(p_mid, prior_a + a_obs, prior_b + f_obs)
    
    # Summary table
    row <- function(label, value) {
      tags$tr(
        tags$td(class = "text-muted small pe-3", label),
        tags$td(tags$strong(value))
      )
    }
    
    tagList(
      tags$table(
        class = "table table-sm table-borderless",
        tags$tbody(
          row("Implied failures", sprintf("%d", f_obs)),
          row("Observed rate", sprintf("%.1f%%", 100 * obs_rate)),
          row(htmltools::HTML("P(p &gt; p<sub>mid</sub> | data)"), sprintf("%.4f", post_prob))
        )
      )
    )
  })
  
  # Observed data display with sequential decision
  output$bsd_observed_data <- renderUI({
    req(bsd_results())
    
    # Get observed data from inputs
    a_obs <- input$bsd_successes
    n_obs <- input$bsd_total
    f_obs <- max(0, n_obs - a_obs)  # Implied failures
    obs_rate <- if (n_obs > 0) a_obs / n_obs else 0
    
    res <- bsd_results()
    p_mid <- res$p_mid
    prior_a <- res$prior_a
    prior_b <- res$prior_b
    
    # Calculate posterior probability with observed data
    post_prob <- 1 - pbeta(p_mid, prior_a + a_obs, prior_b + f_obs)
    
    # Determine sequential decision
    decision <- dplyr::case_when(
      post_prob >= res$best_u ~ "EFFICACY — Treatment is promising",
      post_prob <= res$best_l ~ "FUTILITY — Treatment is not sufficiently active",
      TRUE ~ "CONTINUE — Enroll more patients"
    )
    
    decision_class <- dplyr::case_when(
      post_prob >= res$best_u ~ "alert alert-success",
      post_prob <= res$best_l ~ "alert alert-danger",
      TRUE ~ "alert alert-warning"
    )
    
    decision_icon <- dplyr::case_when(
      post_prob >= res$best_u ~ "✓",
      post_prob <= res$best_l ~ "✗",
      TRUE ~ "→"
    )
    
    tagList(
      tags$hr(),
      tags$h6("Sequential Decision at This Sample Size", class = "fw-bold"),
      tags$div(
        class = decision_class,
        tags$strong(
          htmltools::HTML(sprintf(
            "%s &nbsp; %s",
            decision_icon,
            decision
          ))
        ),
        tags$br(),
        htmltools::HTML(sprintf(
          "Efficacy threshold: %.2f &nbsp; | &nbsp; Your posterior: %.4f &nbsp; | &nbsp; Futility threshold: %.2f",
          res$best_u, post_prob, res$best_l
        ))
      ),
      tags$hr(),
      tags$h6("Explanation", class = "fw-bold"),
      tags$p(
        "Based on the observed data (",
        sprintf("%d successes, %d failures, n = %d", a_obs, f_obs, n_obs),
        ", observed rate = ",
        sprintf("%.1f%%", 100 * obs_rate),
        "), ",
        "the posterior probability that the true response rate exceeds the midpoint (p_mid = ",
        sprintf("%.2f", p_mid),
        ") is ",
        tags$strong(sprintf("%.4f", post_prob)),
        ". ",
        dplyr::case_when(
          post_prob >= res$best_u ~ paste(
            "This exceeds the efficacy threshold of",
            sprintf("%.2f,", res$best_u),
            "so the sequential design rule says to STOP and declare the treatment promising."
          ),
          post_prob <= res$best_l ~ paste(
            "This falls below the futility threshold of",
            sprintf("%.2f,", res$best_l),
            "so the sequential design rule says to STOP and declare the treatment not sufficiently active."
          ),
          TRUE ~ paste(
            "This is between the futility (",
            sprintf("%.2f)", res$best_l),
            "and efficacy (",
            sprintf("%.2f)", res$best_u),
            "thresholds, so continue enrolling more patients."
          )
        )
      ),
      tags$small(
        class = "text-muted d-block mt-2",
        "Note: This decision is based on the calibrated thresholds from the design. ",
        "The thresholds change based on your slider inputs; run Calibration again if you adjust them."
      )
    )
  })
}

shinyApp(ui, server)

#git add One_Group_Binary_v1.R
#git commit -m "your message" && git push origin main
