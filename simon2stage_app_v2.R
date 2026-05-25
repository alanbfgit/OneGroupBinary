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
      geom_line(aes(y = prior_density), color = "grey30", linewidth = 0.9, linetype = "dashed") +
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
}

shinyApp(ui, server)