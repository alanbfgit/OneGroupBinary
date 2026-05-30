library(shiny)
library(bslib)
library(ggplot2)
library(scales)

# ── Helper functions ──────────────────────────────────────────────────────────

# LR martingale E-value: computed on log scale for numerical stability
lr_e_value <- function(a, b, p0, p1) {
  if ((a + b) == 0L) return(1)
  exp(a * log(p1 / p0) + b * log((1 - p1) / (1 - p0)))
}

# One-sided exact binomial p-value: P(X >= a | n, p0)
exact_pval <- function(a, n, p0) {
  if (n == 0L) return(NA_real_)
  binom.test(a, n, p0, alternative = "greater")$p.value
}

# Clopper-Pearson exact 95% CI
exact_ci <- function(a, n) {
  if (n == 0L) return(c(NA_real_, NA_real_))
  as.numeric(binom.test(a, n)$conf.int)
}

make_decision <- function(e, n, eff_thresh, fut_thresh, max_n) {
  if (e >= eff_thresh) return("EFFICACY")
  if (e <= fut_thresh) return("FUTILITY")
  if (n >= max_n)      return("INCONCLUSIVE")
  return("CONTINUE")
}

fmt_pval <- function(p) {
  if (is.na(p)) return("\u2014")
  if (p < 0.001) "< 0.001" else sprintf("%.4f", p)
}

# Shared decision alert UI (used by both tabs)
decision_alert <- function(dec, e, n, a, b, eff_thresh, fut_thresh,
                           first_eff = NA_integer_, first_fut = NA_integer_) {
  cls <- switch(dec,
    EFFICACY     = "alert alert-success",
    FUTILITY     = "alert alert-danger",
    INCONCLUSIVE = "alert alert-warning",
    CONTINUE     = "alert alert-info"
  )
  sym <- switch(dec,
    EFFICACY     = "\u2713",
    FUTILITY     = "\u2717",
    INCONCLUSIVE = "\u26a0",
    CONTINUE     = "\u2192"
  )
  expl <- switch(dec,
    EFFICACY     = sprintf("E \u2265 %.2f (= 1/\u03b1): reject H\u2080. Anytime-valid evidence of efficacy.", eff_thresh),
    FUTILITY     = sprintf("E \u2264 %.4f (= \u03b1): stop for futility. Insufficient evidence of activity.", fut_thresh),
    INCONCLUSIVE = "Maximum N reached without crossing either boundary. Trial is inconclusive.",
    CONTINUE     = sprintf("E is between %.4f (futility) and %.2f (efficacy). Continue enrolling.", fut_thresh, eff_thresh)
  )
  cross_note <- ""
  if (!is.na(first_eff))
    cross_note <- sprintf(" Efficacy boundary first crossed at patient %d.", first_eff)
  else if (!is.na(first_fut))
    cross_note <- sprintf(" Futility boundary first crossed at patient %d.", first_fut)

  tagList(
    tags$div(
      class = cls,
      tags$strong(paste(sym, dec)),
      tags$br(),
      sprintf("n = %d  |  a = %d  |  b = %d  |  E = %.4f", n, a, b, e),
      if (nchar(cross_note) > 0) tags$span(tags$br(), cross_note) else NULL,
      tags$br(),
      tags$small(expl)
    )
  )
}

# Shared results summary table
summary_table_df <- function(a, b, n, rate, ci, e, pval, eff_thresh, fut_thresh, dec) {
  data.frame(
    Parameter = c(
      "Successes (a)", "Failures (b)", "Total (n)",
      "Observed rate", "95% CI (Clopper\u2013Pearson)",
      "E-value", "1/E  (min \u03b1 to reject)",
      "One-sided p-value",
      "Efficacy threshold  (1/\u03b1)", "Futility threshold  (\u03b1)",
      "Decision"
    ),
    Value = c(
      as.character(a),
      as.character(b),
      as.character(n),
      if (is.na(rate)) "\u2014" else sprintf("%.1f%%", 100 * rate),
      if (any(is.na(ci))) "\u2014" else sprintf("(%.3f,  %.3f)", ci[1], ci[2]),
      sprintf("%.4f", e),
      sprintf("%.4f", 1 / e),
      fmt_pval(pval),
      sprintf("%.2f", eff_thresh),
      sprintf("%.4f", fut_thresh),
      dec
    ),
    check.names = FALSE
  )
}

# ── Helper: observed rate threshold that corresponds to an E-value threshold ──
# For a given n and E threshold, find the observed success count a such that
# E(a, n-a, p0, p1) = E_thresh. Returns the observed rate a/n.
# Used to make futility criterion interpretable to users.

rate_for_e_threshold <- function(e_thresh, n, p0, p1) {
  log_r1 <- log(p1 / p0)
  log_r0 <- log((1 - p1) / (1 - p0))
  
  # E = exp(a * log_r1 + (n-a) * log_r0)
  # log(E) = a * log_r1 + (n-a) * log_r0
  # log(E) = a * (log_r1 - log_r0) + n * log_r0
  # a = (log(E) - n * log_r0) / (log_r1 - log_r0)
  
  if (abs(log_r1 - log_r0) < 1e-10) return(NA_real_)  # Degenerate case
  
  a_exact <- (log(e_thresh) - n * log_r0) / (log_r1 - log_r0)
  a_exact <- max(0, min(n, a_exact))
  a_exact / n
}

# ── OC simulation (used by Tab 3) ────────────────────────────────────────────
# Vectorized LR martingale simulator. Returns prob_efficacy, prob_futility,
# prob_inconclusive, and ESS for a given true p under the design parameters.

simulate_e_oc <- function(true_p, p0, p1, alpha, max_n, n_sims = 5000) {
  eff_thresh <- 1 / alpha
  fut_thresh <- alpha
  log_r1     <- log(p1 / p0)
  log_r0     <- log((1 - p1) / (1 - p0))

  obs <- matrix(rbinom(n_sims * max_n, 1L, true_p), nrow = n_sims, ncol = max_n)

  # Row-wise cumulative log E-values via column scan (avoids slow apply)
  log_e <- obs * log_r1 + (1L - obs) * log_r0
  for (j in 2:max_n) log_e[, j] <- log_e[, j - 1L] + log_e[, j]
  e_mat <- exp(log_e)

  # First column where each row crosses a boundary
  cross_first <- function(bool_mat) {
    result <- rep(max_n + 1L, nrow(bool_mat))
    for (j in seq_len(ncol(bool_mat))) {
      pending <- result > max_n
      if (!any(pending)) break
      result[pending & bool_mat[, j]] <- j
    }
    result
  }

  first_eff <- cross_first(e_mat >= eff_thresh)
  first_fut <- cross_first(e_mat <= fut_thresh)

  eff_win <- first_eff < first_fut
  fut_win <- first_fut < first_eff
  stop_n  <- pmin(first_eff, first_fut)
  stop_n[stop_n > max_n] <- max_n

  list(
    prob_efficacy     = mean(eff_win),
    prob_futility     = mean(fut_win),
    prob_inconclusive = mean(!eff_win & !fut_win),
    ess               = mean(stop_n)
  )
}

# ── UI ────────────────────────────────────────────────────────────────────────

ui <- page_sidebar(
  title = "Always-Safe Sequential Test \u2014 LR Martingale (E-Value)",
  theme = bs_theme(version = 5, bootswatch = "flatly"),

  sidebar = sidebar(
    width = 310,
    tags$p(tags$strong("Design Parameters")),
    sliderInput("p0",    htmltools::HTML("Null rate (p\u2080)"),
                min = 0.05, max = 0.60, value = 0.20, step = 0.01),
    sliderInput("p1",    htmltools::HTML("Target rate (p\u2081)"),
                min = 0.10, max = 0.80, value = 0.40, step = 0.01),
    sliderInput("alpha", htmltools::HTML("Type I error (\u03b1)"),
                min = 0.01, max = 0.20, value = 0.05, step = 0.01),
    numericInput("max_n", "Maximum N", value = 50, min = 10, max = 1000, step = 1),
    hr(),
    tags$p(tags$strong("Futility Boundary")),
    numericInput("fut_e", htmltools::HTML("Futility threshold (E \u2264)"),
                 value = 0.05, min = 0.001, max = 0.50, step = 0.001),
    tags$small(tags$em(
      "E-value below which to stop for futility. ",
      "See 'Equivalent Observed Rate' box below."
    )),
    hr(),
    uiOutput("sidebar_boundaries"),
    hr(),
    uiOutput("param_warning")
  ),

  navset_card_tab(

    # ── Tab 1: Snapshot ────────────────────────────────────────────────────────
    nav_panel(
      "Snapshot",
      layout_column_wrap(
        width = 1/2, fill = FALSE,
        card(
          card_header("Observed Data"),
          card_body(
            numericInput("snap_a", "Successes (a)", value = 12, min = 0, max = 10000, step = 1),
            numericInput("snap_b", "Failures (b)",  value = 18, min = 0, max = 10000, step = 1)
          )
        ),
        card(
          card_header("Decision"),
          card_body(uiOutput("snap_decision"))
        )
      ),
      layout_column_wrap(
        width = 1/4, fill = FALSE,
        value_box(
          title    = "Observed Rate",
          value    = textOutput("snap_rate"),
          showcase = bsicons::bs_icon("percent"),
          theme    = "primary"
        ),
        value_box(
          title    = htmltools::HTML("95% CI (Clopper&ndash;Pearson)"),
          value    = textOutput("snap_ci"),
          showcase = bsicons::bs_icon("bar-chart"),
          theme    = "info"
        ),
        value_box(
          title    = "E-Value",
          value    = textOutput("snap_e"),
          showcase = bsicons::bs_icon("lightning-charge"),
          theme    = "warning"
        ),
        value_box(
          title    = htmltools::HTML("1/E &nbsp;(min \u03b1)"),
          value    = textOutput("snap_inv_e"),
          showcase = bsicons::bs_icon("shield-check"),
          theme    = "secondary"
        ),
        value_box(
          title    = "One-Sided p-value",
          value    = textOutput("snap_pval"),
          showcase = bsicons::bs_icon("calculator"),
          theme    = "secondary"
        )
      ),
      card(
        card_header("Results Summary"),
        tableOutput("snap_table"),
        card_footer(uiOutput("snap_footer"))
      ),
      accordion(
        id = "snap_accordion",
        accordion_panel(
          htmltools::HTML("One-Sided p-Value Details &nbsp; <span style='font-size: 0.85rem; color: #7f8c8d;'>(click to expand)</span>"),
          uiOutput("snap_pval_detail")
        )
      )
    ),

    # ── Tab 2: Sequential ──────────────────────────────────────────────────────
    nav_panel(
      "Sequential",
      layout_column_wrap(
        width = 1/2, fill = FALSE,
        card(
          card_header("Data Input"),
          card_body(
            radioButtons("seq_method", "Input method:",
                         choices  = c("Paste sequence" = "text", "Upload CSV" = "csv"),
                         selected = "text", inline = TRUE),
            hr(),
            conditionalPanel(
              "input.seq_method == 'text'",
              tags$p(tags$small(tags$em(
                "Comma-separated 1s and 0s in patient order (1 = success, 0 = failure):"
              ))),
              textAreaInput("seq_text", NULL, value = "",
                            placeholder = "e.g. 1,0,1,1,0,1,0,0,1,...",
                            rows = 4, width = "100%"),
              actionButton("seq_run", "Process Sequence",
                           class = "btn-primary w-100", icon = icon("play"))
            ),
            conditionalPanel(
              "input.seq_method == 'csv'",
              tags$p(tags$small(tags$em(
                "Single column of 1s and 0s, one row per patient.",
                "A header row (non-numeric) is ignored automatically."
              ))),
              fileInput("seq_file", NULL, accept = ".csv",
                        buttonLabel = "Browse...", placeholder = "No file selected"),
              actionButton("seq_run_csv", "Process File",
                           class = "btn-primary w-100", icon = icon("play"))
            )
          )
        ),
        card(
          card_header("Decision at Final Observation"),
          card_body(uiOutput("seq_decision"))
        )
      ),
      layout_column_wrap(
        width = 1/4, fill = FALSE,
        value_box(
          title    = "Final Observed Rate",
          value    = textOutput("seq_rate"),
          showcase = bsicons::bs_icon("percent"),
          theme    = "primary"
        ),
        value_box(
          title    = htmltools::HTML("95% CI (Clopper&ndash;Pearson)"),
          value    = textOutput("seq_ci"),
          showcase = bsicons::bs_icon("bar-chart"),
          theme    = "info"
        ),
        value_box(
          title    = "Final E-Value",
          value    = textOutput("seq_e"),
          showcase = bsicons::bs_icon("lightning-charge"),
          theme    = "warning"
        ),
        value_box(
          title    = htmltools::HTML("1/E &nbsp;(min \u03b1)"),
          value    = textOutput("seq_inv_e"),
          showcase = bsicons::bs_icon("shield-check"),
          theme    = "secondary"
        ),
        value_box(
          title    = "One-Sided p-value",
          value    = textOutput("seq_pval"),
          showcase = bsicons::bs_icon("calculator"),
          theme    = "secondary"
        )
      ),
      card(
        full_screen = TRUE,
        card_header("E-Value Trajectory"),
        plotOutput("seq_plot", height = "450px"),
        card_footer(uiOutput("seq_plot_footer"))
      ),
      card(
        card_header("Results Summary"),
        tableOutput("seq_table"),
        card_footer(uiOutput("seq_footer"))
      ),
      accordion(
        id = "seq_accordion",
        accordion_panel(
          htmltools::HTML("One-Sided p-Value Details &nbsp; <span style='font-size: 0.85rem; color: #7f8c8d;'>(click to expand)</span>"),
          uiOutput("seq_pval_detail")
        )
      )
    ),

    # ── Tab 3: Operating Characteristics ──────────────────────────────────────
    nav_panel(
      "Operating Characteristics",

      layout_column_wrap(
        width = 1/3, fill = FALSE,
        card(
          card_header("Simulation Settings"),
          card_body(
            sliderInput("oc_sims", "Simulations per p value",
                        min = 1000, max = 20000, value = 5000, step = 1000),
            tags$small(tags$em(
              "5,000 is sufficient for smooth display. ",
              "Increase for publication-quality precision."
            )),
            tags$br(), tags$br(),
            actionButton("oc_run", "Run OC Simulation",
                         class = "btn-primary w-100", icon = icon("play"))
          )
        ),
        card(
          card_header(htmltools::HTML("At p\u2080 (under H\u2080)")),
          card_body(uiOutput("oc_h0_summary"))
        ),
        card(
          card_header(htmltools::HTML("At p\u2081 (under H\u2081)")),
          card_body(uiOutput("oc_h1_summary"))
        )
      ),

      layout_column_wrap(
        width = 1/2,
        card(
          full_screen = TRUE,
          card_header("Decision Probabilities vs. True Response Rate"),
          plotOutput("oc_plot_probs"),
          card_footer(uiOutput("oc_plot_footer")),
          style = "min-height: 420px;"
        ),
        card(
          full_screen = TRUE,
          card_header("Expected Sample Size vs. True Response Rate"),
          plotOutput("oc_plot_ess"),
          style = "min-height: 420px;"
        )
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # ── Sidebar ────────────────────────────────────────────────────────────────

  output$sidebar_boundaries <- renderUI({
    eff <- 1 / input$alpha
    fut <- input$fut_e
    
    # Compute equivalent observed rates at key sample sizes
    rates_at_n <- data.frame(
      n = c(10, 20, 30, input$max_n),
      rate = NA_real_
    )
    for (i in seq_len(nrow(rates_at_n))) {
      rates_at_n$rate[i] <- rate_for_e_threshold(fut, rates_at_n$n[i], input$p0, input$p1)
    }
    rates_at_n <- rates_at_n[!is.na(rates_at_n$rate), ]
    
    tagList(
      tags$p(tags$strong("Stopping Boundaries")),
      tags$table(
        class = "table table-sm table-borderless mb-2",
        tags$tbody(
          tags$tr(
            tags$td(class = "text-muted small pe-3",
                    htmltools::HTML("Efficacy &nbsp;(E \u2265)")),
            tags$td(tags$strong(sprintf("%.2f", eff)))
          ),
          tags$tr(
            tags$td(class = "text-muted small pe-3",
                    htmltools::HTML("Futility &nbsp;(E \u2264)")),
            tags$td(tags$strong(sprintf("%.4f", fut)))
          )
        )
      ),
      tags$p(tags$strong("Equivalent Observed Rate (Futility)"), style = "font-size: 0.85rem; margin-top: 0.5rem;"),
      tags$table(
        class = "table table-sm table-borderless mb-0",
        style = "font-size: 0.8rem;",
        tags$tbody(
          lapply(seq_len(nrow(rates_at_n)), function(i) {
            tags$tr(
              tags$td(class = "text-muted pe-2", sprintf("At n = %d:", rates_at_n$n[i])),
              tags$td(tags$strong(sprintf("%.1f%%", 100 * rates_at_n$rate[i])))
            )
          })
        )
      )
    )
  })

  output$param_warning <- renderUI({
    if (input$p1 <= input$p0)
      tags$div(class = "alert alert-danger mb-0",
        tags$strong("Invalid:"),
        htmltools::HTML(" p\u2081 must exceed p\u2080.")
      )
  })

  # ── Tab 1: Snapshot ────────────────────────────────────────────────────────

  snap <- reactive({
    req(input$p1 > input$p0)
    a  <- max(0L, as.integer(input$snap_a))
    b  <- max(0L, as.integer(input$snap_b))
    n  <- a + b
    p0 <- input$p0;  p1 <- input$p1;  al <- input$alpha

    e          <- lr_e_value(a, b, p0, p1)
    rate       <- if (n > 0L) a / n else NA_real_
    pval       <- exact_pval(a, n, p0)
    ci         <- exact_ci(a, n)
    eff_thresh <- 1 / al
    fut_thresh <- input$fut_e
    dec        <- make_decision(e, n, eff_thresh, fut_thresh, input$max_n)

    list(a = a, b = b, n = n, rate = rate, e = e, pval = pval, ci = ci,
         dec = dec, eff_thresh = eff_thresh, fut_thresh = fut_thresh)
  })

  output$snap_rate  <- renderText({
    s <- snap()
    if (is.na(s$rate)) "\u2014" else sprintf("%.1f%%", 100 * s$rate)
  })
  output$snap_ci    <- renderText({
    s <- snap()
    if (any(is.na(s$ci))) "\u2014" else sprintf("(%.3f, %.3f)", s$ci[1], s$ci[2])
  })
  output$snap_e     <- renderText({ sprintf("%.4f", snap()$e) })
  output$snap_inv_e <- renderText({ sprintf("%.4f", 1 / snap()$e) })
  output$snap_pval  <- renderText({ fmt_pval(snap()$pval) })

  output$snap_decision <- renderUI({
    s <- snap()
    decision_alert(s$dec, s$e, s$n, s$a, s$b, s$eff_thresh, s$fut_thresh)
  })

  output$snap_table <- renderTable({
    s <- snap()
    summary_table_df(s$a, s$b, s$n, s$rate, s$ci, s$e,
                     s$pval, s$eff_thresh, s$fut_thresh, s$dec)
  }, striped = TRUE, hover = TRUE, width = "100%")

  output$snap_footer <- renderUI({
    tags$small(tags$em(htmltools::HTML(
      "E-value: LR martingale = (p\u2081/p\u2080)<sup>a</sup>",
      " &times; ((1&minus;p\u2081)/(1&minus;p\u2080))<sup>b</sup>. &ensp;",
      "95% CI: Clopper&ndash;Pearson exact. &ensp;",
      "1/E is the minimum \u03b1 at which you would reject H\u2080 with the current data."
    )))
  })

  output$snap_pval_detail <- renderUI({
    s <- snap()
    tagList(
      tags$h6("What is this?", class = "fw-bold mt-3 mb-2"),
      tags$p(
        "The one-sided p-value is the probability of observing ",
        tags$strong(sprintf("%d or more", s$a)),
        " successes if the true response rate were exactly ",
        tags$strong(sprintf("p\u2080 = %.2f%%", 100 * input$p0)),
        ". It is a classical hypothesis test (exact binomial test) with:",
        tags$br(),
        "Null hypothesis: p = p\u2080", tags$br(),
        "Alternative hypothesis: p > p\u2080"
      ),
      tags$h6("Interpretation", class = "fw-bold mt-3 mb-2"),
      tags$ul(
        tags$li(
          tags$strong("Very small p-value (e.g., < 0.001):"),
          " Strong evidence against the null. Your data are unlikely if the null were true."
        ),
        tags$li(
          tags$strong("Moderate p-value (e.g., 0.01\u20130.05):"),
          " Suggestive evidence against the null, consistent with the E-value approach."
        ),
        tags$li(
          tags$strong("Large p-value (e.g., > 0.10):"),
          " Weak evidence against the null. The data are consistent with p = p\u2080."
        )
      ),
      tags$h6("How does it relate to the E-value?", class = "fw-bold mt-3 mb-2"),
      tags$p(
        "Both the p-value and E-value measure evidence against H\u2080, but differently:",
        tags$br(),
        tags$strong("E-value:"),
        " A likelihood ratio that can be computed sequentially and remains valid even if you peek at data multiple times.",
        tags$br(),
        tags$strong("p-value:"),
        " A classical fixed-sample test. It is straightforward to interpret but does not permit peeking without inflating Type I error.",
        tags$br(), tags$br(),
        "At a given sample size, a smaller p-value and a larger E-value both indicate stronger evidence."
      ),
      tags$h6("Confidence Interval", class = "fw-bold mt-3 mb-2"),
      tags$p(
        "The 95% confidence interval (Clopper\u2013Pearson exact) gives a range of response rates ",
        "that are consistent with your observed data. If the interval excludes p\u2080, that is mild evidence ",
        "against the null."
      )
    )
  })

  # ── Tab 2: Sequential ──────────────────────────────────────────────────────

  seq_obs <- reactiveVal(NULL)

  # Parse and validate a character vector of "0"/"1" tokens
  validate_binary <- function(parts) {
    parts <- parts[nchar(trimws(parts)) > 0]
    vals  <- suppressWarnings(as.integer(parts))
    if (length(vals) == 0 || any(is.na(vals)) || any(!vals %in% c(0L, 1L)))
      return(NULL)
    vals
  }

  # Text input processing
  observeEvent(input$seq_run, {
    raw   <- gsub("[[:space:]]", "", input$seq_text)
    parts <- strsplit(raw, ",")[[1]]
    v     <- validate_binary(parts)
    if (is.null(v)) {
      showNotification(
        "Could not parse input. Enter only 0s and 1s, comma-separated.",
        type = "error", duration = 6
      )
    }
    seq_obs(v)
  })

  # CSV input processing
  observeEvent(input$seq_run_csv, {
    req(input$seq_file)
    v <- tryCatch({
      # Read all lines, collapse, split on commas/newlines
      raw   <- paste(readLines(input$seq_file$datapath, warn = FALSE), collapse = ",")
      raw   <- gsub("[[:space:]]", "", raw)
      parts <- strsplit(raw, ",")[[1]]
      # Keep only tokens that are literally "0" or "1" (silently drops text headers)
      parts <- parts[parts %in% c("0", "1")]
      validate_binary(parts)
    }, error = function(e) NULL)

    if (is.null(v)) {
      showNotification(
        "Could not read file. CSV must contain a column of 0s and 1s.",
        type = "error", duration = 6
      )
    }
    seq_obs(v)
  })

  # E-value trajectory
  seq_traj <- reactive({
    req(seq_obs(), input$p1 > input$p0)
    obs   <- seq_obs()
    n     <- length(obs)
    p0    <- input$p0;  p1 <- input$p1
    cum_a <- cumsum(obs)
    cum_b <- seq_len(n) - cum_a
    log_e <- cum_a * log(p1 / p0) + cum_b * log((1 - p1) / (1 - p0))
    data.frame(patient = seq_len(n), cum_a = cum_a, cum_b = cum_b, e = exp(log_e))
  })

  # Final-observation summary statistics
  seq_sum <- reactive({
    traj       <- seq_traj()
    n          <- nrow(traj)
    a          <- traj$cum_a[n]
    b          <- traj$cum_b[n]
    e          <- traj$e[n]
    eff_thresh <- 1 / input$alpha
    fut_thresh <- input$fut_e

    first_eff <- { idx <- which(traj$e >= eff_thresh); if (length(idx)) min(idx) else NA_integer_ }
    first_fut <- { idx <- which(traj$e <= fut_thresh); if (length(idx)) min(idx) else NA_integer_ }

    list(
      a = a, b = b, n = n, rate = a / n, e = e,
      pval = exact_pval(a, n, input$p0),
      ci   = exact_ci(a, n),
      dec  = make_decision(e, n, eff_thresh, fut_thresh, input$max_n),
      eff_thresh = eff_thresh, fut_thresh = fut_thresh,
      first_eff = first_eff, first_fut = first_fut
    )
  })

  output$seq_rate   <- renderText({ s <- seq_sum(); sprintf("%.1f%%", 100 * s$rate) })
  output$seq_ci     <- renderText({
    s <- seq_sum()
    if (any(is.na(s$ci))) "\u2014" else sprintf("(%.3f, %.3f)", s$ci[1], s$ci[2])
  })
  output$seq_e      <- renderText({ sprintf("%.4f", seq_sum()$e) })
  output$seq_inv_e  <- renderText({ sprintf("%.4f", 1 / seq_sum()$e) })
  output$seq_pval   <- renderText({ fmt_pval(seq_sum()$pval) })

  output$seq_decision <- renderUI({
    s <- seq_sum()
    decision_alert(s$dec, s$e, s$n, s$a, s$b, s$eff_thresh, s$fut_thresh,
                   s$first_eff, s$first_fut)
  })

  output$seq_plot <- renderPlot({
    traj       <- seq_traj()
    s          <- seq_sum()
    eff_thresh <- s$eff_thresh
    fut_thresh <- s$fut_thresh

    plt <- ggplot(traj, aes(x = patient, y = e)) +
      geom_hline(yintercept = 1,          color = "grey60",  linewidth = 0.5, linetype = "dotted") +
      geom_hline(yintercept = eff_thresh, color = "#1e8449", linewidth = 1.0, linetype = "dashed") +
      geom_hline(yintercept = fut_thresh, color = "#c0392b", linewidth = 1.0, linetype = "dashed") +
      geom_line( color = "#2980b9", linewidth = 1.1) +
      geom_point(color = "#2980b9", size = 1.8, alpha = 0.7) +
      annotate("label",
               x = 1, y = eff_thresh,
               label = sprintf("Efficacy: E \u2265 %.2f  (1/\u03b1)", eff_thresh),
               hjust = 0, vjust = -0.4,
               color = "#1e8449", fill = "white", label.color = "#1e8449", size = 3.5) +
      annotate("label",
               x = 1, y = fut_thresh,
               label = sprintf("Futility: E \u2264 %.4f  (\u03b1)", fut_thresh),
               hjust = 0, vjust = 1.4,
               color = "#c0392b", fill = "white", label.color = "#c0392b", size = 3.5) +
      scale_y_log10(labels = scales::label_number(big.mark = ",")) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
      labs(
        x       = "Patient Number",
        y       = "E-Value (log scale)",
        caption = sprintf(
          "LR Martingale  |  p\u2080 = %.2f,  p\u2081 = %.2f  |  \u03b1 = %.2f  |  1/\u03b1 = %.2f",
          input$p0, input$p1, input$alpha, eff_thresh
        )
      ) +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(),
            plot.caption     = element_text(color = "grey50"))

    # Mark first efficacy crossing
    if (!is.na(s$first_eff)) {
      plt <- plt +
        geom_vline(xintercept = s$first_eff, color = "#1e8449",
                   linewidth = 0.8, linetype = "dotted") +
        annotate("text",
                 x = s$first_eff, y = eff_thresh,
                 label = sprintf("n = %d", s$first_eff),
                 hjust = -0.2, vjust = 0.5,
                 color = "#1e8449", size = 3.5, fontface = "bold")
    }

    # Mark first futility crossing
    if (!is.na(s$first_fut)) {
      plt <- plt +
        geom_vline(xintercept = s$first_fut, color = "#c0392b",
                   linewidth = 0.8, linetype = "dotted") +
        annotate("text",
                 x = s$first_fut, y = fut_thresh,
                 label = sprintf("n = %d", s$first_fut),
                 hjust = -0.2, vjust = 0.5,
                 color = "#c0392b", size = 3.5, fontface = "bold")
    }

    plt
  })

  output$seq_plot_footer <- renderUI({
    s <- seq_sum()
    cross_note <- ""
    if (!is.na(s$first_eff))
      cross_note <- sprintf(
        " Efficacy boundary (E \u2265 %.2f) first crossed at patient %d.", s$eff_thresh, s$first_eff)
    else if (!is.na(s$first_fut))
      cross_note <- sprintf(
        " Futility boundary (E \u2264 %.4f) first crossed at patient %d.", s$fut_thresh, s$first_fut)

    tags$small(tags$em(htmltools::HTML(paste0(
      "Green dashed line = efficacy threshold (1/\u03b1 = ",
      sprintf("%.2f", s$eff_thresh), "). &ensp;",
      "Red dashed line = futility threshold (\u03b1 = ",
      sprintf("%.4f", s$fut_thresh), "). &ensp;",
      "Grey dotted line = E = 1 (no evidence).",
      cross_note
    ))))
  })

  output$seq_table <- renderTable({
    s <- seq_sum()
    summary_table_df(s$a, s$b, s$n, s$rate, s$ci, s$e,
                     s$pval, s$eff_thresh, s$fut_thresh, s$dec)
  }, striped = TRUE, hover = TRUE, width = "100%")

  output$seq_footer <- renderUI({
    tags$small(tags$em(htmltools::HTML(
      "E-value: LR martingale = (p\u2081/p\u2080)<sup>a</sup>",
      " &times; ((1&minus;p\u2081)/(1&minus;p\u2080))<sup>b</sup>. &ensp;",
      "95% CI: Clopper&ndash;Pearson exact. &ensp;",
      "1/E is the minimum \u03b1 at which you would reject H\u2080 with the current data."
    )))
  })

  output$seq_pval_detail <- renderUI({
    s <- seq_sum()
    tagList(
      tags$h6("What is this?", class = "fw-bold mt-3 mb-2"),
      tags$p(
        "The one-sided p-value is the probability of observing ",
        tags$strong(sprintf("%d or more", s$a)),
        " successes if the true response rate were exactly ",
        tags$strong(sprintf("p\u2080 = %.2f%%", 100 * input$p0)),
        ". It is a classical hypothesis test (exact binomial test) with:",
        tags$br(),
        "Null hypothesis: p = p\u2080", tags$br(),
        "Alternative hypothesis: p > p\u2080"
      ),
      tags$h6("Interpretation", class = "fw-bold mt-3 mb-2"),
      tags$ul(
        tags$li(
          tags$strong("Very small p-value (e.g., < 0.001):"),
          " Strong evidence against the null. Your data are unlikely if the null were true."
        ),
        tags$li(
          tags$strong("Moderate p-value (e.g., 0.01\u20130.05):"),
          " Suggestive evidence against the null, consistent with the E-value approach."
        ),
        tags$li(
          tags$strong("Large p-value (e.g., > 0.10):"),
          " Weak evidence against the null. The data are consistent with p = p\u2080."
        )
      ),
      tags$h6("How does it relate to the E-value?", class = "fw-bold mt-3 mb-2"),
      tags$p(
        "Both the p-value and E-value measure evidence against H\u2080, but differently:",
        tags$br(),
        tags$strong("E-value:"),
        " A likelihood ratio that can be computed sequentially and remains valid even if you peek at data multiple times.",
        tags$br(),
        tags$strong("p-value:"),
        " A classical fixed-sample test. It is straightforward to interpret but does not permit peeking without inflating Type I error.",
        tags$br(), tags$br(),
        "At a given sample size, a smaller p-value and a larger E-value both indicate stronger evidence."
      ),
      tags$h6("Confidence Interval", class = "fw-bold mt-3 mb-2"),
      tags$p(
        "The 95% confidence interval (Clopper\u2013Pearson exact) gives a range of response rates ",
        "that are consistent with your observed data. If the interval excludes p\u2080, that is mild evidence ",
        "against the null."
      )
    )
  })

  # ── Tab 3: Operating Characteristics ────────────────────────────────────────

  oc_results <- eventReactive(input$oc_run, {
    req(input$p1 > input$p0)

    p0    <- input$p0
    p1    <- input$p1
    al    <- input$alpha
    max_n <- input$max_n
    n_s   <- input$oc_sims

    # Grid spanning from below p0 to above p1; always include p0 and p1 exactly
    p_lo   <- max(0.02, p0 - 0.15)
    p_hi   <- min(0.97, p1 + 0.15)
    p_grid <- sort(unique(c(seq(p_lo, p_hi, by = 0.025), p0, p1)))

    withProgress(message = "Simulating OC curves \u2014 please wait...", value = 0, {
      rows <- lapply(seq_along(p_grid), function(i) {
        incProgress(1 / length(p_grid),
                    detail = sprintf("p = %.3f  (%d / %d)", p_grid[i], i, length(p_grid)))
        res <- simulate_e_oc(p_grid[i], p0, p1, al, max_n, n_sims = n_s)
        data.frame(
          true_p            = p_grid[i],
          prob_efficacy     = res$prob_efficacy,
          prob_futility     = res$prob_futility,
          prob_inconclusive = res$prob_inconclusive,
          ess               = res$ess
        )
      })
    })

    do.call(rbind, rows)
  })

  # Helper: extract row nearest to a target p value
  oc_at <- function(df, p_val) df[which.min(abs(df$true_p - p_val)), ]

  # Summary at p0
  output$oc_h0_summary <- renderUI({
    row <- oc_at(oc_results(), input$p0)
    mk  <- function(label, val)
      tags$tr(tags$td(class = "text-muted small pe-3", label),
               tags$td(tags$strong(val)))
    tagList(
      tags$table(
        class = "table table-sm table-borderless mb-1",
        tags$tbody(
          mk("Type I Error",      sprintf("%.4f", row$prob_efficacy)),
          mk("P(futility stop)",  sprintf("%.4f", row$prob_futility)),
          mk("P(inconclusive)",   sprintf("%.4f", row$prob_inconclusive)),
          mk("E[N]",              sprintf("%.1f",  row$ess))
        )
      ),
      tags$small(class = "text-muted fst-italic",
                 sprintf("Efficacy boundary = 1/\u03b1 = %.2f  |  Futility boundary = E \u2264 %.4f",
                         1 / input$alpha, input$fut_e))
    )
  })

  # Summary at p1
  output$oc_h1_summary <- renderUI({
    row <- oc_at(oc_results(), input$p1)
    mk  <- function(label, val)
      tags$tr(tags$td(class = "text-muted small pe-3", label),
               tags$td(tags$strong(val)))
    tagList(
      tags$table(
        class = "table table-sm table-borderless mb-1",
        tags$tbody(
          mk("Power",             sprintf("%.4f", row$prob_efficacy)),
          mk("P(futility stop)",  sprintf("%.4f", row$prob_futility)),
          mk("P(inconclusive)",   sprintf("%.4f", row$prob_inconclusive)),
          mk("E[N]",              sprintf("%.1f",  row$ess))
        )
      ),
      tags$small(class = "text-muted fst-italic",
                 sprintf("Efficacy boundary = 1/\u03b1 = %.2f  |  Futility boundary = E \u2264 %.4f",
                         1 / input$alpha, input$fut_e))
    )
  })

  # Decision probabilities plot
  output$oc_plot_probs <- renderPlot({
    df   <- oc_results()
    p0   <- input$p0
    p1   <- input$p1

    # Reshape to long for ggplot
    long <- rbind(
      data.frame(true_p = df$true_p, prob = df$prob_efficacy,
                 outcome = "Efficacy (reject H\u2080)"),
      data.frame(true_p = df$true_p, prob = df$prob_futility,
                 outcome = "Futility stop"),
      data.frame(true_p = df$true_p, prob = df$prob_inconclusive,
                 outcome = "Inconclusive at max N")
    )
    long$outcome <- factor(long$outcome,
                           levels = c("Efficacy (reject H\u2080)",
                                      "Futility stop",
                                      "Inconclusive at max N"))

    cols <- c("Efficacy (reject H\u2080)" = "#1e8449",
              "Futility stop"             = "#c0392b",
              "Inconclusive at max N"     = "#7f8c8d")

    ggplot(long, aes(x = true_p, y = prob, color = outcome)) +
      geom_vline(xintercept = p0, color = "grey40", linewidth = 0.7, linetype = "dashed") +
      geom_vline(xintercept = p1, color = "grey40", linewidth = 0.7, linetype = "dashed") +
      annotate("text", x = p0, y = 1.02, label = sprintf("p\u2080 = %.2f", p0),
               hjust = 0.5, vjust = 0, color = "grey30", size = 3.5) +
      annotate("text", x = p1, y = 1.02, label = sprintf("p\u2081 = %.2f", p1),
               hjust = 0.5, vjust = 0, color = "grey30", size = 3.5) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 2, alpha = 0.7) +
      scale_color_manual(values = cols, name = NULL) +
      scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
      scale_y_continuous(labels = scales::label_percent(accuracy = 1),
                         limits = c(0, 1.08), breaks = seq(0, 1, by = 0.2)) +
      labs(
        x       = "True Response Rate (p)",
        y       = "Probability",
        caption = sprintf(
          "LR Martingale  |  p\u2080 = %.2f, p\u2081 = %.2f  |  \u03b1 = %.2f  |  N = %d  |  %s sims per point",
          p0, p1, input$alpha, input$max_n,
          formatC(input$oc_sims, format = "d", big.mark = ",")
        )
      ) +
      theme_minimal(base_size = 13) +
      theme(legend.position   = "bottom",
            panel.grid.minor  = element_blank(),
            plot.caption      = element_text(color = "grey50"))
  })

  output$oc_plot_footer <- renderUI({
    tags$small(tags$em(
      "Green = probability of declaring efficacy (= Type I error under p\u2080, power under p\u2081). ",
      "Red = probability of stopping for futility. ",
      "Grey = probability of reaching maximum N without a decision."
    ))
  })

  # ESS plot
  output$oc_plot_ess <- renderPlot({
    df    <- oc_results()
    p0    <- input$p0
    p1    <- input$p1
    max_n <- input$max_n

    ggplot(df, aes(x = true_p, y = ess)) +
      geom_hline(yintercept = max_n, color = "grey60",
                 linewidth = 0.6, linetype = "dotted") +
      geom_vline(xintercept = p0, color = "grey40",
                 linewidth = 0.7, linetype = "dashed") +
      geom_vline(xintercept = p1, color = "grey40",
                 linewidth = 0.7, linetype = "dashed") +
      annotate("text", x = p0, y = max_n * 1.02,
               label = sprintf("p\u2080 = %.2f", p0),
               hjust = 0.5, vjust = 0, color = "grey30", size = 3.5) +
      annotate("text", x = p1, y = max_n * 1.02,
               label = sprintf("p\u2081 = %.2f", p1),
               hjust = 0.5, vjust = 0, color = "grey30", size = 3.5) +
      annotate("text", x = min(df$true_p), y = max_n,
               label = sprintf("Max N = %d", max_n),
               hjust = 0, vjust = -0.4, color = "grey50", size = 3.2) +
      geom_line(color = "#2980b9", linewidth = 1.2) +
      geom_point(color = "#2980b9", size = 2, alpha = 0.7) +
      scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
      scale_y_continuous(limits = c(0, max_n * 1.08),
                         breaks = scales::pretty_breaks(n = 6)) +
      labs(
        x       = "True Response Rate (p)",
        y       = "Expected Sample Size  E[N]",
        caption = sprintf(
          "LR Martingale  |  p\u2080 = %.2f, p\u2081 = %.2f  |  \u03b1 = %.2f  |  N = %d  |  %s sims per point",
          p0, p1, input$alpha, input$max_n,
          formatC(input$oc_sims, format = "d", big.mark = ",")
        )
      ) +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(),
            plot.caption     = element_text(color = "grey50"))
  })
}

shinyApp(ui, server)
