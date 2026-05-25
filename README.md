# OneGroupBinary

**Simon's Two-Stage Design with Bayesian Posterior Analysis**  
*A Shiny application for single-arm Phase II oncology trials*

---

## Overview

This Shiny app implements **Simon's Two-Stage design** for single-arm Phase II clinical trials with a binary endpoint (e.g., tumour response). It is intended for settings where the primary goal is to decide whether a treatment's response rate is promising enough to warrant further study.

The app provides:

- **Optimal and MiniMax** two-stage design parameters via the `clinfun` package
- **Bayesian posterior analysis** of the observed response rate using a conjugate Beta–Binomial model
- Interactive decision rules with posterior probabilities at each design boundary

---

## Design Background

Simon's Two-Stage design (Simon, 1989) controls both:

- **Type I error (α):** the probability of declaring an inactive treatment promising
- **Type II error (β):** the probability of failing to detect an active treatment

Two criteria are supported:

| Criterion | Minimises |
|-----------|-----------|
| **Optimal** | Expected sample size under *p₀* (useful when early stopping is likely) |
| **MiniMax** | Maximum total sample size |

### Stopping Rule

1. Enroll **n₁** patients in Stage 1.  
   - If responses ≤ **r₁**, stop early — the treatment is not sufficiently active.  
   - Otherwise, continue to Stage 2.
2. Enroll the remaining patients (total **N**).  
   - If total responses > **r**, declare the treatment **promising**.  
   - Otherwise, conclude it does not meet the target activity level.

---

## Bayesian Posterior Analysis

Alongside the frequentist Simon design, the app fits a **Beta–Binomial conjugate model**:

$$
p \sim \text{Beta}(\alpha_0, \beta_0) \quad \text{(prior)}
$$

$$
p \mid a, n \sim \text{Beta}(\alpha_0 + a,\ \beta_0 + n - a) \quad \text{(posterior)}
$$

where *a* is the number of observed responses and *n* is the total number of observations.

The default prior `Beta(1, 1)` is uniform (non-informative). Custom priors can encode historical data or clinical belief.

Key posterior summaries reported:

- **P(p ≤ p₀ | data)** — probability the true rate is undesirably low
- **P(p ≥ p₁ | data)** — probability the true rate meets the target

---

## App Features

- **Simon Design tab:** value boxes for stage-1 and final thresholds, full design table, and plain-language decision rules with posterior probabilities at each boundary.
- **Posterior Distribution tab:** interactive sliders for observed data (n, successes), shaded posterior density plot with regions colour-coded relative to p₀ and p₁, overlaid prior density (dashed line), and posterior probability value boxes.

---

## Installation

```r
# Install required packages if not already present
install.packages(c("shiny", "bslib", "bsicons", "ggplot2", "dplyr",
                   "scales", "clinfun"))
```

---

## Usage

```r
shiny::runApp("simon2stage_app_v2.R")
```

Or open `simon2stage_app_v2.R` in RStudio and click **Run App**.

### Input Parameters

| Parameter | Description |
|-----------|-------------|
| **p₀** | Undesirable (null) response rate |
| **p₁** | Acceptable (target) response rate |
| **α** | Type I error rate |
| **1 − β** | Power |
| **Design criterion** | Optimal or MiniMax |
| **Beta prior (α, β)** | Prior on the response rate; default Beta(1,1) is uniform |

---

## Dependencies

| Package | Purpose |
|---------|---------|
| `shiny` | Web application framework |
| `bslib` | Bootstrap 5 UI components and theming |
| `bsicons` | Bootstrap icons |
| `clinfun` | Simon two-stage design calculations (`ph2simon`) |
| `ggplot2` | Posterior density plot |
| `dplyr` | Data wrangling inside plot |
| `scales` | Axis label formatting |

---

## Reference

Simon R. (1989). Optimal two-stage designs for phase II clinical trials. *Controlled Clinical Trials*, **10**(1), 1–10. https://doi.org/10.1016/0197-2456(89)90015-9

---

## Author

Alan Forsythe — Forsythe and Bear LLC, 2026
