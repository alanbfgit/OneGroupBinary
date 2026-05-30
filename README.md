# One-Group Binary Response Trial Design Tool

A comprehensive Shiny application for designing and analyzing single-arm clinical trials with binary response outcomes (e.g., Phase II oncology trials). This tool implements multiple design methods and provides interactive visualizations for design operating characteristics and Bayesian posterior probability calculations.

## Features

### 1. **Simon's Two-Stage Design** (Tab 1)
Classical frequentist design for single-arm trials with binary endpoints.

- **Parameters:**
  - Undesirable response rate (p₀): null hypothesis
  - Acceptable response rate (p₁): alternative hypothesis
  - Type I error (α)
  - Power (1 - β)
  
- **Design Options:**
  - **Optimal**: Minimizes expected sample size under p₀
  - **MiniMax**: Minimizes maximum sample size (n)
  
- **Outputs:**
  - Stage 1 and total sample sizes
  - Stopping thresholds (r₁ and r)
  - Operating characteristics (EN, PET)
  - Decision rules with posterior probability calculations

---

### 2. **Posterior Distribution Visualization** (Tab 2)
Interactive Bayesian analysis of observed data using a Beta prior.

- **Inputs:**
  - Total observations (n)
  - Observed successes (a)
  - Beta prior parameters (α, β) — default Beta(1,1) is uniform
  
- **Outputs:**
  - Posterior Beta distribution
  - P(p ≤ p₀ | data) and P(p ≥ p₁ | data)
  - Visualization with prior overlay, thresholds, and region shading

---

### 3. **Thall & Simon Bayesian Sequential Design** (Tab 3)
Adaptive Bayesian design with interim analyses and stopping boundaries.

- **Parameters:**
  - Maximum sample size (N)
  - Cohort size for interim looks
  - Futility threshold: P(p > p₀ | data)
  - Efficacy threshold: P(p > p₀ | data)
  
- **Outputs:**
  - Stopping boundaries table
  - Operating characteristics from 5,000 simulations:
    - Type I error and power
    - Expected sample size (ESS) under null and alternative
    - Probability of early stopping

---

### 4. **Bayesian Sequential Design with Threshold Calibration** (Tab 4)
Fully calibrated Bayesian design using grid search to find optimal stopping thresholds.

- **Process:**
  - Automatically searches grid of efficacy (u) and futility (l) thresholds
  - Calibrates to match target Type I error and power
  - Performs 5,000 confirmatory simulations
  
- **Outputs:**
  - Calibrated boundary plot with expected trajectories
  - Operating characteristics under H₀ and H₁
  - Real-time decision tool for observed data:
    - Sequential decision (stop for efficacy, futility, or continue)
    - Posterior probability calculation
    - Visual boundary comparison

---

## Installation

### Requirements
- R 4.0+
- Shiny, bslib (Bootstrap 5 UI framework)
- clinfun (for Simon's two-stage design calculations)
- ggplot2, dplyr, scales
- bsicons (Bootstrap icons)

### Setup
```r
# Install required packages
pkgs <- c("shiny", "bslib", "clinfun", "ggplot2", "dplyr", "scales", "bsicons")
install.packages(pkgs)

# Run the app
shiny::runApp("One_Group_Binary_v1.R")
```

Alternatively, use RStudio to open the file and click **"Run App"**.

---

## Usage Guide

### Basic Workflow

1. **Set Design Parameters** (Sidebar)
   - Specify p₀ (undesirable rate), p₁ (acceptable rate)
   - Set α (Type I error) and power (1 - β)
   - Choose prior: Beta(α, β) for posterior calculations
   
2. **Select Design Method** (Tabs 1-4)
   - **Tab 1**: Use Simon's design for classical frequentist planning
   - **Tab 2**: Visualize posterior with sample data to understand Bayesian inference
   - **Tab 3**: Set up adaptive trial with fixed interim looks
   - **Tab 4**: Calibrate thresholds for fully adaptive monitoring

3. **Interpret Results**
   - Review decision rules and operating characteristics
   - Use decision boundary plots to understand sample paths
   - Monitor real observed data using the sequential decision tool (Tab 4)

### Example: Phase II Lung Cancer Trial
- **p₀ = 0.20** (undesirable response rate)
- **p₁ = 0.40** (acceptable response rate)
- **α = 0.10, Power = 0.90**
- **Prior**: Beta(1, 1) (uninformative)

**Simon Optimal Design Result:**
- Stage 1: Enroll 23 patients; stop if ≤ 3 responses
- Stage 2: Enroll 19 more patients (total 42)
- Declare promising if > 10 total responses

---

## Key Concepts

### Posterior Probability Interpretation
The posterior probability P(p > p₀ | data) combines:
- **Prior belief** about the true response rate (Beta prior)
- **Observed data** (number of successes and failures)

This provides a direct probability statement about the parameter, useful for decision-making.

### Operating Characteristics
- **Type I Error**: False positive rate (probability of success when treatment is inactive, p = p₀)
- **Power**: True positive rate (probability of success when treatment is active, p = p₁)
- **ESS**: Expected number of patients enrolled (useful for planning)
- **PET**: Probability of early termination under p₀

### Bayesian Adaptive Designs
Provide efficiency through:
- **Continuous monitoring**: Evaluate after each patient
- **Early stopping**: Declare efficacy or futility without enrolling all N patients
- **Flexible thresholds**: Boundaries adjust to maintain Type I error control

---

## Technical Details

### Simon's Two-Stage Design
Calculated using `clinfun::ph2simon()`. Returns all designs meeting power and Type I error constraints; app selects either Optimal (min EN) or MiniMax (min N).

### Bayesian Sequential Design Calibration (Tab 4)
- **Grid search**: Tests 77 combinations of (u, l) thresholds
- **Each cell**: 2,000 simulations under both H₀ and H₁
- **Selection criterion**: Minimize |Type I error - α| + |Power - (1-β)|
- **Confirmation**: Final 5,000 simulations at best (u, l)

---

## File Information

- **File**: `One_Group_Binary_v1.R`
- **Language**: R with Shiny
- **UI Framework**: Bootstrap 5 (bslib)
- **Lines of Code**: ~1,180
- **Author**: Genelux - Forsythe 2026

---

## References

1. **Simon R** (1989). Optimal two-stage designs for phase II clinical trials. Controlled Clinical Trials, 10(1), 1-10.

2. **Thall PF & Simon R** (1994). Practical Bayesian guidelines for phase IIb clinical trials. Biometrics, 50(2), 337-349.

3. **Fleiss JL, Levin B, Paik MC** (2003). Statistical Methods for Rates and Proportions (3rd ed.). Wiley.

---

## License

This tool is provided for educational and research purposes. For regulatory submissions, consult with biostatisticians and regulatory affairs specialists.

---

## Support & Contributions

For questions, bug reports, or feature requests, please open an issue or contact the development team.

**Repository**: [alanbfgit](https://github.com/alanbfgit/One-Group-Binary)

---

## Citation

If you use this tool in your research, please cite:

> Forsythe A. (2026). One-Group Binary Response Trial Design Tool. R Shiny Application. Genelux.

