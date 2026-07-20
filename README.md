# ZIPXG: Zero Inflated Poisson XGamma Regression Model

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

This repository contains the complete R code for the paper:

**"A New Zero Inflated Poisson XGamma Distribution and its Regression Model: Estimation, Simulation and Application"**  
by Yasin Asar and Caner Tanış  
*The Journal of Supercomputing*

---

## Overview

The ZIPXG model is a novel mixture distribution designed to handle count data characterized by excessive zeros and significant overdispersion. This repository provides all code necessary to reproduce the simulation studies, real-data analysis, and figures presented in the paper.

The methodology employs a hybrid GA+BFGS optimization strategy to handle the non-convex likelihood surface, implemented using parallel computing for large-scale Monte Carlo simulations.

---

## Repository Structure
ZIPXG/
├── README.md # This file
├── session_info.txt # R session information for reproducibility
├── LICENSE # MIT License
│
├── multiple_start_diagnostic_results.csv # Results from convergence diagnostics
├── scenario_C_ZINB_results.csv # Results from misspecification simulation
│
├── sim4_1.R # Simulation for ZIPXG distribution (Section 4.1)
├── sim4_2.R # Main simulation for regression model (Section 4.2)
├── sim4_2_dgp_zinb.R # Misspecification scenario (Scenario C)
├── AffairsDataAnalysis.R # Affairs data analysis (Section 5)
├── VuongTests_Rootograms.R # Vuong tests and rootograms
├── SensitivityAnalysis.R # GA hyperparameter sensitivity analysis
├── OptimizerEffect.R # Disentangling optimizer effect
├── multiple_start_diagnostic.R # Global convergence diagnostics (Section 4.2.3)
├── benchmark.R # Parallel performance benchmarking
├── timings.R # Model fitting timings
└── test_likelihood.R # Unit test for numerical stability


---

## File Descriptions and Corresponding Paper Sections

### Core Simulation Files

| File | Description | Paper Section |
|------|-------------|---------------|
| `sim4_1.R` | Simulation for the ZIPXG distribution parameters ($\theta$ and $\omega$). Evaluates Bias, MSE, and MRE across sample sizes. | Section 4.1 |
| `sim4_2.R` | Main simulation for the ZIPXG regression model (Scenarios A and B). Includes out-of-sample MSE and AIC/BIC comparisons. | Section 4.2 |
| `sim4_2_dgp_zinb.R` | Misspecification scenario (Scenario C). Data generated from ZINB distribution to test robustness. | Section 4.2 (Scenario C) |

### Real Data Analysis

| File | Description | Paper Section |
|------|-------------|---------------|
| `AffairsDataAnalysis.R` | Full analysis of the Affairs dataset. Includes variable selection, model fitting, AIC/BIC comparison, coefficient tables, and out-of-sample MSE. | Section 5 |
| `VuongTests_Rootograms.R` | Vuong non-nested hypothesis tests and rootogram visualizations for all models. | Section 5 |

### Diagnostic and Supplementary Analyses

| File | Description | Paper Section |
|------|-------------|---------------|
| `SensitivityAnalysis.R` | Sensitivity analysis of GA hyperparameters (light/moderate/heavy settings). | Appendix A |
| `OptimizerEffect.R` | Disentangles optimizer effect vs. model effect. Compares ZIPXG and ZINB using the same GA+BFGS optimizer. | Section 4.2.2 |
| `multiple_start_diagnostic.R` | Global convergence diagnostics. Runs 100 independent optimizations with different random seeds. | Section 4.2.3 |
| `benchmark.R` | Parallel performance benchmarking. Measures speedup and efficiency across 1, 2, 4, 8, and 12 cores. | Section 4.2.1 |
| `timings.R` | Model fitting timings for all six models (ZIP, ZINB, ZGeo, HP, HNB, ZIPXG). | Appendix C |
| `test_likelihood.R` | Unit test verifying numerical stability of the corrected log-likelihood function. | Code Availability Statement |

---

## System Requirements

### Hardware

The simulations were performed on a MacBook Pro with Apple M4 Pro processor (12 physical cores) and 24 GB RAM. The code can run on any system with R installed, but parallel execution benefits from multiple cores.

### Software

- **R version:** 4.5.2 or higher
- **Operating System:** macOS / Linux / Windows

### Required R Packages

```r
install.packages(c(
    "MASS",           # Multivariate normal generation
    "pscl",           # Zero-inflated and hurdle models
    "doParallel",     # Parallel computing
    "foreach",        # Parallel loops
    "GA",             # Genetic algorithm
    "numDeriv",       # Numerical derivatives
    "knitr",          # LaTeX table generation
    "kableExtra",     # LaTeX table styling
    "dplyr",          # Data manipulation
    "ggplot2",        # Plotting
    "gridExtra",      # Grid plotting
    "AER",            # Affairs dataset
    "glmmTMB",        # Alternative GLMM fitting (optional)
    "countreg"        # Rootograms (optional)
))
```

### Session Information

R version 4.5.2 (2025-10-31)
Platform: aarch64-apple-darwin20

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] dplyr_1.1.4         pscl_1.5.9          numDeriv_2016.8-1.1 kableExtra_1.4.0   
 [5] knitr_1.50          ggplot2_4.0.1       doParallel_1.0.17   GA_3.2.5           
 [9] iterators_1.0.14    foreach_1.5.2       MASS_7.3-65
 
 
Full session information is provided in session_info.txt

### Reproducing the Results
 