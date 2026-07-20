# ZIPXG: Zero Inflated Poisson XGamma Regression Model

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21452901.svg)](https://doi.org/10.5281/zenodo.21452901)

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

├── multiple_start_diagnostic_results.csv # Results from convergence diagnostics

├── scenario_C_ZINB_results.csv # Results from misspecification simulation

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

#### Step 1: Clone the Repository

```r
git clone https://github.com/yasinasar/ZIPXG.git
cd ZIPXG
```


#### Step 2: Install Required Packages

```r
install.packages(c("MASS", "pscl", "doParallel", "foreach", "GA", 
                   "numDeriv", "knitr", "kableExtra", "dplyr", 
                   "ggplot2", "gridExtra", "AER", "glmmTMB"))
# Optional for rootograms:
install.packages("countreg", repos = "https://R-Forge.R-project.org")
```

#### Step 3: Run Analyses

##### 3.1 Simulation for ZIPXG Distribution (Section 4.1)

```r
source("sim4_1.R")
```

Output: ZIPX_Simulation_Grid.pdf (Figure 1)

#####  3.2 Main Simulation for Regression Model (Section 4.2)

Scenario A and B:

```r
source("sim4_2.R")
```

Output: Tables 4, 5, 17, and 18

#####   3.3 Misspecification Scenario (Scenario C)

```r
source("sim4_2_dgp_zinb.R")
```

Output: Tables 6 and 19, scenario_C_ZINB_results.csv

#####  3.4 Sensitivity Analysis (Appendix A)

```r
source("SensitivityAnalysis.R")
```

Output: Table 14


#####  3.5 Optimizer Effect (Section 4.2.2)

```r
source("OptimizerEffect.R")
```

Expected Console Output:

===== Disentangling optimizer effect =====

Mean MSE (ZIPXG): 2.2594

Mean MSE (ZINB): 3.0171

Difference (ZINB - ZIPXG): 0.7577


#####  3.6 Multiple-Start Diagnostic (Section 4.2.3)

```r
source("multiple_start_diagnostic.R")
```

Outputs:

multiple_start_diagnostic_results.csv

multiple_start_loglik_distribution.pdf (Figure 2)

multiple_start_loglik_distribution.eps

Output: Table 7 


#####  3.7 Parallel Performance Benchmarking (Section 4.2.1)

```r
source("benchmark.R")
```

Outputs:

performance_benchmark_results.csv

LaTeX table (speedup and efficiency) printed to console


#####   3.8 Model Fitting Timings (Appendix C)

```r
source("timings.R")
```

Output: Average computational times for each model printed to console (Table 16)

#####   3.9 Affairs Data Analysis (Section 5)

```r
source("AffairsDataAnalysis.R")
```

Outputs: Tables 10, 11 and 12 printed to console.

#####   3.10 Vuong Tests and Rootograms (Section 5)

```r
source("VuongTests_Rootograms.R")
```

Vuong test results (Table 13) printed to console

rootograms_all.eps (Figure 4)

#####   3.11 Numerical Stability Unit Test

```r
source("test_likelihood.R")
```

Expected Output: ALL TESTS PASSED


## Troubleshooting

### Common Issues

#### 1. Parallel execution fails on Windowss

Use makeCluster(no_cores, type = "PSOCK") instead of the default, or set no_cores <- 1 for serial execution.

#### 2. GA package version mismatch

The code was tested with GA version 3.2.5. If you encounter issues, try:

```r
install.packages("GA", version = "3.2.5")
```

#### 3. Memory issues with large simulations

Reduce R (number of replications) in sim4_2.R from 10000 to 1000 for testing. Increase available memory or reduce parallel workers.

#### 4. Missing countreg package

Rootograms require countreg which may not be available on CRAN. Install from R‑Forge:

```r
install.packages("countreg", repos = "https://R-Forge.R-project.org")
```


## Citation

If you use this code in your research, please cite:

@article{asar2026zipxg,

  title={A New Zero Inflated Poisson XGamma Distribution and its Regression Model: Estimation, Simulation and Application},
  
  author={Asar, Yasin and Tanış, Caner},
  
  journal={The Journal of Supercomputing},
  
  year={2026}
}

This version of the code is permanently archived on Zenodo:

https://zenodo.org/badge/DOI/10.5281/zenodo.21452901.svg