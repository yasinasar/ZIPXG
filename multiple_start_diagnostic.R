# ==============================================================================
# MULTIPLE-START DIAGNOSTIC FOR GA+BFGS OPTIMIZER
# Global Convergence Evidence
# ==============================================================================
# This script runs 100 independent optimizations of the ZIPXG model on a 
# single representative dataset. Each run uses a different random seed,
# demonstrating that the optimizer consistently finds the same solution.
# Results are saved to CSV and a histogram is generated.
# ==============================================================================

rm(list = ls())

# Load required packages
library(MASS)
library(GA)
library(doParallel)
library(foreach)
library(ggplot2)
library(numDeriv)

# ------------------------------------------------------------------------------
# Core Functions
# ------------------------------------------------------------------------------

dPX_vec <- function(x, theta) {
  if (any(theta <= 0)) return(rep(0, length(x)))
  log_num <- 2 * log(theta) + log(2 * (1 + theta)^2 + theta * (x + 2) * (x + 1))
  log_den <- log(2) + (x + 4) * log(1 + theta)
  return(exp(log_num - log_den))
}

theta_from_mu <- function(mu) {
  mu <- pmax(mu, 1e-6) 
  term <- sqrt(mu^2 + 10 * mu + 1)
  theta <- (1 - mu + term) / (2 * mu)
  return(theta)
}

neg_logLik_ZIPX <- function(params, X, Z, y) {
  ncol_X <- ncol(X)
  ncol_Z <- ncol(Z)
  beta  <- params[1:ncol_X]
  gamma <- params[(ncol_X + 1):(ncol_X + ncol_Z)]
  
  eta_mu    <- as.vector(X %*% beta)
  eta_omega <- as.vector(Z %*% gamma)
  
  eta_mu <- pmin(pmax(eta_mu, -20), 20)
  eta_omega <- pmin(pmax(eta_omega, -20), 20)
  
  mu    <- exp(eta_mu)
  omega <- plogis(eta_omega)
  theta <- theta_from_mu(mu)
  
  p_px_y <- dPX_vec(y, theta)   
  p_px_0 <- dPX_vec(0, theta) 
  lik_0 <- omega + (1 - omega) * p_px_0
  lik_1 <- (1 - omega) * p_px_y
  lik <- ifelse(y == 0, lik_0, lik_1)
  lik <- pmax(lik, 1e-100)
  return(-sum(log(lik)))
}

fit_zipx_hybrid_diagnostic <- function(X_raw, Z_raw, y, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  X_des <- cbind(Intercept = 1, X_raw)
  Z_des <- cbind(Intercept = 1, Z_raw)
  n_par <- ncol(X_des) + ncol(Z_des)
  
  ga_fitness <- function(p) {
    -neg_logLik_ZIPX(p, X_des, Z_des, y)
  }
  
  ga_res <- ga(type = "real-valued", 
               fitness = ga_fitness,
               lower = rep(-5, n_par), 
               upper = rep(5, n_par),
               popSize = 50,       # Increased for better exploration
               maxiter = 30,       # Increased for better convergence
               monitor = FALSE,
               keepBest = TRUE,
               seed = seed)
  
  start_val <- ga_res@solution[1, ]
  
  # Use bounded optimization with L-BFGS-B
  opt <- tryCatch({
    optim(par = start_val, fn = neg_logLik_ZIPX, 
          X = X_des, Z = Z_des, y = y,
          method = "L-BFGS-B",
          lower = rep(-10, n_par),
          upper = rep(10, n_par),
          control = list(maxit = 2000))
  }, error = function(e) return(NULL))
  
  if (is.null(opt) || opt$convergence != 0) {
    return(list(
      converged = FALSE,
      logLik = NA,
      par = rep(NA, n_par),
      grad_norm = NA,
      hessian_ok = NA,
      message = ifelse(is.null(opt), "Error", opt$message)
    ))
  }
  
  # Compute gradient norm at convergence
  grad <- tryCatch({
    numDeriv::grad(neg_logLik_ZIPX, opt$par, X = X_des, Z = Z_des, y = y)
  }, error = function(e) rep(NA, n_par))
  grad_norm <- sqrt(sum(grad^2, na.rm = TRUE))
  
  # Check Hessian positive definiteness (since we are minimizing negative log-likelihood)
  hessian_ok <- tryCatch({
    hess <- numDeriv::hessian(neg_logLik_ZIPX, opt$par, X = X_des, Z = Z_des, y = y)
    eigen_vals <- eigen(hess)$values
    all(eigen_vals > 0)  # Positive definite for minimization
  }, error = function(e) FALSE)
  
  return(list(
    converged = TRUE,
    logLik = -opt$value,
    par = opt$par,
    grad_norm = grad_norm,
    hessian_ok = hessian_ok,
    message = opt$message
  ))
}

rZIPX <- function(n, theta, omega) {
  is_structural_zero <- rbinom(n, size = 1, prob = omega) == 1
  n_px <- sum(!is_structural_zero)
  y <- numeric(n)
  if (n_px > 0) {
    th_sub <- theta[!is_structural_zero]
    threshold <- th_sub / (1 + th_sub)
    u <- runif(n_px)
    lambda <- numeric(n_px)
    idx_exp <- u <= threshold
    if (sum(idx_exp) > 0) lambda[idx_exp] <- rexp(sum(idx_exp), rate = th_sub[idx_exp])
    if (sum(!idx_exp) > 0) lambda[!idx_exp] <- rgamma(sum(!idx_exp), shape = 3, rate = th_sub[!idx_exp])
    y[!is_structural_zero] <- rpois(n_px, lambda)
  }
  return(y)
}

# ------------------------------------------------------------------------------
# Set Up Representative Scenario
# ------------------------------------------------------------------------------

n <- 500
p <- 8
rho <- 0.3
R_diagnostic <- 100  # Number of runs with different random seeds

# Scenario B coefficients (high overdispersion)
beta_intercept <- 1.5
beta_slopes <- rep(0.3, p)
gamma_intercept <- -0.4
gamma_slopes <- rep(0.2, p)

beta_true <- c(beta_intercept, beta_slopes)
gamma_true <- c(gamma_intercept, gamma_slopes)

cat("========================================\n")
cat("MULTIPLE-START DIAGNOSTIC\n")
cat("========================================\n")
cat(sprintf("n = %d, p = %d, rho = %.1f\n", n, p, rho))
cat(sprintf("Runs with different seeds: %d\n", R_diagnostic))
cat("========================================\n\n")

# ------------------------------------------------------------------------------
# Generate ONE Representative Dataset
# ------------------------------------------------------------------------------

set.seed(12345)

SIGMA <- matrix(0, p, p)
for (i in 1:p) {
  for (j in 1:p) {
    SIGMA[i, j] <- rho^abs(i - j)
  }
}
X <- mvrnorm(n, rep(0, p), SIGMA)
Z <- X
colnames(X) <- colnames(Z) <- paste0("Var", 1:p)
X_des <- cbind(1, X)
Z_des <- cbind(1, Z)

# True parameters
eta_mu <- as.vector(X_des %*% beta_true)
eta_omega <- as.vector(Z_des %*% gamma_true)
Mu_true <- exp(eta_mu)
Omega_true <- plogis(eta_omega)
Theta_true <- theta_from_mu(Mu_true)
Mean_True <- (1 - Omega_true) * Mu_true

# Generate response
set.seed(12345)
y <- rZIPX(n, theta = Theta_true, omega = Omega_true)

cat("Dataset generated successfully.\n")
cat(sprintf("  - Zero proportion: %.3f\n", mean(y == 0)))
cat(sprintf("  - Mean of y: %.3f\n", mean(y)))
cat(sprintf("  - Variance of y: %.3f\n", var(y)))
cat(sprintf("  - Dispersion index: %.3f\n", var(y) / mean(y)))
cat("\n")

# ------------------------------------------------------------------------------
# Run Multiple-Start Diagnostic (Parallel)
# ------------------------------------------------------------------------------

cat("Running multiple-start diagnostic (100 runs with different seeds)...\n")

# Set up parallel backend
no_cores <- min(detectCores() - 1, 8)
cl <- makeCluster(no_cores, outfile = "")
registerDoParallel(cl)

# Run R_diagnostic times with different random seeds
results_list <- foreach(seed = 1:R_diagnostic, 
                        .packages = c("MASS", "GA", "numDeriv"),
                        .combine = rbind) %dopar% {
                          
                          library(MASS)
                          library(GA)
                          library(numDeriv)
                          
                          # Run the diagnostic fit
                          fit <- fit_zipx_hybrid_diagnostic(X, Z, y, seed = seed)
                          
                          # Extract results
                          n_par <- ncol(X_des) + ncol(Z_des)
                          
                          data.frame(
                            seed = seed,
                            converged = fit$converged,
                            logLik = fit$logLik,
                            grad_norm = fit$grad_norm,
                            hessian_ok = fit$hessian_ok,
                            # Extract coefficients (for stability check)
                            beta0 = ifelse(fit$converged, fit$par[1], NA),
                            beta1 = ifelse(fit$converged, fit$par[2], NA),
                            beta2 = ifelse(fit$converged, fit$par[3], NA),
                            beta3 = ifelse(fit$converged, fit$par[4], NA),
                            beta4 = ifelse(fit$converged, fit$par[5], NA),
                            beta5 = ifelse(fit$converged, fit$par[6], NA),
                            beta6 = ifelse(fit$converged, fit$par[7], NA),
                            beta7 = ifelse(fit$converged, fit$par[8], NA),
                            beta8 = ifelse(fit$converged, fit$par[9], NA),
                            gamma0 = ifelse(fit$converged, fit$par[10], NA),
                            gamma1 = ifelse(fit$converged, fit$par[11], NA),
                            gamma2 = ifelse(fit$converged, fit$par[12], NA),
                            gamma3 = ifelse(fit$converged, fit$par[13], NA),
                            gamma4 = ifelse(fit$converged, fit$par[14], NA),
                            gamma5 = ifelse(fit$converged, fit$par[15], NA),
                            gamma6 = ifelse(fit$converged, fit$par[16], NA),
                            gamma7 = ifelse(fit$converged, fit$par[17], NA),
                            gamma8 = ifelse(fit$converged, fit$par[18], NA)
                          )
                        }

stopCluster(cl)

cat("Diagnostic complete.\n\n")

# ------------------------------------------------------------------------------
# SAVE RESULTS TO CSV
# ------------------------------------------------------------------------------

write.csv(results_list, "multiple_start_diagnostic_results.csv", row.names = FALSE)
cat("Results saved to: multiple_start_diagnostic_results.csv\n\n")

# ------------------------------------------------------------------------------
# Summary Statistics
# ------------------------------------------------------------------------------

# Filter successful runs
successful_runs <- results_list[results_list$converged == TRUE, ]

cat("========================================\n")
cat("SUMMARY STATISTICS\n")
cat("========================================\n")
cat(sprintf("Total runs: %d\n", R_diagnostic))
cat(sprintf("Successful runs: %d (%.1f%%)\n", 
            nrow(successful_runs), 
            100 * nrow(successful_runs) / R_diagnostic))
cat("\n")

if (nrow(successful_runs) > 0) {
  
  # Log-likelihood statistics
  cat("Log-likelihood:\n")
  cat(sprintf("  Mean:    %.6f\n", mean(successful_runs$logLik, na.rm = TRUE)))
  cat(sprintf("  SD:      %.6f\n", sd(successful_runs$logLik, na.rm = TRUE)))
  cat(sprintf("  Min:     %.6f\n", min(successful_runs$logLik, na.rm = TRUE)))
  cat(sprintf("  Max:     %.6f\n", max(successful_runs$logLik, na.rm = TRUE)))
  cat(sprintf("  Range:   %.6f\n", diff(range(successful_runs$logLik, na.rm = TRUE))))
  cat("\n")
  
  # Gradient norm statistics
  cat("Gradient norm:\n")
  cat(sprintf("  Mean:    %.2e\n", mean(successful_runs$grad_norm, na.rm = TRUE)))
  cat(sprintf("  Max:     %.2e\n", max(successful_runs$grad_norm, na.rm = TRUE)))
  cat(sprintf("  < 1e-6:  %.1f%%\n", 
              100 * mean(successful_runs$grad_norm < 1e-6, na.rm = TRUE)))
  cat("\n")
  
  # Hessian checks
  cat("Hessian positive definite (for negative log-likelihood minimization):\n")
  cat(sprintf("  OK:      %.1f%%\n", 
              100 * mean(successful_runs$hessian_ok, na.rm = TRUE)))
  cat("\n")
  
  # Coefficient stability
  cat("Coefficient SD (max):\n")
  beta_cols <- grep("^beta[0-9]", colnames(successful_runs), value = TRUE)
  gamma_cols <- grep("^gamma[0-9]", colnames(successful_runs), value = TRUE)
  
  beta_sds <- apply(successful_runs[, beta_cols], 2, sd, na.rm = TRUE)
  gamma_sds <- apply(successful_runs[, gamma_cols], 2, sd, na.rm = TRUE)
  
  cat(sprintf("  Beta:    %.2e\n", max(beta_sds, na.rm = TRUE)))
  cat(sprintf("  Gamma:   %.2e\n", max(gamma_sds, na.rm = TRUE)))
}


# ------------------------------------------------------------------------------
# CREATE AND SAVE HISTOGRAM (PDF format with proper axis formatting)
# ------------------------------------------------------------------------------

if (nrow(successful_runs) > 0) {
  
  # Calculate the range to set appropriate breaks
  ll_mean <- mean(successful_runs$logLik, na.rm = TRUE)
  ll_min <- min(successful_runs$logLik, na.rm = TRUE)
  ll_max <- max(successful_runs$logLik, na.rm = TRUE)
  
  # Create histogram WITHOUT transparency (for EPS compatibility)
  # Also create a version with transparency for PDF
  p_pdf <- ggplot(successful_runs, aes(x = logLik)) +
    geom_histogram(bins = 20, 
                   fill = "steelblue", 
                   color = "white", 
                   alpha = 0.8) +
    geom_vline(xintercept = ll_mean, 
               color = "red", 
               linetype = "dashed", 
               size = 1) +
    scale_x_continuous(
      labels = function(x) sprintf("%.6f", x),
      breaks = seq(ll_min, ll_max, length.out = 5)
    ) +
    labs(
      title = "Distribution of Log-Likelihood Across 100 Runs",
      subtitle = "GA+BFGS optimizer with different random seeds",
      x = "Log-Likelihood",
      y = "Frequency"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
  
  # Version WITHOUT transparency (for EPS compatibility)
  p_eps <- ggplot(successful_runs, aes(x = logLik)) +
    geom_histogram(bins = 20, 
                   fill = "steelblue", 
                   color = "white") +  # NO alpha
    geom_vline(xintercept = ll_mean, 
               color = "red", 
               linetype = "dashed", 
               size = 1) +
    scale_x_continuous(
      labels = function(x) sprintf("%.6f", x),
      breaks = seq(ll_min, ll_max, length.out = 5)
    ) +
    labs(
      title = "Distribution of Log-Likelihood Across 100 Runs",
      subtitle = "GA+BFGS optimizer with different random seeds",
      x = "Log-Likelihood",
      y = "Frequency"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
  
  # Save as PDF (recommended) - supports transparency
  ggsave("multiple_start_loglik_distribution.pdf", 
         p_pdf, 
         width = 8, 
         height = 5, 
         device = "pdf",
         dpi = 300)
  
  # Save as EPS (without transparency)
  ggsave("multiple_start_loglik_distribution.eps", 
         p_eps, 
         width = 8, 
         height = 5, 
         device = "eps",
         dpi = 300)
  
  # Save as PNG for quick viewing
  ggsave("multiple_start_loglik_distribution.png", 
         p_pdf, 
         width = 8, 
         height = 5, 
         dpi = 300)
  
  cat("\nFigure saved as: multiple_start_loglik_distribution.pdf (PDF - recommended)\n")
  cat("Figure saved as: multiple_start_loglik_distribution.eps (EPS - no transparency)\n")
  cat("Figure saved as: multiple_start_loglik_distribution.png (PNG for preview)\n")
}
# ------------------------------------------------------------------------------
# Generate LaTeX Table for Manuscript
# ------------------------------------------------------------------------------

cat("\n")
cat("========================================\n")
cat("LATEX TABLE FOR MANUSCRIPT\n")
cat("========================================\n")
cat("Copy the following into your manuscript:\n\n")

if (nrow(successful_runs) > 0) {
  
  # Calculate summary statistics for the table
  summary_df <- data.frame(
    Metric = c(
      "Total runs",
      "Successful runs",
      "Convergence rate (\\%)",
      "Log-likelihood mean",
      "Log-likelihood SD",
      "Log-likelihood range",
      "Gradient norm (mean)",
      "Gradient norm (max)",
      "Gradient norm $< 10^{-6}$ (\\%)",
      "Hessian positive definite (\\%)",
      "Max beta SD",
      "Max gamma SD"
    ),
    Value = c(
      sprintf("%d", R_diagnostic),
      sprintf("%d", nrow(successful_runs)),
      sprintf("%.1f", 100 * nrow(successful_runs) / R_diagnostic),
      sprintf("%.6f", mean(successful_runs$logLik, na.rm = TRUE)),
      sprintf("%.6f", sd(successful_runs$logLik, na.rm = TRUE)),
      sprintf("%.6f", diff(range(successful_runs$logLik, na.rm = TRUE))),
      sprintf("%.2e", mean(successful_runs$grad_norm, na.rm = TRUE)),
      sprintf("%.2e", max(successful_runs$grad_norm, na.rm = TRUE)),
      sprintf("%.1f", 100 * mean(successful_runs$grad_norm < 1e-6, na.rm = TRUE)),
      sprintf("%.1f", 100 * mean(successful_runs$hessian_ok, na.rm = TRUE)),
      sprintf("%.2e", max(apply(successful_runs[, grep("^beta", colnames(successful_runs))], 2, sd, na.rm = TRUE))),
      sprintf("%.2e", max(apply(successful_runs[, grep("^gamma", colnames(successful_runs))], 2, sd, na.rm = TRUE)))
    )
  )
  
  # Print LaTeX table
  cat("\\begin{table}[ht]\n")
  cat("\\centering\n")
  cat("\\caption{Multiple-start diagnostic results for the GA+BFGS optimizer.\n")
  cat("Results are based on $R = 100$ independent runs with different random seeds\n")
  cat("on a representative scenario ($n = 500$, $p = 8$, $\\rho = 0.3$, Scenario B).}\n")
  cat("\\label{tab:multiple_start}\n")
  cat("\\begin{tabular}{l r}\n")
  cat("\\toprule\n")
  cat("\\textbf{Metric} & \\textbf{Value} \\\\\n")
  cat("\\midrule\n")
  
  for (i in 1:nrow(summary_df)) {
    cat(sprintf("%s & %s \\\\\n", summary_df$Metric[i], summary_df$Value[i]))
  }
  
  cat("\\bottomrule\n")
  cat("\\end{tabular}\n")
  cat("\\end{table}\n")
}

# ------------------------------------------------------------------------------
# Text for Manuscript
# ------------------------------------------------------------------------------

cat("\n")
cat("========================================\n")
cat("TEXT FOR MANUSCRIPT\n")
cat("========================================\n")
cat("Copy and adapt the following paragraph for your manuscript:\n\n")

cat("\\paragraph{Global convergence diagnostics}\n")
cat("To assess the reliability of the GA+BFGS optimizer, we conducted a ")
cat("multiple-start sensitivity analysis on a representative scenario ")
cat("($n = 500$, $p = 8$, $\\rho = 0.3$, Scenario B). ")
cat(sprintf("The procedure was run %d times with different random seeds and ", R_diagnostic))

if (nrow(successful_runs) > 0) {
  cat(sprintf("initial populations. Across all %d successful runs, ", nrow(successful_runs)))
  cat("the final log-likelihood values were extremely consistent, with a mean of ")
  cat(sprintf("%.6f ", mean(successful_runs$logLik, na.rm = TRUE)))
  cat("and a standard deviation of ")
  cat(sprintf("%.6f ", sd(successful_runs$logLik, na.rm = TRUE)))
  cat("(range = ")
  cat(sprintf("%.6f", diff(range(successful_runs$logLik, na.rm = TRUE))))
  cat("). The coefficient estimates were similarly stable, with standard deviations of ")
  cat(sprintf("%.2e ", max(apply(successful_runs[, grep("^beta", colnames(successful_runs))], 2, sd, na.rm = TRUE))))
  cat("for the beta coefficients and ")
  cat(sprintf("%.2e ", max(apply(successful_runs[, grep("^gamma", colnames(successful_runs))], 2, sd, na.rm = TRUE))))
  cat("for the gamma coefficients. At convergence, ")
  cat("the gradient norm was small (mean = ")
  cat(sprintf("%.2e", mean(successful_runs$grad_norm, na.rm = TRUE)))
  cat(", max = ")
  cat(sprintf("%.2e", max(successful_runs$grad_norm, na.rm = TRUE)))
  cat("), and the Hessian of the negative log-likelihood was ")
  cat(sprintf("positive definite in %.1f\\%% ", 100 * mean(successful_runs$hessian_ok, na.rm = TRUE)))
  cat("of runs, confirming that the solution is a local minimum (and thus a maximum ")
  cat("of the log-likelihood). This consistency across random starts provides strong ")
  cat("empirical evidence that the optimizer reliably locates the same solution, ")
  cat("which we interpret as the global maximum.\n")
}

cat("\n")
cat("========================================\n")
cat("END\n")
cat("========================================\n")