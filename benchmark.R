# ==============================================================================
# PERFORMANCE BENCHMARKING FOR PARALLEL COMPUTATION
# ==============================================================================
# This script measures serial and parallel runtime for the simulation framework.
# It uses the same core functions as the main simulation.
# ==============================================================================

rm(list = ls())

# Load required packages
library(MASS)
library(pscl)
library(doParallel)
library(foreach)
library(GA)

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
  ncol_X <- ncol(X); ncol_Z <- ncol(Z)
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

fit_zipx_hybrid <- function(X_raw, Z_raw, y) {
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
               popSize = 20,
               maxiter = 15,
               monitor = FALSE,
               keepBest = TRUE)
  
  start_val <- ga_res@solution[1, ]
  
  opt <- tryCatch({
    optim(par = start_val, fn = neg_logLik_ZIPX, X = X_des, Z = Z_des, y = y,
          method = "BFGS", control = list(maxit = 1000))
  }, error = function(e) return(NULL))
  
  if(is.null(opt)) return(NULL)
  
  logLik_val <- -opt$value
  n_obs <- length(y)
  
  return(list(par = opt$par, 
              logLik = logLik_val, 
              AIC = 2 * n_par - 2 * logLik_val, 
              BIC = n_par * log(n_obs) - 2 * logLik_val))
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

calc_mse_oos_robust <- function(pred_vals, true_vals) {
  if(any(is.na(pred_vals)) || any(pred_vals > 10000)) return(NA)
  return(mean((pred_vals - true_vals)^2, na.rm = TRUE))
}

# ------------------------------------------------------------------------------
# Set up the representative scenario (n=500, p=8, rho=0.3, Scenario B)
# ------------------------------------------------------------------------------

set_overdispersion_coefs <- function(p) {
  beta_intercept <- 1.5
  beta_slopes <- rep(0.3, p)
  gamma_intercept <- -0.4
  gamma_slopes <- rep(0.2, p)
  beta <- c(beta_intercept, beta_slopes)
  gamma <- c(gamma_intercept, gamma_slopes)
  return(list(beta = matrix(beta, ncol = 1), 
              gamma = matrix(gamma, ncol = 1)))
}

# Fixed parameters for the benchmark scenario
n <- 500
p <- 8
rho <- 0.3
R_bench <- 100  # Number of replications for benchmarking (sufficient for timing)

# Generate covariates once (same for all replications)
SIGMA <- matrix(0, p, p)
for (i in 1:p) {
  for (j in 1:p) {
    SIGMA[i, j] <- rho^abs(i - j)
  }
}
X <- mvrnorm(n, rep(0, p), SIGMA)
Z <- X
colnames(X) <- colnames(Z) <- paste0("Var", 1:p)

coefs <- set_overdispersion_coefs(p)
beta.count <- coefs$beta
beta.zero  <- coefs$gamma

X_des <- cbind(1, X)
Z_des <- cbind(1, Z)
eta_mu <- as.vector(X_des %*% beta.count)
eta_omega <- as.vector(Z_des %*% beta.zero)
Mu_true <- exp(eta_mu)
Omega_true <- plogis(eta_omega)
Mean_True <- (1 - Omega_true) * Mu_true

# Pre-generate all random seeds for reproducibility
set.seed(12345)
seeds <- sample(1:1e6, R_bench)

# ------------------------------------------------------------------------------
# Benchmark function: runs one replication and returns out-of-sample MSE for ZIPXG
# (We only need to time the full fitting process; we don't need MSE for benchmarking,
# but we keep the same structure to mimic the actual simulation.)
# ------------------------------------------------------------------------------

run_one_replication <- function(k, X, Z, y, Mean_True) {
  set.seed(seeds[k])
  
  # Train-test split (80/20)
  n_train <- floor(n * 0.8)
  train_idx <- sample(1:n, n_train, replace = FALSE)
  test_idx <- setdiff(1:n, train_idx)
  
  X_train <- X[train_idx, , drop = FALSE]
  X_test  <- X[test_idx, , drop = FALSE]
  Z_train <- Z[train_idx, , drop = FALSE]
  Z_test  <- Z[test_idx, , drop = FALSE]
  y_train <- y[train_idx]
  y_test  <- y[test_idx]
  Mean_True_test <- Mean_True[test_idx]
  
  df_train <- data.frame(y = y_train, X_train, Z_train)
  form <- as.formula(paste("y ~", paste(colnames(X_train), collapse = "+"), "|", 
                           paste(colnames(Z_train), collapse = "+")))
  
  # Fit ZIPXG (the most computationally expensive model)
  res_zipx <- fit_zipx_hybrid(X_train, Z_train, y_train)
  if (is.null(res_zipx)) return(NA)
  
  X_test_des <- cbind(Intercept = 1, X_test)
  Z_test_des <- cbind(Intercept = 1, Z_test)
  ncol_X <- ncol(X_test_des)
  beta_hat  <- res_zipx$par[1:ncol_X]
  gamma_hat <- res_zipx$par[(ncol_X + 1):length(res_zipx$par)]
  pred_zipx <- as.vector((1 - plogis(Z_test_des %*% gamma_hat)) * exp(X_test_des %*% beta_hat))
  
  mse <- calc_mse_oos_robust(pred_zipx, Mean_True_test)
  return(mse)
}

# ------------------------------------------------------------------------------
# Measure runtime for different core counts
# ------------------------------------------------------------------------------

core_counts <- c(1, 2, 4, 8, 12)
results <- data.frame(Cores = core_counts,
                      Runtime = NA,
                      Speedup = NA,
                      Efficiency = NA)

for (i in seq_along(core_counts)) {
  nc <- core_counts[i]
  cat("\nTesting with", nc, "cores...\n")
  
  # Generate fresh response for each replication (using the same seeds)
  # We'll generate y inside the loop to ensure each replication gets its own y.
  # However, to keep timing consistent, we should generate y inside run_one_replication.
  # We'll adapt the run_one_replication to generate y on the fly.
  
  # Revised run function that generates y inside
  run_one_replication_with_y <- function(k) {
    set.seed(seeds[k])
    # Generate y from ZIPXG using the precomputed Theta_true
    Theta_true <- theta_from_mu(Mu_true)
    y <- rZIPX(n, theta = Theta_true, omega = Omega_true)
    if (sum(y == 0) == 0 || sum(y == 0) == n) return(NA)
    
    # Rest of the fitting (same as before)
    n_train <- floor(n * 0.8)
    train_idx <- sample(1:n, n_train, replace = FALSE)
    test_idx <- setdiff(1:n, train_idx)
    
    X_train <- X[train_idx, , drop = FALSE]
    X_test  <- X[test_idx, , drop = FALSE]
    Z_train <- Z[train_idx, , drop = FALSE]
    Z_test  <- Z[test_idx, , drop = FALSE]
    y_train <- y[train_idx]
    Mean_True_test <- Mean_True[test_idx]
    
    df_train <- data.frame(y = y_train, X_train, Z_train)
    form <- as.formula(paste("y ~", paste(colnames(X_train), collapse = "+"), "|", 
                             paste(colnames(Z_train), collapse = "+")))
    
    res_zipx <- fit_zipx_hybrid(X_train, Z_train, y_train)
    if (is.null(res_zipx)) return(NA)
    
    X_test_des <- cbind(Intercept = 1, X_test)
    Z_test_des <- cbind(Intercept = 1, Z_test)
    ncol_X <- ncol(X_test_des)
    beta_hat  <- res_zipx$par[1:ncol_X]
    gamma_hat <- res_zipx$par[(ncol_X + 1):length(res_zipx$par)]
    pred_zipx <- as.vector((1 - plogis(Z_test_des %*% gamma_hat)) * exp(X_test_des %*% beta_hat))
    return(calc_mse_oos_robust(pred_zipx, Mean_True_test))
  }
  
  # Set up parallel backend
  if (nc == 1) {
    # For serial, we use a simple for loop (to avoid overhead of foreach)
    start_time <- Sys.time()
    mse_vals <- numeric(R_bench)
    for (k in 1:R_bench) {
      mse_vals[k] <- run_one_replication_with_y(k)
    }
    end_time <- Sys.time()
  } else {
    cl <- makeCluster(nc, outfile = "")
    registerDoParallel(cl)
    start_time <- Sys.time()
    mse_vals <- foreach(k = 1:R_bench, .combine = c, 
                        .packages = c("MASS", "pscl", "GA")) %dopar% {
                          # Load required libraries in each worker
                          library(MASS)
                          library(pscl)
                          library(GA)
                          run_one_replication_with_y(k)
                        }
    end_time <- Sys.time()
    stopCluster(cl)
  }
  
  time_elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  results$Runtime[i] <- time_elapsed
  cat("  Runtime:", round(time_elapsed, 2), "seconds\n")
}

# Compute speedup and efficiency relative to 1 core
results$Speedup <- results$Runtime[1] / results$Runtime
results$Efficiency <- (results$Speedup / results$Cores) * 100

# ------------------------------------------------------------------------------
# Print results and generate LaTeX table
# ------------------------------------------------------------------------------

cat("\n\n========================================\n")
cat("PERFORMANCE BENCHMARK RESULTS\n")
cat("========================================\n")
print(results)

# Generate LaTeX table
cat("\n\nLaTeX table for manuscript:\n\n")
cat("\\begin{table}[ht]\n")
cat("\\centering\n")
cat("\\caption{Computational performance of the parallel simulation framework.\n")
cat("Results are based on $R =", R_bench, "$ replications of a representative scenario\n")
cat("($n = 500$, $p = 8$, $\\rho = 0.3$, Scenario B). Hardware: MacBook Pro\n")
cat("with Apple M4 Pro processor (12 physical cores) and 24 GB RAM.}\n")
cat("\\label{tab:performance}\n")
cat("\\begin{tabular}{c|c|c|c}\n")
cat("\\toprule\n")
cat("\\textbf{Cores} & \\textbf{Runtime (s)} & \\textbf{Speedup} & \\textbf{Efficiency (\\%)} \\\\\n")
cat("\\midrule\n")
for (i in 1:nrow(results)) {
  cat(sprintf("%d & %.1f & %.2f & %.1f \\\\\n", 
              results$Cores[i], results$Runtime[i], results$Speedup[i], results$Efficiency[i]))
}
cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\end{table}\n")

# ------------------------------------------------------------------------------
# Save results to CSV for reference
# ------------------------------------------------------------------------------
write.csv(results, "performance_benchmark_results.csv", row.names = FALSE)
cat("\nResults saved to performance_benchmark_results.csv\n")