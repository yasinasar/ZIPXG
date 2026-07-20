# ==============================================================================
# SIMULATION FOR SCENARIO C: ZINB-AS-DGP (MISSPECIFICATION)
# ==============================================================================
# This script replicates the main simulation but generates data from a ZINB
# distribution instead of ZIPXG. All other settings (coefficients, sample sizes,
# train/test split) are identical to Scenario B.
# ==============================================================================

rm(list = ls())

library(MASS)
library(pscl)
library(doParallel)
library(foreach)
library(knitr)
library(kableExtra)
library(dplyr)
library(GA)

# ------------------------------------------------------------------------------
# Core Functions (unchanged)
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# NEW: ZINB Data Generator (for misspecification)
# ------------------------------------------------------------------------------
rZINB <- function(n, mu, omega, theta_NB = 0.5) {
  # mu: vector of NB means (positive)
  # omega: vector of zero-inflation probabilities
  # theta_NB: dispersion parameter (size) of NB (constant)
  is_structural_zero <- rbinom(n, size = 1, prob = omega) == 1
  n_nb <- sum(!is_structural_zero)
  y <- numeric(n)
  if (n_nb > 0) {
    mu_sub <- mu[!is_structural_zero]
    y[!is_structural_zero] <- rnbinom(n_nb, size = theta_NB, mu = mu_sub)
  }
  return(y)
}

# ------------------------------------------------------------------------------
# Helper functions (same as before)
# ------------------------------------------------------------------------------
calc_mse_oos_robust <- function(pred_vals, true_vals) {
  if(any(is.na(pred_vals)) || any(pred_vals > 10000)) return(NA)
  return(mean((pred_vals - true_vals)^2, na.rm = TRUE))
}

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

# ------------------------------------------------------------------------------
# Simulation parameters (use fewer replications for quick check; set R=10000 for final)
# ------------------------------------------------------------------------------
R <- 1000          # Number of replications (increase to 10000 for final)
n_vec <- c(100, 250, 500, 1000)
p_vec <- c(4, 8, 12)
rho_vec <- c(0.3, 0.5)

# Choose DGP: "ZIPXG" or "ZINB"
dgp <- "ZINB"   # Set to "ZINB" for misspecification scenario

scenarios <- expand.grid(n = n_vec, p = p_vec, rho = rho_vec)
final_results_list <- list()

no_cores <- detectCores() - 1 
cl <- makeCluster(no_cores, outfile = "") 
registerDoParallel(cl)

# ------------------------------------------------------------------------------
# Main simulation loop
# ------------------------------------------------------------------------------
for (s in 1:nrow(scenarios)) {
  curr_n <- scenarios$n[s]
  curr_p <- scenarios$p[s]
  curr_rho <- scenarios$rho[s]
  cat(sprintf("\n[Scenario %d/%d] n=%d, p=%d, rho=%.1f (DGP: %s)\n", 
              s, nrow(scenarios), curr_n, curr_p, curr_rho, dgp))
  
  coefs <- set_overdispersion_coefs(curr_p)
  beta.count <- coefs$beta
  beta.zero  <- coefs$gamma
  
  SIGMA <- matrix(0, curr_p, curr_p)
  for (i in 1:curr_p) {
    for (j in 1:curr_p) {
      SIGMA[i, j] <- curr_rho^abs(i - j)
    }
  }
  
  sim_res <- foreach(k = 1:R, .combine = rbind, 
                     .packages = c("MASS", "pscl", "GA")) %dopar% {
                       set.seed(s * 10000 + k) 
                       
                       # Generate covariates
                       X <- mvrnorm(curr_n, rep(0, curr_p), SIGMA)
                       Z <- X 
                       colnames(X) <- colnames(Z) <- paste0("Var", 1:curr_p)
                       X_des <- cbind(1, X)
                       Z_des <- cbind(1, Z)
                       
                       # True parameters
                       eta_mu <- as.vector(X_des %*% beta.count)
                       eta_omega <- as.vector(Z_des %*% beta.zero)
                       Mu_true <- exp(eta_mu)
                       Omega_true <- plogis(eta_omega)
                       Mean_True <- (1 - Omega_true) * Mu_true
                       
                       # Generate response based on chosen DGP
                       if (dgp == "ZIPXG") {
                         Theta_true <- theta_from_mu(Mu_true)
                         y <- rZIPX(curr_n, theta = Theta_true, omega = Omega_true)
                       } else if (dgp == "ZINB") {
                         theta_NB <- 0.5   # Fixed dispersion parameter (can be adjusted)
                         y <- rZINB(curr_n, mu = Mu_true, omega = Omega_true, theta_NB = theta_NB)
                       } else {
                         stop("Unknown DGP")
                       }
                       
                       # Skip if all zeros or no zeros
                       if(sum(y == 0) == 0 || sum(y == 0) == curr_n) return(rep(NA, 30))
                       
                       # Train-test split
                       n_train <- floor(curr_n * 0.8)
                       train_idx <- sample(1:curr_n, n_train, replace = FALSE)
                       test_idx <- setdiff(1:curr_n, train_idx)
                       
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
                       
                       out <- tryCatch({
                         # Fit models on training data
                         m1 <- zeroinfl(form, data = df_train, dist = "poisson")
                         pred_m1 <- predict(m1, newdata = data.frame(X_test, Z_test), type = "response")
                         mse_oos_m1 <- calc_mse_oos_robust(pred_m1, Mean_True_test)
                         
                         m2 <- zeroinfl(form, data = df_train, dist = "negbin")
                         pred_m2 <- predict(m2, newdata = data.frame(X_test, Z_test), type = "response")
                         mse_oos_m2 <- calc_mse_oos_robust(pred_m2, Mean_True_test)
                         
                         m3 <- tryCatch(zeroinfl(form, data = df_train, dist = "geometric"), error = function(e) NULL)
                         if(!is.null(m3)) {
                           pred_m3 <- predict(m3, newdata = data.frame(X_test, Z_test), type = "response")
                           mse_oos_m3 <- calc_mse_oos_robust(pred_m3, Mean_True_test)
                         } else {
                           mse_oos_m3 <- NA
                         }
                         
                         m4 <- hurdle(form, data = df_train, dist = "poisson")
                         pred_m4 <- predict(m4, newdata = data.frame(X_test, Z_test), type = "response")
                         mse_oos_m4 <- calc_mse_oos_robust(pred_m4, Mean_True_test)
                         
                         m5 <- hurdle(form, data = df_train, dist = "negbin")
                         pred_m5 <- predict(m5, newdata = data.frame(X_test, Z_test), type = "response")
                         mse_oos_m5 <- calc_mse_oos_robust(pred_m5, Mean_True_test)
                         
                         # ZIPXG model
                         res_zipx <- fit_zipx_hybrid(X_train, Z_train, y_train)
                         if(!is.null(res_zipx)) {
                           X_test_des <- cbind(Intercept = 1, X_test)
                           Z_test_des <- cbind(Intercept = 1, Z_test)
                           ncol_X <- ncol(X_test_des)
                           beta_hat  <- res_zipx$par[1:ncol_X]
                           gamma_hat <- res_zipx$par[(ncol_X + 1):length(res_zipx$par)]
                           pred_zipx <- as.vector((1 - plogis(Z_test_des %*% gamma_hat)) * exp(X_test_des %*% beta_hat))
                           mse_oos_zipx <- calc_mse_oos_robust(pred_zipx, Mean_True_test)
                           aic_zipx <- res_zipx$AIC
                           bic_zipx <- res_zipx$BIC
                         } else {
                           mse_oos_zipx <- NA
                           aic_zipx <- NA
                           bic_zipx <- NA
                         }
                         
                         # In-sample AIC/BIC for competitors
                         aic_m1 <- AIC(m1); bic_m1 <- BIC(m1)
                         aic_m2 <- AIC(m2); bic_m2 <- BIC(m2)
                         aic_m3 <- if(!is.null(m3)) AIC(m3) else NA; bic_m3 <- if(!is.null(m3)) BIC(m3) else NA
                         aic_m4 <- AIC(m4); bic_m4 <- BIC(m4)
                         aic_m5 <- AIC(m5); bic_m5 <- BIC(m5)
                         
                         c(aic_m1, aic_m2, aic_m3, aic_m4, aic_m5, aic_zipx,
                           bic_m1, bic_m2, bic_m3, bic_m4, bic_m5, bic_zipx,
                           mse_oos_m1, mse_oos_m2, mse_oos_m3, mse_oos_m4, mse_oos_m5, mse_oos_zipx)
                         
                       }, error = function(e) return(rep(NA, 30)))
                       
                       return(out)
                     }
  
  sim_clean <- as.data.frame(sim_res)
  means <- colMeans(sim_clean, na.rm = TRUE)
  scenario_summary <- c(curr_n, curr_p, curr_rho, nrow(na.omit(sim_clean)), means)
  final_results_list[[s]] <- scenario_summary
}

stopCluster(cl)

# ------------------------------------------------------------------------------
# Process results
# ------------------------------------------------------------------------------
results_df <- as.data.frame(do.call(rbind, final_results_list))
colnames(results_df) <- c("n", "p", "rho", "ValidRuns",
                          paste0("AIC_", c("ZIP","ZINB","ZGeo","HP","HNB","ZIPX")),
                          paste0("BIC_", c("ZIP","ZINB","ZGeo","HP","HNB","ZIPX")),
                          paste0("MSE_OOS_", c("ZIP","ZINB","ZGeo","HP","HNB","ZIPX")))

# ------------------------------------------------------------------------------
# Generate LaTeX Table for Scenario C (Misspecification)
# ------------------------------------------------------------------------------
bold_best <- function(df, columns) {
  df[columns] <- t(apply(df[columns], 1, function(x) {
    best <- min(x, na.rm = TRUE)
    ifelse(x == best, paste0("\\textbf{", sprintf("%.4f", x), "}"), sprintf("%.4f", x))
  }))
  return(df)
}

mse_oos_cols <- 16:21
table_mse <- results_df[, c("n", "p", "rho", 
                            "MSE_OOS_ZIP", "MSE_OOS_ZINB", "MSE_OOS_ZGeo",
                            "MSE_OOS_HP", "MSE_OOS_HNB", "MSE_OOS_ZIPX")]
table_mse <- bold_best(table_mse, 4:9)
colnames(table_mse)[4:9] <- c("ZIP", "ZINB", "ZGeo", "HP", "HNB", "ZIPXG")

latex_code_mse <- kable(table_mse, format = "latex", booktabs = TRUE, escape = FALSE,
                        align = c("c", "c", "c", rep("r", 6)),
                        caption = paste("Out-of-sample MSE under misspecification (ZINB-as-DGP, Scenario C).",
                                        "Best values are highlighted in bold. Results are averaged over $R =", R, "$ replications."),
                        linesep = c("", "", "", "\\addlinespace")) %>%
  add_header_above(c("Parameters" = 3, "Out-of-Sample MSE" = 6)) %>%
  kable_styling(latex_options = c("scale_down", "hold_position"))

cat(latex_code_mse)

# ------------------------------------------------------------------------------
# AIC/BIC table
# ------------------------------------------------------------------------------
aic_cols <- 5:10
bic_cols <- 11:16

table_aic_bic <- results_df[, c("n", "p", "rho", 
                                paste0("AIC_", c("ZIP","ZINB","ZGeo","HP","HNB","ZIPXG")),
                                paste0("BIC_", c("ZIP","ZINB","ZGeo","HP","HNB","ZIPXG")))]
colnames(table_aic_bic)[4:9] <- c("ZIP", "ZINB", "ZGeo", "HP", "HNB", "ZIPXG")
colnames(table_aic_bic)[10:15] <- c("ZIP", "ZINB", "ZGeo", "HP", "HNB", "ZIPXG")

# ------------------------------------------------------------------------------
# Helper function to bold the minimum value in each row
# ------------------------------------------------------------------------------

bold_minimum <- function(x) {
  if (all(is.na(x))) return(rep("NA", length(x)))
  best <- min(x, na.rm = TRUE)
  sapply(x, function(val) {
    if (is.na(val)) return("NA")
    if (abs(val - best) < 1e-8) {
      return(paste0("\\textbf{", sprintf("%.4f", val), "}"))
    } else {
      return(sprintf("%.4f", val))
    }
  })
}


# ==============================================================================
# COMBINED AIC/BIC TABLE WITH BOLDED MINIMUMS 
# ==============================================================================

# Extract AIC columns (5:10) and BIC columns (11:16)
aic_cols <- 5:10
bic_cols <- 11:16

# Create a copy of the numeric data for formatting
combined_numeric <- results_df[, c(1, 2, 3, aic_cols, bic_cols)]

# Apply bolding to AIC columns (positions 4-9)
combined_numeric[, 4:9] <- t(apply(combined_numeric[, 4:9], 1, bold_minimum))

# Apply bolding to BIC columns (positions 10-15)
combined_numeric[, 10:15] <- t(apply(combined_numeric[, 10:15], 1, bold_minimum))

# Set column names
colnames(combined_numeric) <- c("n", "p", "$\\rho$",
                                "ZIP", "ZINB", "ZGeo", "HP", "HNB", "ZIPXG",
                                "ZIP", "ZINB", "ZGeo", "HP", "HNB", "ZIPXG")

# Generate combined LaTeX table
latex_table <- kable(combined_numeric, 
                     format = "latex", 
                     booktabs = TRUE, 
                     escape = FALSE,
                     align = c("c", "c", "c", rep("r", 12)),
                     caption = "Model comparison: AIC and BIC values (descriptive only). Best values per row are highlighted in bold. Formal model comparison should rely on out-of-sample MSE (Table~1). Results are averaged over $R = 10{,}000$ Monte Carlo replications.",
                     linesep = c("", "", "", "\\addlinespace")) %>%
  add_header_above(c("Parameters" = 3, 
                     "AIC (Lower is Better)" = 6,
                     "BIC (Lower is Better)" = 6)) %>%
  kable_styling(latex_options = c("scale_down", "hold_position"))

# Print the table
cat("\n\n")
cat("-" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("- COMBINED AIC/BIC TABLE (with bolded minimums)\n")
cat("-" %>% rep(80) %>% paste(collapse = ""), "\n\n")
cat(latex_table)

# ------------------------------------------------------------------------------
# Also save results to CSV for further use
# ------------------------------------------------------------------------------
write.csv(results_df, "scenario_C_ZINB_results.csv", row.names = FALSE)
cat("\nResults saved to scenario_C_ZINB_results.csv\n")