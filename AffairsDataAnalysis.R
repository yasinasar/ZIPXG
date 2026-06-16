# ==============================================================================
# ANALYSIS OF AFFAIRS DATA
# - AIC/BIC/LogLik
# - Coefficient tables with standard errors
# - Out-of-sample MSE
# - Variable selection justification (stepwise AIC)
# ==============================================================================
# Revised: 2026-06-16
# ==============================================================================
rm(list = ls())
# 1. Load Libraries ------------------------------------------------------------
if(!require(glmmTMB)) install.packages("glmmTMB")
if(!require(GA)) install.packages("GA")
if(!require(pscl)) install.packages("pscl")
if(!require(AER)) install.packages("AER")
if(!require(numDeriv)) install.packages("numDeriv")
if(!require(knitr)) install.packages("knitr")
if(!require(kableExtra)) install.packages("kableExtra")
if(!require(dplyr)) install.packages("dplyr")
if(!require(ggplot2)) install.packages("ggplot2")

library(glmmTMB)
library(GA)
library(pscl)
library(AER)
library(numDeriv)
library(knitr)
library(kableExtra)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------------------------
# Core Functions for ZIPXG
# ------------------------------------------------------------------------------
dPX_vec <- function(x, theta) {
  if (any(theta <= 0)) return(rep(0, length(x)))
  log_num <- 2 * log(theta) + log(2 * (1 + theta)^2 + theta * (x + 2) * (x + 1))
  log_den <- log(2) + (x + 4) * log(1 + theta)
  return(exp(log_num - log_den))
}

dPX_safe <- function(x, theta) {
  theta <- pmax(theta, 1e-8)
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

# Hybrid fit for simulation (light settings)
fit_zipx_hybrid <- function(X_raw, Z_raw, y, popSize = 20, maxiter = 15) {
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
               popSize = popSize,
               maxiter = maxiter,
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
# ------------------------------------------------------------------------------
# Helper function to produce a LaTeX table for a standard model (zeroinfl or hurdle)
print_std_model <- function(model, name, formula_obj = final_formula, data = df) {
  if(is.null(model)) {
    cat(paste("Model", name, "failed to converge.\n"))
    return(NULL)
  }
  summ <- summary(model)
  c_mat <- summ$coefficients$count
  z_mat <- summ$coefficients$zero
  # Extract coefficients, SE, etc.
  est <- c(c_mat[,1], z_mat[,1])
  se <- c(c_mat[,2], z_mat[,2])
  z_val <- est / se
  p_val <- 2 * (1 - pnorm(abs(z_val)))
  
  n_count <- nrow(c_mat)
  n_zero <- nrow(z_mat)
  
  # Build data frame
  count_names <- rownames(c_mat)
  zero_names <- rownames(z_mat)
  
  df_count <- data.frame(
    #Component = "Count",
    Covariate = count_names,
    MLE = sprintf("%.4f", est[1:n_count]),
    SE = sprintf("%.4f", se[1:n_count]),
    CI = sprintf("(%.4f, %.4f)", est[1:n_count] - 1.96*se[1:n_count], est[1:n_count] + 1.96*se[1:n_count]),
    z = sprintf("%.4f", z_val[1:n_count]),
    p = sprintf("%.4f", p_val[1:n_count])
  )
  df_zero <- data.frame(
    #Component = "Zero-inflation",
    Covariate = zero_names,
    MLE = sprintf("%.4f", est[(n_count+1):(n_count+n_zero)]),
    SE = sprintf("%.4f", se[(n_count+1):(n_count+n_zero)]),
    CI = sprintf("(%.4f, %.4f)", est[(n_count+1):(n_count+n_zero)] - 1.96*se[(n_count+1):(n_count+n_zero)], 
                 est[(n_count+1):(n_count+n_zero)] + 1.96*se[(n_count+1):(n_count+n_zero)]),
    z = sprintf("%.4f", z_val[(n_count+1):(n_count+n_zero)]),
    p = sprintf("%.4f", p_val[(n_count+1):(n_count+n_zero)])
  )
  df_out <- rbind(df_count, df_zero)
  
  cat(paste0("\n% --- Table for ", name, " ---\n"))
  kbl(df_out, format = "latex", booktabs = TRUE, escape = FALSE, row.names = FALSE,
      col.names = c("Covariate", "MLE", "SE", "95\\% CI", "$z$", "$p$-value"),
      align = "lcccccc", caption = paste("Coefficient estimates for the", name, "model (Affairs data)")) %>%
    kable_styling(latex_options = "hold_position") %>%
    print()
}

# ------------------------------------------------------------------------------
# Load Data and Basic Overdispersion Check
# ------------------------------------------------------------------------------
data("Affairs", package = "AER")
df <- Affairs
y <- df$affairs

cat("===== Overdispersion Check =====\n")
cat("Mean of affairs:", mean(y), "\n")
cat("Variance of affairs:", var(y), "\n")
cat("Variance/Mean ratio:", var(y)/mean(y), "\n\n")

# ------------------------------------------------------------------------------
# Variable Selection 
# ------------------------------------------------------------------------------
# We perform stepwise AIC on a zero-inflated Poisson model (full model) to select
# the final predictors. This is a standard data-driven approach.
# The result is the formula used throughout the analysis.
# Model 1
full_form <- affairs ~ gender + age + yearsmarried + children + religiousness + education + occupation + rating |
  gender + age + yearsmarried + children + religiousness + education + occupation + rating
# ------------------------------------------------------------------------------
# Stepwise selection using `step` on a zeroinfl model (may take a few seconds)
# Note: `step` works with `zeroinfl` objects.
cat("Performing stepwise AIC variable selection (this may take a moment)...\n")
set.seed(562026)
full_model <- zeroinfl(full_form, data = df, dist = "poisson")
step_model <- step(full_model, direction = "backward", trace = TRUE)
final_formula <- formula(step_model)
cat("Final model formula:\n")
print(final_formula)

# Model 2
# best model from step() function with zeroinfl with ZIP:
# final_formula <- affairs ~ gender + age + yearsmarried + children + religiousness + rating |
#                            gender + age + yearsmarried + children + religiousness + rating
# ------------------------------------------------------------------------------
# Model 3
final_formula <- affairs ~ yearsmarried + religiousness + rating |
  religiousness + rating
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Fitting models on the data with selected predictor variables
# ------------------------------------------------------------------------------
m_zip_full  <- zeroinfl(final_formula, data = df, dist = "poisson")
m_zinb_full <- zeroinfl(final_formula, data = df, dist = "negbin")
m_zgeo_full <- tryCatch(zeroinfl(final_formula, data = df, dist = "geometric"), error = function(e) NULL)
m_hp_full   <- hurdle(final_formula, data = df, dist = "poisson")
m_hnb_full  <- hurdle(final_formula, data = df, dist = "negbin")

# Prepare covariate matrices without intercept (for ZIPXG)
# Extract count and zero parts from final_formula (two‑part formula)
f_char <- as.character(final_formula)

# f_char[3] contains the RHS with '|', e.g. "age + yearsmarried + ... | age + ..."
parts <- strsplit(f_char[3], " \\| ")[[1]]
count_part <- as.formula(paste("~", parts[1]))
zero_part  <- as.formula(paste("~", parts[2]))
X_full <- model.matrix(count_part, data = df)[, -1, drop = FALSE]
Z_full <- model.matrix(zero_part, data = df)[, -1, drop = FALSE]
y_full <- df$affairs

fit_zipx_full <- fit_zipx_hybrid(X_full, Z_full, y_full, 
                                 popSize = 50, maxiter = 40)

# ------------------------------------------------------------------------------
# LogLik/AIC/BIC Table
# ------------------------------------------------------------------------------
logLik_zipx <- fit_zipx_full$logLik
aic_zipx <- fit_zipx_full$AIC
bic_zipx <- fit_zipx_full$BIC

cat("\n===== LogLik/AIC/BIC Table =====\n")

tab <- data.frame(
  Model = c("ZIP", "ZINB", "ZGeo", "HP", "HNB", "ZIPXG"),
  LogLik = c(logLik(m_zip_full), logLik(m_zinb_full), 
             if (!is.null(m_zgeo_full)) logLik(m_zgeo_full) else NA, 
             logLik(m_hp_full), logLik(m_hnb_full), logLik_zipx),
  AIC = c(AIC(m_zip_full), AIC(m_zinb_full), 
          if (!is.null(m_zgeo_full)) AIC(m_zgeo_full) else NA, 
          AIC(m_hp_full), AIC(m_hnb_full), aic_zipx),
  BIC = c(BIC(m_zip_full), BIC(m_zinb_full), 
          if (!is.null(m_zgeo_full)) BIC(m_zgeo_full) else NA, 
          BIC(m_hp_full), BIC(m_hnb_full), bic_zipx)
)
print(tab)
# ------------------------------------------------------------------------------
# For ZIPXG, compute bootstrap SEs with B=1000 replications
cat("\n=== Bootstrapping standard errors for ZIPXG (B=1000) ===\n")
B <- 1000
boot_mat <- matrix(NA, nrow = B, ncol = length(fit_zipx_full$par))
set.seed(456)
for (i in 1:B) {
  idx <- sample(1:nrow(df), replace = TRUE)
  X_boot <- X_full[idx, , drop = FALSE]
  Z_boot <- Z_full[idx, , drop = FALSE]
  y_boot <- y[idx]
  fit_boot <- fit_zipx_hybrid(X_boot, Z_boot, y_boot, popSize = 50, maxiter = 40)
  if (!is.null(fit_boot)) {
    boot_mat[i, ] <- fit_boot$par
  }
}
boot_se <- apply(boot_mat, 2, sd, na.rm = TRUE)

# Build ZIPXG coefficient table 
X_names <- colnames(X_full)   # count covariate names (no intercept)
Z_names <- colnames(Z_full)   # zero covariate names (no intercept)

# Add intercepts
count_names <- c("(Intercept)", X_names)   
zero_names  <- c("(Intercept)", Z_names)   

all_names <- c(paste0(count_names), paste0(zero_names))  

est_zipx <- fit_zipx_full$par
boot_se <- apply(boot_mat, 2, sd, na.rm = TRUE)  

z_val_zipx <- est_zipx / boot_se
p_val_zipx <- 2 * (1 - pnorm(abs(z_val_zipx)))

# ------------------------------------------------------------------------------
# Asyptotic standard errors can be computed if the invertibility of the observed Fisher information is possible:
# library(numDeriv)
# # Compute Hessian at the MLE (negative log-likelihood)
# hess <- hessian(function(p) neg_logLik_ZIPX(p, X_full, Z_full, y), fit_zipx_full$par)
# # Observed Fisher information = -Hessian (since we have negative log-likelihood)
# fisher_info <- -hess
# # Asymptotic variance-covariance matrix
# vcov <- solve(fisher_info)
# # Standard errors
# asymp_se <- sqrt(diag(vcov))
# ------------------------------------------------------------------------------
# ZIPXG results
df_zipx <- data.frame(
  #Component = c(rep("Count", length(count_names)), rep("Zero-inflation", length(zero_names))),
  Covariate = all_names,
  MLE = sprintf("%.4f", est_zipx),
  SE = sprintf("%.4f", boot_se),
  CI = sprintf("(%.4f, %.4f)", est_zipx - 1.96*boot_se, est_zipx + 1.96*boot_se),
  z = sprintf("%.4f", z_val_zipx),
  p = sprintf("%.4f", p_val_zipx)
)

# ------------------------------------------------------------------------------
# Coefficients and p-values of the models to see on consol
# ------------------------------------------------------------------------------
cat("\n===== ZIP Model Coefficients =====\n")
print(summary(m_zip_full))

cat("\n===== ZINB Model Coefficients =====\n")
print(summary(m_zinb_full))

cat("\n===== ZGeo Model Coefficients =====\n")
print(summary(m_zgeo_full))

cat("\n===== HP Model Coefficients =====\n")
print(summary(m_hp_full))

cat("\n===== HNB Model Coefficients =====\n")
print(summary(m_hnb_full))

cat("\n===== ZIPXG Model Coefficients =====\n")
print(df_zipx)

# ------------------------------------------------------------------------------
# Out-of-Sample MSE using 1000 replications
# ------------------------------------------------------------------------------
cat("\n===== Out-of-Sample MSE (1000 random 80/20 splits) =====\n")

set.seed(2026)
R_mse <- 1000  # number of random splits (can increase to 500)
mse_results <- data.frame(ZIP = numeric(R_mse), ZINB = numeric(R_mse), 
                          ZGeo = numeric(R_mse), HP = numeric(R_mse), 
                          HNB = numeric(R_mse), ZIPXG = numeric(R_mse),
                          valid_ZIPXG = logical(R_mse))


for (r in 1:R_mse) {
  # Split indices
  n <- nrow(df)
  train_idx <- sample(1:n, floor(0.8 * n))
  test_idx <- setdiff(1:n, train_idx)
  
  X_train <- X_full[train_idx, , drop = FALSE]
  X_test  <- X_full[test_idx, , drop = FALSE]
  Z_train <- Z_full[train_idx, , drop = FALSE]
  Z_test  <- Z_full[test_idx, , drop = FALSE]
  y_train <- y_full[train_idx]
  y_test  <- y_full[test_idx]
  
  df_train <- df[train_idx, ]
  df_test  <- df[test_idx, ]
  
  # Standard models (using formula)
  form <- final_formula
  m_zip  <- zeroinfl(form, data = df_train, dist = "poisson")
  m_zinb <- zeroinfl(form, data = df_train, dist = "negbin")
  m_zgeo <- tryCatch(zeroinfl(form, data = df_train, dist = "geometric"), error = function(e) NULL)
  m_hp   <- hurdle(form, data = df_train, dist = "poisson")
  m_hnb  <- hurdle(form, data = df_train, dist = "negbin")
  
  pred_zip  <- predict(m_zip, newdata = df_test, type = "response")
  pred_zinb <- predict(m_zinb, newdata = df_test, type = "response")
  pred_zgeo <- if (!is.null(m_zgeo)) predict(m_zgeo, newdata = df_test, type = "response") else NA
  pred_hp   <- predict(m_hp, newdata = df_test, type = "response")
  pred_hnb  <- predict(m_hnb, newdata = df_test, type = "response")
  
  # ZIPXG: fit on training
  # Use more intensive GA settings for real data (popSize=50, maxiter=40)
  fit_zipx <- fit_zipx_hybrid(X_train, Z_train, y_train, 
                              popSize = 50, maxiter = 40)
  pred_zipx <- NA
  if (!is.null(fit_zipx)) {
    X_test_des <- cbind(Intercept = 1, X_test)
    Z_test_des <- cbind(Intercept = 1, Z_test)
    nX <- ncol(X_test_des)
    beta_hat  <- fit_zipx$par[1:nX]
    gamma_hat <- fit_zipx$par[(nX+1):length(fit_zipx$par)]
    
    # Compute predictions with safety clamping
    eta_mu <- X_test_des %*% beta_hat
    eta_omega <- Z_test_des %*% gamma_hat
    #eta_mu <- pmin(pmax(eta_mu, -20), 20)        # prevent overflow
    #eta_omega <- pmin(pmax(eta_omega, -20), 20)
    mu_pred <- exp(eta_mu)
    omega_pred <- plogis(eta_omega)
    pred_zipx <- as.vector((1 - omega_pred) * mu_pred)
    # Cap extreme values
    pred_zipx[pred_zipx > 1e6] <- NA
    if (any(is.na(pred_zipx)) || any(is.infinite(pred_zipx))) pred_zipx <- NA
  }
  
  # Compute MSE (only for valid predictions)
  mse_zipx <- if (!is.null(pred_zipx) && all(is.finite(pred_zipx)) && !any(is.na(pred_zipx))) {
    mean((y_test - pred_zipx)^2, na.rm = TRUE)
  } else NA
  
  mse_results[r, ] <- c(
    mean((y_test - pred_zip)^2, na.rm = TRUE),
    mean((y_test - pred_zinb)^2, na.rm = TRUE),
    if (!is.null(pred_zgeo)) mean((y_test - pred_zgeo)^2, na.rm = TRUE) else NA,
    mean((y_test - pred_hp)^2, na.rm = TRUE),
    mean((y_test - pred_hnb)^2, na.rm = TRUE),
    mse_zipx,
    !is.na(mse_zipx)
  )
}

# Average MSE (ignoring invalid ZIPXG replications)
avg_mse <- colMeans(mse_results[, 1:6], na.rm = TRUE)
valid_frac <- mean(mse_results$valid_ZIPXG, na.rm = TRUE)
cat("Proportion of valid ZIPXG replications:", valid_frac, "\n")
cat("Average out-of-sample MSE:\n")
print(round(avg_mse, 4))

# ------------------------------------------------------------------------------
# Generate latex tables ready for paper
# ------------------------------------------------------------------------------
# AIC-BIC Table 
# ------------------------------------------------------------------------------

# tab <- tab[order(tab$AIC), ] # Sort by AIC
zipx_row <- which(tab$Model == "ZIPXG")
# Print LaTeX table
cat("\n===== LaTeX Table for AIC/BIC =====\n")
kbl(tab, format = "latex", booktabs = TRUE, digits = 2,
    row.names = FALSE,
    caption = "Model Selection Criteria: Log-Likelihood, AIC, and BIC for Affairs Data",
    col.names = c("Model", "LogLik", "AIC", "BIC"),
    align = "lccc") %>%
  kable_styling(latex_options = "hold_position", font_size = 10) %>%
  row_spec(zipx_row, bold = TRUE) %>%
  print()

# ------------------------------------------------------------------------------
# Coefficient Tables with Standard Errors
# ------------------------------------------------------------------------------

cat("\n===== Generating Coefficient Tables =====\n")

cat("\n===== Tables for standard models =====\n")
print_std_model(m_zip_full, "ZIP")
print_std_model(m_zinb_full, "ZINB")
print_std_model(m_zgeo_full, "ZGeo")
print_std_model(m_hp_full, "HP")
print_std_model(m_hnb_full, "HNB")

cat("\n===== Table for ZIPXG =====\n")
kbl(df_zipx, format = "latex", booktabs = TRUE, escape = FALSE, row.names = FALSE,
    col.names = c("Covariate", "MLE", "SE", "95\\% CI", "$z$", "$p$-value"),
    align = "lcccccc", caption = "Coefficient estimates for the ZIPXG model (Affairs data, bootstrap SEs)") %>%
  kable_styling(latex_options = "hold_position") %>%
  print()

# ------------------------------------------------------------------------------
cat("\n===== Analysis Complete =====\n")

