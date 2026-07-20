# ==============================================================================
# VUONG CLOSENESS TEST
# ==============================================================================
library(tidyverse)
library(knitr)
library(kableExtra)
library(countreg)

# Observation-wise log-likelihood for ZIPXG
ll_obs_ZIPXG <- function(params, X_raw, Z_raw, y) {
  # Add intercepts to match fit_zipx_hybrid
  X_des <- cbind(Intercept = 1, X_raw)
  Z_des <- cbind(Intercept = 1, Z_raw)
  
  ncol_X <- ncol(X_des)
  ncol_Z <- ncol(Z_des)
  beta  <- params[1:ncol_X]
  gamma <- params[(ncol_X + 1):(ncol_X + ncol_Z)]
  
  eta_mu    <- as.vector(X_des %*% beta)
  eta_omega <- as.vector(Z_des %*% gamma)
  
  # Safety clamps matching neg_logLik_ZIPX
  #eta_mu <- pmin(pmax(eta_mu, -20), 20)
  #eta_omega <- pmin(pmax(eta_omega, -20), 20)
  
  mu    <- exp(eta_mu)
  omega <- plogis(eta_omega)
  theta <- theta_from_mu(mu)
  
  p_px_y <- dPX_vec(y, theta)   
  p_px_0 <- dPX_vec(0, theta) 
  
  lik_0 <- omega + (1 - omega) * p_px_0
  lik_1 <- (1 - omega) * p_px_y
  lik <- ifelse(y == 0, lik_0, lik_1)
  
  lik <- pmax(lik, 1e-100) # Safety floor
  
  # Return the VECTOR of individual log-likelihoods, not the sum
  return(log(lik))
}

# Observation-wise log-likelihood for pscl models
ll_obs_pscl <- function(model, y) {
  model_type <- class(model)[1]  # "zeroinfl" or "hurdle"
  
  # FIX: Safely extract the distribution type
  dist_type <- if (is.list(model$dist)) model$dist$count else model$dist
  
  # Extract matrices and coefficients
  X_count <- model.matrix(model, model = "count")
  Z_zero  <- model.matrix(model, model = "zero")
  
  beta_count <- coef(model, model = "count")
  gamma_zero <- coef(model, model = "zero")
  
  # Calculate mu and zero probabilities manually
  mu <- exp(as.vector(X_count %*% beta_count))
  zero_pred <- plogis(as.vector(Z_zero %*% gamma_zero))
  
  # Calculate base distribution probabilities
  if (dist_type == "poisson") {
    p0 <- dpois(0, lambda = mu)
    py <- dpois(y, lambda = mu)
  } else if (dist_type %in% c("negbin", "geometric")) {
    theta <- if(dist_type == "geometric") 1 else model$theta
    p0 <- dnbinom(0, size = theta, mu = mu)
    py <- dnbinom(y, size = theta, mu = mu)
  } else {
    stop("Distribution not supported by this extractor.")
  }
  
  # Calculate final likelihoods based on model structure
  if (model_type == "zeroinfl") {
    omega <- zero_pred 
    lik_0 <- omega + (1 - omega) * p0
    lik_y <- (1 - omega) * py
    lik <- ifelse(y == 0, lik_0, lik_y)
    
  } else if (model_type == "hurdle") {
    pi_nonzero <- zero_pred 
    lik_0 <- 1 - pi_nonzero
    lik_y <- pi_nonzero * (py / (1 - p0)) 
    lik <- ifelse(y == 0, lik_0, lik_y)
  }
  
  return(log(pmax(lik, 1e-100)))
}

# Main Vuong Test Function
vuong_test <- function(ll1, ll2, name1 = "Model 1", name2 = "Model 2") {
  n <- length(ll1)
  
  # Calculate difference in observation-wise log-likelihoods
  m <- ll1 - ll2
  
  # Vuong Statistic components
  mean_m <- mean(m)
  sd_m <- sd(m)
  
  # Calculate V-statistic
  vuong_stat <- sqrt(n) * mean_m / sd_m
  
  # Two-tailed p-value
  p_value <- 2 * pnorm(-abs(vuong_stat))
  
  cat("\n====================================================\n")
  cat(" Vuong Non-Nested Hypothesis Test \n")
  cat("====================================================\n")
  cat(sprintf("Model 1: %s\n", name1))
  cat(sprintf("Model 2: %s\n", name2))
  cat("----------------------------------------------------\n")
  cat(sprintf("Test Statistic (V) : %.4f\n", vuong_stat))
  cat(sprintf("p-value            : %.4e\n", p_value))
  cat("----------------------------------------------------\n")
  
  # Decision Rule (alpha = 0.05)
  if (vuong_stat > 1.96) {
    cat(sprintf("Conclusion: [%s] is strictly preferred.\n", name1))
  } else if (vuong_stat < -1.96) {
    cat(sprintf("Conclusion: [%s] is strictly preferred.\n", name2))
  } else {
    cat("Conclusion: No significant difference between models.\n")
  }
  cat("====================================================\n")
  
  return(invisible(list(statistic = vuong_stat, p_value = p_value)))
}

# ------------------------------------------------------------------------------
# Test results
# ------------------------------------------------------------------------------
cat("\n===== Running Vuong Tests =====\n")

# Extract vector for ZIPXG model
ll_zipx <- ll_obs_ZIPXG(params = fit_zipx_full$par, 
                        X_raw = X_full, 
                        Z_raw = Z_full, 
                        y = y_full)

# Extract vectors for competing models
ll_zip <- ll_obs_pscl(m_zip_full, y_full)
ll_zinb <- ll_obs_pscl(m_zinb_full, y_full)
ll_zgeo <- ll_obs_pscl(m_zgeo_full, y_full)
ll_hnb  <- ll_obs_pscl(m_hnb_full, y_full)
ll_hp  <- ll_obs_pscl(m_hp_full, y_full)

# Compare ZIPXG vs ZIP
vuong_test(ll1 = ll_zipx, 
           ll2 = ll_zip, 
           name1 = "ZIPXG", 
           name2 = "ZIP")

# Compare ZIPXG vs ZINB
vuong_test(ll1 = ll_zipx, 
           ll2 = ll_zinb, 
           name1 = "ZIPXG", 
           name2 = "ZINB")

# Compare ZIPXG vs HP
vuong_test(ll1 = ll_zipx, 
           ll2 = ll_hp, 
           name1 = "ZIPXG", 
           name2 = "HP")

# Compare ZIPXG vs HNB
vuong_test(ll1 = ll_zipx, 
           ll2 = ll_hnb, 
           name1 = "ZIPXG", 
           name2 = "HNB")

# Compare ZIPXG vs ZGeo
vuong_test(ll1 = ll_zipx, 
           ll2 = ll_zgeo, 
           name1 = "ZIPXG", 
           name2 = "ZGeo")


cat("\n===== VERIFICATION: LOG-LIKELIHOOD SUMMATION =====\n")
cat(sprintf("ZIPXG Sum : %12.4f | Model LogLik: %12.4f | Diff: %g\n", 
            sum(ll_zipx), fit_zipx_full$logLik, abs(sum(ll_zipx) - fit_zipx_full$logLik)))

cat(sprintf("ZINB Sum  : %12.4f | Model LogLik: %12.4f | Diff: %g\n", 
            sum(ll_zinb), as.numeric(logLik(m_zinb_full)), abs(sum(ll_zinb) - as.numeric(logLik(m_zinb_full)))))

cat(sprintf("HNB Sum   : %12.4f | Model LogLik: %12.4f | Diff: %g\n", 
            sum(ll_hnb), as.numeric(logLik(m_hnb_full)), abs(sum(ll_hnb) - as.numeric(logLik(m_hnb_full)))))

cat(sprintf("ZIP Sum   : %12.4f | Model LogLik: %12.4f | Diff: %g\n", 
            sum(ll_zip), as.numeric(logLik(m_zip_full)), abs(sum(ll_zip) - as.numeric(logLik(m_zip_full)))))

cat(sprintf("HP Sum    : %12.4f | Model LogLik: %12.4f | Diff: %g\n", 
            sum(ll_hp), as.numeric(logLik(m_hp_full)), abs(sum(ll_hp) - as.numeric(logLik(m_hp_full)))))

if (!is.null(m_zgeo_full)) {
  cat(sprintf("ZGeo Sum  : %12.4f | Model LogLik: %12.4f | Diff: %g\n", 
              sum(ll_zgeo), as.numeric(logLik(m_zgeo_full)), abs(sum(ll_zgeo) - as.numeric(logLik(m_zgeo_full)))))
}
cat("====================================================\n")
# ------------------------------------------------------------------------------
# Create a data frame of Vuong test results
vuong_results <- data.frame(
  Comparison = c("ZIPXG vs ZIP", "ZIPXG vs ZINB", "ZIPXG vs ZGeo", "ZIPXG vs HP", "ZIPXG vs HNB"),
  V_Statistic = c(5.1632, 4.0170, 2.0318, 5.1551, 3.5670),
  p_value = c("2.43e-07", "5.90e-05", "0.0422", "2.54e-07", "3.61e-04"),
  Preferred = rep("ZIPXG", 5)
)
# ------------------------------------------------------------------------------
# LaTeX table
kbl(vuong_results, format = "latex", booktabs = TRUE,
    caption = "Vuong test results for non‑nested model comparisons (Affairs data)",
    col.names = c("Comparison", "V statistic", "$p$-value", "Preferred model"),
    align = "lccr") %>%
  kable_styling(latex_options = "hold_position") %>%
  print()


# ------------------------------------------------------------------------------
# Rootograms
# ------------------------------------------------------------------------------
# Compute observed frequencies
obs_freqs <- table(df$affairs)
max_count <- max(df$affairs)
obs_vec <- rep(0, max_count + 1)
obs_vec[as.numeric(names(obs_freqs)) + 1] <- as.vector(obs_freqs)

# Expected frequencies for ZIPXG (manual)
X_des <- cbind(Intercept = 1, X_full)
Z_des <- cbind(Intercept = 1, Z_full)
beta_hat <- fit_zipx_full$par[1:ncol(X_des)]
gamma_hat <- fit_zipx_full$par[(ncol(X_des)+1):length(fit_zipx_full$par)]

mu_zipx <- exp(X_des %*% beta_hat)
omega_zipx <- plogis(Z_des %*% gamma_hat)
theta_zipx <- theta_from_mu(mu_zipx)

# Expected frequencies for ZIPXG 
exp_freqs_zipx <- numeric(max_count + 1)
for (k in 0:max_count) {
  if (k == 0) {
    prob_k <- omega_zipx + (1 - omega_zipx) * dPX_safe(0, theta_zipx)
  } else {
    prob_k <- (1 - omega_zipx) * dPX_safe(k, theta_zipx)
  }
  exp_freqs_zipx[k+1] <- sum(prob_k)
}

# Rootogram function 
manual_rootogram <- function(observed, expected, main = "", 
                             xlab = "Count", ylab = "sqrt(Frequency)",
                             legend_pos = "topright") {
  max_len <- max(length(observed), length(expected))
  observed <- c(observed, rep(0, max_len - length(observed)))
  expected <- c(expected, rep(0, max_len - length(expected)))
  
  obs_sqrt <- sqrt(observed)
  exp_sqrt <- sqrt(expected)
  ylim <- c(0, max(obs_sqrt, exp_sqrt, na.rm = TRUE) * 1.05)
  x_vals <- 0:(max_len - 1)
  
  plot(x_vals, exp_sqrt, type = "b", pch = 16, lty = 2, col = "blue",
       xlab = xlab, ylab = ylab, main = main, ylim = ylim)
  for(i in 1:length(x_vals)) {
    segments(x_vals[i], exp_sqrt[i], x_vals[i], obs_sqrt[i], lwd = 2, col = "red")
    points(x_vals[i], obs_sqrt[i], pch = 18, col = "red")
  }
  legend(legend_pos, legend = c("Expected (sqrt)", "Observed (sqrt)"), 
         col = c("blue", "red"), lty = c(2, 1), pch = c(16, 18), bty = "n")
}

# Helper function to get expected frequencies for standard models
get_fitted_freqs <- function(model, max_count) {
  model_type <- class(model)[1]
  dist_type <- if(is.list(model$dist)) model$dist$count else model$dist
  
  X_count <- model.matrix(model, model = "count")
  Z_zero  <- model.matrix(model, model = "zero")
  beta_count <- coef(model, model = "count")
  gamma_zero <- coef(model, model = "zero")
  
  mu <- exp(as.vector(X_count %*% beta_count))
  zero_pred <- plogis(as.vector(Z_zero %*% gamma_zero))
  
  exp_freqs <- numeric(max_count + 1)
  
  for (k in 0:max_count) {
    if (dist_type == "poisson") {
      p_k <- dpois(k, lambda = mu)
      p_0 <- dpois(0, lambda = mu)
    } else if (dist_type %in% c("negbin", "geometric")) {
      theta <- if(dist_type == "geometric") 1 else model$theta
      p_k <- dnbinom(k, size = theta, mu = mu)
      p_0 <- dnbinom(0, size = theta, mu = mu)
    }
    
    if (model_type == "zeroinfl") {
      prob_k <- if (k == 0) zero_pred + (1 - zero_pred) * p_0 else (1 - zero_pred) * p_k
    } else if (model_type == "hurdle") {
      prob_k <- if (k == 0) (1 - zero_pred) else zero_pred * (p_k / (1 - p_0))
    }
    exp_freqs[k + 1] <- sum(prob_k)
  }
  
  return(exp_freqs)
}

# List of all models
models_list <- list(
  ZIP = m_zip_full,
  ZINB = m_zinb_full,
  ZGeo = m_zgeo_full,
  HP = m_hp_full,
  HNB = m_hnb_full
)

# Initialize EPS graphics device settings
setEPS()
postscript("rootograms_all.eps", width = 9, height = 7, bg = "white")
par(mfrow = c(2, 3), mar = c(4,4,3,2))

# Generate plots for standard models
for(name in names(models_list)) {
  m <- models_list[[name]]
  if(!is.null(m)) {
    exp_freqs <- get_fitted_freqs(m, max_count)
    if(length(exp_freqs) < length(obs_vec))
      exp_freqs <- c(exp_freqs, rep(0, length(obs_vec) - length(exp_freqs)))
    manual_rootogram(obs_vec, exp_freqs, main = name, legend_pos = "topright")
  } else {
    plot(0, type = "n", main = name, xlab = "", ylab = "")
    text(0, 0, "Model not available", col = "red")
  }
}

# Add ZIPXG in the 6th slot
manual_rootogram(obs_vec, exp_freqs_zipx, main = "ZIPXG", legend_pos = "topright")

dev.off()
cat("Rootograms saved to rootograms_all.eps\n")
