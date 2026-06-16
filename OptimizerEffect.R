# ==============================================================================
# Disentangle optimizer effect: Compare ZIPXG vs ZINB with GA+BFGS for both
# Scenario: high overdispersion (Scenario B), n=500, p=8, rho=0.3
# ==============================================================================

rm(list = ls())
library(MASS)
library(pscl)
library(GA)

# ------------------------------------------------------------------------------
# Core functions 
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

# Hybrid fit for ZIPXG 
fit_zipx_hybrid <- function(X_raw, Z_raw, y) {
  X_des <- cbind(Intercept = 1, X_raw)
  Z_des <- cbind(Intercept = 1, Z_raw)
  n_par <- ncol(X_des) + ncol(Z_des)
  ga_fitness <- function(p) -neg_logLik_ZIPX(p, X_des, Z_des, y)
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
  }, error = function(e) NULL)
  if(is.null(opt)) return(NULL)
  logLik_val <- -opt$value
  list(par = opt$par, logLik = logLik_val,
       AIC = 2*n_par - 2*logLik_val,
       BIC = n_par*log(length(y)) - 2*logLik_val)
}

# ZINB negative log-likelihood (parametrisation: mu, theta (dispersion), omega)
# NegBin pmf: Gamma(y+theta)/(Gamma(theta)*y!) * (theta/(mu+theta))^theta * (mu/(mu+theta))^y
neg_logLik_ZINB <- function(params, X, Z, y) {
  ncol_X <- ncol(X); ncol_Z <- ncol(Z)
  beta  <- params[1:ncol_X]
  gamma <- params[(ncol_X + 1):(ncol_X + ncol_Z)]
  log_theta <- params[ncol_X + ncol_Z + 1]   # extra parameter for dispersion
  theta <- exp(log_theta)  # ensure positivity
  eta_mu <- as.vector(X %*% beta)
  eta_omega <- as.vector(Z %*% gamma)
  mu <- exp(eta_mu)
  omega <- plogis(eta_omega)
  # Negative binomial likelihood for non-zero part
  log_nb <- function(y, mu, theta) {
    # returns log(P(Y=y)) for NB(mu, theta) where Var = mu + mu^2/theta
    # Using standard formula
    log_choose <- lgamma(y + theta) - lgamma(theta) - lgamma(y + 1)
    log_choose + theta * log(theta/(mu+theta)) + y * log(mu/(mu+theta))
  }
  # Likelihood for zero-inflated NB
  lik <- numeric(length(y))
  for(i in seq_along(y)) {
    if(y[i] == 0) {
      lik[i] <- omega[i] + (1 - omega[i]) * ( (theta/(mu[i]+theta))^theta )
    } else {
      lik[i] <- (1 - omega[i]) * exp(log_nb(y[i], mu[i], theta))
    }
  }
  lik <- pmax(lik, 1e-100)
  -sum(log(lik))
}

# Hybrid fit for ZINB using GA+BFGS
fit_zinb_hybrid <- function(X_raw, Z_raw, y) {
  X_des <- cbind(Intercept = 1, X_raw)
  Z_des <- cbind(Intercept = 1, Z_raw)
  n_par <- ncol(X_des) + ncol(Z_des) + 1   # +1 for log_theta
  # Parameter bounds: beta, gamma, log_theta
  lower <- rep(-5, n_par)
  upper <- rep(5, n_par)
  # log_theta can be wider, e.g., -10 to 10
  lower[length(lower)] <- -10
  upper[length(upper)] <- 10
  ga_fitness <- function(p) -neg_logLik_ZINB(p, X_des, Z_des, y)
  ga_res <- ga(type = "real-valued",
               fitness = ga_fitness,
               lower = lower,
               upper = upper,
               popSize = 20,
               maxiter = 15,
               monitor = FALSE,
               keepBest = TRUE)
  start_val <- ga_res@solution[1, ]
  opt <- tryCatch({
    optim(par = start_val, fn = neg_logLik_ZINB, X = X_des, Z = Z_des, y = y,
          method = "BFGS", control = list(maxit = 1000))
  }, error = function(e) NULL)
  if(is.null(opt)) return(NULL)
  logLik_val <- -opt$value
  n_obs <- length(y)
  n_par_total <- n_par
  list(par = opt$par, logLik = logLik_val,
       AIC = 2*n_par_total - 2*logLik_val,
       BIC = n_par_total*log(n_obs) - 2*logLik_val)
}

# ------------------------------------------------------------------------------
# Simulation parameters (high overdispersion scenario)
# ------------------------------------------------------------------------------
set.seed(123)
n <- 500
p <- 8
rho <- 0.3
R_small <- 1000   # number of replications 

# True coefficients (Scenario B)
beta0 <- 1.5
beta_slopes <- rep(0.3, p)
gamma0 <- -0.4
gamma_slopes <- rep(0.2, p)
beta_true <- matrix(c(beta0, beta_slopes), ncol = 1)
gamma_true <- matrix(c(gamma0, gamma_slopes), ncol = 1)

# Covariance matrix
SIGMA <- matrix(0, p, p)
for(i in 1:p) for(j in 1:p) SIGMA[i,j] <- rho^abs(i-j)

# Storage
mse_zipx <- numeric(R_small)
mse_zinb <- numeric(R_small)

for(r in 1:R_small) {
  if(r %% 10 == 0) cat("Replication", r, "\n")
  set.seed(2025 * 1000 + r)
  X <- mvrnorm(n, rep(0,p), SIGMA)
  Z <- X
  X_des <- cbind(1, X)
  Z_des <- cbind(1, Z)
  mu_true <- exp(X_des %*% beta_true)
  omega_true <- plogis(Z_des %*% gamma_true)
  theta_true <- theta_from_mu(mu_true)
  
  # Generate ZIPXG data
  y <- rZIPX(n, theta = theta_true, omega = omega_true)
  if(sum(y==0)==0 || sum(y==0)==n) next
  
  # Train-test split
  train_idx <- sample(1:n, floor(0.8*n))
  test_idx <- setdiff(1:n, train_idx)
  X_train <- X[train_idx,]; X_test <- X[test_idx,]
  Z_train <- Z[train_idx,]; Z_test <- Z[test_idx,]
  y_train <- y[train_idx]; y_test <- y[test_idx]
  true_mean_test <- (1 - omega_true[test_idx]) * mu_true[test_idx]
  
  # Fit ZIPXG with GA+BFGS
  fit_zipx <- fit_zipx_hybrid(X_train, Z_train, y_train)
  if(!is.null(fit_zipx)) {
    X_test_des <- cbind(1, X_test)
    Z_test_des <- cbind(1, Z_test)
    nX <- ncol(X_test_des)
    beta_hat <- fit_zipx$par[1:nX]
    gamma_hat <- fit_zipx$par[(nX+1):length(fit_zipx$par)]
    pred_zipx <- (1 - plogis(Z_test_des %*% gamma_hat)) * exp(X_test_des %*% beta_hat)
    mse_zipx[r] <- mean((pred_zipx - true_mean_test)^2)
  } else {
    mse_zipx[r] <- NA
  }
  
  # Fit ZINB with GA+BFGS
  fit_zinb <- fit_zinb_hybrid(X_train, Z_train, y_train)
  if(!is.null(fit_zinb)) {
    X_test_des <- cbind(1, X_test)
    Z_test_des <- cbind(1, Z_test)
    nX <- ncol(X_test_des)
    beta_hat <- fit_zinb$par[1:nX]
    gamma_hat <- fit_zinb$par[(nX+1):(nX + ncol(Z_test_des))]
    log_theta <- fit_zinb$par[nX + ncol(Z_test_des) + 1]
    theta_hat <- exp(log_theta)
    mu_pred <- exp(X_test_des %*% beta_hat)
    omega_pred <- plogis(Z_test_des %*% gamma_hat)
    # Conditional mean for ZINB: (1-omega_pred) * mu_pred
    pred_zinb <- (1 - omega_pred) * mu_pred
    mse_zinb[r] <- mean((pred_zinb - true_mean_test)^2)
  } else {
    mse_zinb[r] <- NA
  }
}

# Results
cat("\n===== Disentangling optimizer effect =====\n")
cat("Mean MSE (ZIPXG):", mean(mse_zipx, na.rm=TRUE), "\n")
cat("Mean MSE (ZINB):", mean(mse_zinb, na.rm=TRUE), "\n")
cat("Difference (ZINB - ZIPXG):", mean(mse_zinb, na.rm=TRUE) - mean(mse_zipx, na.rm=TRUE), "\n")