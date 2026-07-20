# ==============================================================================
# Sensitivity analysis: Compare different GA settings
# Scenario: n=500, p=8, rho=0.3, high overdispersion (Scenario B)
# ==============================================================================

# Define GA settings to test
ga_settings <- list(
  light = list(popSize = 20, maxiter = 15),
  moderate = list(popSize = 30, maxiter = 30),
  heavy = list(popSize = 50, maxiter = 50)
)

# Modified hybrid fitting function with adjustable GA parameters
fit_zipx_hybrid_ga <- function(X_raw, Z_raw, y, popSize, maxiter) {
  X_des <- cbind(Intercept = 1, X_raw)
  Z_des <- cbind(Intercept = 1, Z_raw)
  n_par <- ncol(X_des) + ncol(Z_des)
  
  ga_fitness <- function(p) -neg_logLik_ZIPX(p, X_des, Z_des, y)
  
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
  }, error = function(e) NULL)
  
  if(is.null(opt)) return(NULL)
  logLik_val <- -opt$value
  list(par = opt$par, logLik = logLik_val,
       AIC = 2*n_par - 2*logLik_val,
       BIC = n_par*log(length(y)) - 2*logLik_val)
}

# Set up scenario parameters
n <- 500         # sample size
p <- 8           # number of predictors
rho <- 0.3       # degree of correlation
R_sens <- 1000   # number of replications 
set.seed(2025)

# Create high overdispersion coefficients (Scenario B)
coefs <- set_overdispersion_coefs(p)
beta.count <- coefs$beta
beta.zero  <- coefs$gamma

# Covariance matrix
SIGMA <- matrix(0, p, p)
for(i in 1:p) for(j in 1:p) SIGMA[i,j] <- rho^abs(i-j)

# Store results for each GA setting
results_list <- list()

for(set_name in names(ga_settings)) {
  cat("\nRunning GA setting:", set_name, "\n")
  pop <- ga_settings[[set_name]]$popSize
  maxit <- ga_settings[[set_name]]$maxiter
  
  mse_zipx <- numeric(R_sens)
  
  for(r in 1:R_sens) {
    set.seed(2025 * 1000 + r)
    X <- mvrnorm(n, rep(0,p), SIGMA)
    Z <- X
    colnames(X) <- colnames(Z) <- paste0("Var",1:p)
    X_des <- cbind(1,X); Z_des <- cbind(1,Z)
    
    eta_mu <- X_des %*% beta.count
    eta_omega <- Z_des %*% beta.zero
    Mu_true <- exp(eta_mu)
    Omega_true <- plogis(eta_omega)
    Theta_true <- theta_from_mu(Mu_true)
    Mean_true <- (1 - Omega_true) * Mu_true
    
    y <- rZIPX(n, theta = Theta_true, omega = Omega_true)
    if(sum(y==0)==0 || sum(y==0)==n) next
    
    # Train-test split
    train_idx <- sample(1:n, floor(0.8*n))
    test_idx <- setdiff(1:n, train_idx)
    X_train <- X[train_idx,]; X_test <- X[test_idx,]
    Z_train <- Z[train_idx,]; Z_test <- Z[test_idx,]
    y_train <- y[train_idx]; y_test <- y[test_idx]
    Mean_true_test <- Mean_true[test_idx]
    
    # Fit ZIPX with specific GA settings
    fit <- tryCatch(
      fit_zipx_hybrid_ga(X_train, Z_train, y_train, popSize = pop, maxiter = maxit),
      error = function(e) NULL)
    if(is.null(fit)) next
    
    X_test_des <- cbind(1, X_test)
    Z_test_des <- cbind(1, Z_test)
    nX <- ncol(X_test_des)
    beta_hat <- fit$par[1:nX]
    gamma_hat <- fit$par[(nX+1):length(fit$par)]
    pred_zipx <- as.vector((1 - plogis(Z_test_des %*% gamma_hat)) * exp(X_test_des %*% beta_hat))
    mse_zipx[r] <- mean((pred_zipx - Mean_true_test)^2)
  }
  results_list[[set_name]] <- mse_zipx
}

# Summarise
summary_table <- data.frame(
  Setting = names(ga_settings),
  Mean_MSE = sapply(results_list, mean, na.rm=TRUE),
  SD_MSE = sapply(results_list, sd, na.rm=TRUE)
)
print(summary_table)
