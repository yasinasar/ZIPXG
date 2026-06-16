# ==============================================================================
# Timing comparison for a single scenario
# ==============================================================================
library(MASS)
library(pscl)
library(GA)
library(microbenchmark) 

# Use the existing functions (assumed already loaded)
# Source the core functions or copy them here
# ==============================================================================
# Comparison scenario
set.seed(123)
n <- 500; p <- 8; rho <- 0.3
R_timing <- 1000   # number of replications for timing
coefs <- set_overdispersion_coefs(p)
beta.count <- coefs$beta
beta.zero <- coefs$gamma

SIGMA <- matrix(0, p, p)
for(i in 1:p) for(j in 1:p) SIGMA[i,j] <- rho^abs(i-j)

times <- data.frame(ZIP=NA, ZINB=NA, ZGeo=NA, HP=NA, HNB=NA, ZIPXG=NA)

for(r in 1:R_timing) {
  # generate data 
  X <- mvrnorm(n, rep(0,p), SIGMA); Z <- X
  colnames(X) <- colnames(Z) <- paste0("Var",1:p)
  X_des <- cbind(1,X); Z_des <- cbind(1,Z)
  eta_mu <- X_des %*% beta.count
  eta_omega <- Z_des %*% beta.zero
  Mu_true <- exp(eta_mu)
  Omega_true <- plogis(eta_omega)
  Theta_true <- theta_from_mu(Mu_true)
  Mean_True <- (1 - Omega_true) * Mu_true
  y <- rZIPX(n, theta = Theta_true, omega = Omega_true)
  
  # train-test split
  train_idx <- sample(1:n, floor(0.8*n))
  X_train <- X[train_idx,]; X_test <- X[-train_idx,]
  Z_train <- Z[train_idx,]; Z_test <- Z[-train_idx,]
  y_train <- y[train_idx]
  df_train <- data.frame(y=y_train, X_train, Z_train)
  form <- as.formula(paste("y ~", paste(colnames(X_train), collapse = "+"), "|",
                           paste(colnames(Z_train), collapse = "+")))
  
  # time ZIP
  t_zip <- system.time({ m1 <- zeroinfl(form, data=df_train, dist="poisson") })[3]
  # time ZINB
  t_zinb <- system.time({ m2 <- zeroinfl(form, data=df_train, dist="negbin") })[3]
  # time ZGeo
  t_zgeo <- system.time({ m3 <- tryCatch(zeroinfl(form, data=df_train, dist="geometric"), error=function(e) NULL) })[3]
  # time HP
  t_hp <- system.time({ m4 <- hurdle(form, data=df_train, dist="poisson") })[3]
  # time HNB
  t_hnb <- system.time({ m5 <- hurdle(form, data=df_train, dist="negbin") })[3]
  # time ZIPXG
  t_zipxg <- system.time({ res <- fit_zipx_hybrid(X_train, Z_train, y_train) })[3]
  
  times[r,] <- c(t_zip, t_zinb, t_zgeo, t_hp, t_hnb, t_zipxg)
}

# Average times
avg_times <- colMeans(times, na.rm=TRUE)
print(avg_times)