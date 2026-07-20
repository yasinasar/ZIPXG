# ==============================================================================
# Simulation design 1: for section 4.1
# ==============================================================================

rm(list = ls())

library(doParallel)
library(foreach)
library(ggplot2)
library(gridExtra) # Required to combine the 4 plots

# -----------------------------------------------------------------------------#
# CORE FUNCTIONS
# -----------------------------------------------------------------------------#
dPX <- function(x, theta) {
  num <- theta^2 * (2*(1+theta)^2 + theta*(x+2)*(x+1))
  den <- 2 * (1+theta)^(x+4)
  return(num / den)
}

rPX <- function(n, theta) {
  threshold <- theta / (1 + theta)
  u <- runif(n)
  lambda <- numeric(n)
  idx_exp <- u <= threshold
  n_exp <- sum(idx_exp)
  idx_gamma <- !idx_exp
  n_gamma <- sum(idx_gamma)
  if (n_exp > 0) lambda[idx_exp] <- rexp(n_exp, rate = theta)
  if (n_gamma > 0) lambda[idx_gamma] <- rgamma(n_gamma, shape = 3, rate = theta)
  y <- rpois(n, lambda = lambda)
  return(y)
}

rZIPX <- function(n, theta, omega) {
  is_structural_zero <- rbinom(n, size = 1, prob = omega) == 1
  count_px <- sum(!is_structural_zero)
  result <- numeric(n) 
  if (count_px > 0) {
    px_values <- rPX(n = count_px, theta = theta)
    result[!is_structural_zero] <- px_values
  }
  return(result)
}

logLik_ZIPX_dist <- function(params, x) {
  theta <- params[1]
  omega <- params[2]
  if (theta <= 0) return(-Inf)
  if (omega < 0 || omega >= 1) return(-Inf)
  log_theta <- log(theta)
  log_1_plus_theta <- log(1 + theta)
  inner_term <- 2 * (1 + theta)^2 + theta * (x + 2) * (x + 1)
  log_px_prob <- 2 * log_theta + log(inner_term) - log(2) - (x + 4) * log_1_plus_theta
  px_prob <- exp(log_px_prob)
  probs <- numeric(length(x))
  idx_zero <- (x == 0)
  if (any(idx_zero)) probs[idx_zero] <- omega + (1 - omega) * px_prob[idx_zero]
  idx_pos <- (x > 0)
  if (any(idx_pos)) probs[idx_pos] <- (1 - omega) * px_prob[idx_pos]
  log_lik <- sum(log(probs + 1e-300))
  return(-log_lik)
}

# -----------------------------------------------------------------------------#
# SIMULATION FUNCTION
# -----------------------------------------------------------------------------#

run_simulation_scenario <- function(theta_true, omega_true, nvec, R, cl) {
  
  RESULT.mse <- matrix(NA, nrow = length(nvec), ncol = 2)
  RESULT.bias <- matrix(NA, nrow = length(nvec), ncol = 2)
  RESULT.mre  <- matrix(NA, nrow = length(nvec), ncol = 2)
  
  for (i in 1:length(nvec)) {
    n <- nvec[i]
    
    # Run Parallel Loop
    results_matrix <- foreach(r = 1:R, .combine = rbind, 
                              .export = c("rZIPX", "rPX", "logLik_ZIPX_dist")) %dopar% {
                                simulated_data <- rZIPX(n = n, theta = theta_true, omega = omega_true)
                                tryCatch({
                                  fit <- optim(par = c(0.1, 0.1), fn = logLik_ZIPX_dist, x = simulated_data, 
                                               method = "L-BFGS-B", lower = c(0.001, 0), upper = c(Inf, 0.99))
                                  return(c(fit$par[1], fit$par[2]))
                                }, error = function(e) return(c(NA, NA)))
                              }
    
    results_matrix <- na.omit(results_matrix)
    est_theta <- results_matrix[, 1]
    est_omega <- results_matrix[, 2]
    
    # Calculate Metrics
    RESULT.mse[i, ]  <- c(mean((est_theta - theta_true)^2), mean((est_omega - omega_true)^2))
    RESULT.bias[i, ] <- c(mean(est_theta) - theta_true, mean(est_omega) - omega_true)
    RESULT.mre[i, ]  <- c(mean(abs((est_theta - theta_true) / theta_true)), 
                          mean(abs((est_omega - omega_true) / omega_true)))
  }
  
  # Return data frame for plotting
  df_long <- rbind(
    data.frame(n=nvec, Param="theta", Value=RESULT.mse[,1], Metric="MSE"),
    data.frame(n=nvec, Param="omega", Value=RESULT.mse[,2], Metric="MSE"),
    data.frame(n=nvec, Param="theta", Value=RESULT.bias[,1], Metric="Bias"),
    data.frame(n=nvec, Param="omega", Value=RESULT.bias[,2], Metric="Bias"),
    data.frame(n=nvec, Param="theta", Value=RESULT.mre[,1], Metric="MRE"),
    data.frame(n=nvec, Param="omega", Value=RESULT.mre[,2], Metric="MRE")
  )
  df_long$Metric <- factor(df_long$Metric, levels = c("Bias", "MSE", "MRE"))
  
  return(df_long)
}

# -----------------------------------------------------------------------------#
# EXECUTION
# -----------------------------------------------------------------------------#

# Parallel Setup
num_cores <- parallel::detectCores() - 1 
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat("Running multi-scenario simulation on", num_cores, "cores...\n")

# Global Config
R <- 10000 
nvec <- c(100, 200, 400, 600, 800, 1000)

# Define 4 Scenarios (Modify these values as needed)
scenarios <- list(
  list(theta = 0.9, omega = 0.3), 
  list(theta = 0.9, omega = 0.4),
  list(theta = 1.5, omega = 0.6), 
  list(theta = 2.0, omega = 0.6),  
  list(theta = 0.5, omega = 0.1),
  list(theta = 3.0, omega = 0.7)
)

plot_list <- list()

# Loop through scenarios
for (k in 1:length(scenarios)) {
  
  th <- scenarios[[k]]$theta
  om <- scenarios[[k]]$omega
  
  cat(sprintf("Running Scenario %d: theta=%.1f, omega=%.1f\n", k, th, om))
  
  # Run Simulation
  df_res <- run_simulation_scenario(th, om, nvec, R, cl)
  
  # Create Title Label
  plot_title <- bquote(paste("Scenario ", .(k), ": ", theta == .(th), ", ", omega == .(om)))
  
  # Create Plot
  p <- ggplot(df_res, aes(x = n, y = Value, color = Param, shape = Param)) +
    geom_line(linewidth = 0.8) +        
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
    facet_wrap(~Metric, scales = "free_y", nrow = 1) +
    scale_color_manual(values = c("theta" = "#1f77b4", "omega" = "#d62728"),
                       labels = c("theta" = expression(theta), "omega" = expression(omega))) +
    scale_shape_manual(values = c("theta" = 19, "omega" = 17),
                       labels = c("theta" = expression(theta), "omega" = expression(omega))) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      strip.background = element_rect(fill = "gray95"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(title = plot_title, x = "Sample Size (n)", y = "")
  
  plot_list[[k]] <- p
}

stopCluster(cl)

# -----------------------------------------------------------------------------#
# 4. ARRANGE AND DISPLAY
# -----------------------------------------------------------------------------#

cat("Generating final 2x2 grid plot...\n")

# Combine the 4 plots into one grid
final_plot <- grid.arrange(
  plot_list[[1]], plot_list[[2]], 
  plot_list[[3]], plot_list[[4]],
  plot_list[[5]], plot_list[[6]],
  ncol = 2,
  top = "Simulation Performance: Bias, MSE, and MRE for Different Parameter Sets"
)

# Optional: Save to file
# ggsave("ZIPX_Simulation_Grid.pdf", final_plot, width = 12, height = 8)