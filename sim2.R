#-------------------------------------------------------------------------------
# Case 1:

# Simple random sampling
# Population v+m

# TODO
# - plot is weird, fix it
# - runs so slowly due to biased estimator 
# - calculate denominator on each iteration + population
# - add plot for both a) denominator and b) sample mean
# ------------------------------------------------------------------------------
set.seed(seed)

# Simulation parameters
n_simulations <- 1000
N_A <- 3200  # Population size for frame population U_A
N_B <- 2400  # Population size for target population U_B
n_A <- 800   # Sample size from U_A

message("GWSM Unbiasedness Simulation")
message("========")
message("Population sizes: N_A = ", N_A, ", N_B = ", N_B)
message("Sample size: n_A = ", n_A)
message("Number of simulations: ", n_simulations)

# Create link matrix Theta_AB
Theta_AB <- Matrix(0, nrow = N_A, ncol = N_B, sparse = TRUE)
vplusm_min <- 10
vplusm_max <- 30
for (j in 1:N_A) {
  n_links <- sample(vplusm_min:vplusm_max, 1)
  linked_units <- sample(N_B, n_links)
  Theta_AB[j, linked_units] <- 1
}
# Ensure every unit in U_B is linked to at least one unit in U_A
for (i in 1:N_B) {
  # If col sum is empty, i is not linked to any j. In this case
  # randomly select j and set (i,j) to one
  if (sum(Theta_AB[, i]) == 0) {
    j <- sample(N_A, 1)
    Theta_AB[j, i] <- 1
  }
}


theta_i_plus <- colSums(Theta_AB)

# Create standardized link matrix Theta_AB_tilde
# Following equation in paper: Theta_AB_tilde[j,i] = Theta_AB[j,i] / theta_i^+
theta_i_plus <- colSums(Theta_AB)
Theta_AB_tilde <- Theta_AB
for (i in 1:N_B) {
  if (theta_i_plus[i] > 0) {
    Theta_AB_tilde[, i] <- Theta_AB_tilde[, i] / theta_i_plus[i]
  }
}

# Create variable of interest Y for target population U_B
p_y1 <- 0.95
Y_B <- rbinom(N_B, size = 1, prob = p_y1)
true_total_Y_B <- sum(Y_B)

# Sampling design (Simple Random Sampling Without Replacement)
pi_A <- rep(n_A / N_A, N_A)  # Inclusion probabilities

# --- Standard GWSM estimator (population-based)

gwsm_estimator <- function(sample_A, Theta_AB_tilde, pi_A, Y_B) {
  N_A <- nrow(Theta_AB_tilde)
  t_A <- sparseVector(i = sample_A, x = 1, length = N_A)
  weight_vec <- t_A / pi_A
  W <- as.numeric(crossprod(Theta_AB_tilde, weight_vec))
  Y_B_hat <- sum(W * Y_B)
  return(Y_B_hat)
}

# --- Modified GWSM estimator (sample-based)
gwsm_estimator_sample <- function(sample_A, Theta_AB, pi_A, Y_B, N_B) {
  # Subset Theta_AB to sampled units
  Theta_AB_sample <- Theta_AB[sample_A, , drop = FALSE]
  
  # Compute sample-based theta_i_plus (vectorized)
  theta_i_plus_sample <- colSums(Theta_AB_sample)
  
  # Avoid loop: create Theta_AB_tilde_sample using vectorized division
  denom <- theta_i_plus_sample
  denom[denom == 0] <- 1  # For safety to avoid division by zero
  
  # Matrix divide: each column j divided by denom[j]
  Theta_AB_tilde_sample <- sweep(Theta_AB_sample, 2, denom, "/")
  zero_cols <- which(theta_i_plus_sample == 0)
  for (col in zero_cols) {
    Theta_AB_tilde_sample[, col] <- 0
  }
  
  # SRSWOR: weight is 1 / pi_A for all sample_A
  weight <- 1 / pi_A[1]
  W_sample <- as.numeric(crossprod(Theta_AB_tilde_sample, rep(weight, length(sample_A))))
  Y_B_hat_sample <- sum(W_sample * Y_B)
  
  return(list(Y_B_hat_sample = Y_B_hat_sample,
              theta_i_plus_sample = theta_i_plus_sample,
              W_sample = W_sample))
}

# Run simulation
estimates_pop <- numeric(n_simulations)
estimates_sample <- numeric(n_simulations)
theta_sample_sanity <- matrix(0, nrow = n_simulations, ncol = N_B)
sampled_A_size <- numeric(n_simulations)
sampled_B_size <- numeric(n_simulations)

simulations <- seq(10, n_simulations, by = 10)
mean_estimates_pop <- numeric(length(simulations))
mean_estimates_sample <- numeric(length(simulations))

message("Step 6: Running simulation...")
for (sim in 1:n_simulations) {
  if (sim %% 100 == 0) message("Progress: ", sim)
  sample_A <- sample(N_A, n_A, replace = FALSE)
  sampled_A_size[sim] <- length(sample_A)
  linked_B <- which(colSums(Theta_AB[sample_A, , drop = FALSE]) > 0)
  sampled_B_size[sim] <- length(linked_B)
  # GWSM with population theta
  estimates_pop[sim] <- gwsm_estimator(sample_A, Theta_AB_tilde, pi_A, Y_B)
  # GWSM with sample theta
  out <- gwsm_estimator_sample(sample_A, Theta_AB, pi_A, Y_B, N_B)
  estimates_sample[sim] <- out$Y_B_hat_sample
  theta_sample_sanity[sim, ] <- out$theta_i_plus_sample
  
  if (sim %% 10 == 0) {
    mean_estimates_pop[sim/10] <- mean(estimates_pop[1:sim])
    mean_estimates_sample[sim/10] <- mean(estimates_sample[1:sim])
  }
}

# Sanity check for theta_i_plus_sample <= theta_i_plus
all_sane <- all(theta_sample_sanity <= matrix(rep(theta_i_plus, n_simulations),
                                              nrow = n_simulations, byrow = TRUE))
message("Sanity check: All sample-based theta_i_plus <= population theta_i_plus? ", all_sane)

# Results summary
message(sprintf("True total Y_B: %d", true_total_Y_B))
message(sprintf("Mean estimate (population GWSM): %0.2f", mean(estimates_pop)))
message(sprintf("Mean estimate (sample-based GWSM): %0.2f", mean(estimates_sample)))

# Plot results
data_plot2 <- data.frame(
  Simulation = simulations,
  Mean_Estimate_Pop = mean_estimates_pop,
  Mean_Estimate_Sample = mean_estimates_sample
)

# Convert to long format for plotting both mean estimates
data_plot2_long <- data_plot2 %>%
  pivot_longer(cols = c("Mean_Estimate_Pop", "Mean_Estimate_Sample"),
               names_to = "Estimator", values_to = "Mean_Estimate")

plt2 <- ggplot(data_plot2_long, aes(x = Simulation, y = Mean_Estimate, color = Estimator, linetype = Estimator)) +
  geom_line(size = 1) +
  geom_hline(aes(yintercept = true_total_Y_B, color = "True Total Y_B", linetype = "True Total Y_B"), size = 1) +
  scale_color_manual(values = c("Mean_Estimate_Pop" = "blue",
                                "Mean_Estimate_Sample" = "orange",
                                "True Total Y_B" = "green")) +
  scale_linetype_manual(values = c("Mean_Estimate_Pop" = "dashed",
                                   "Mean_Estimate_Sample" = "dashed",
                                   "True Total Y_B" = "solid")) +
  labs(title = "Sim 2: Sample v+m | Mean Estimate of Y_B over Simulations",
       x = "Simulation Number",
       y = "Mean Estimate of Y_B",
       color = "Line",
       linetype = "Line") +
  theme_minimal() +
  theme(text = element_text(size = 14))

print(plt2)

ggsave(file.path("output", "plt2.png")) 