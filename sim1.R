#-------------------------------------------------------------------------------
# Case 1:

# Simple random sampling
# Population v+m
# ------------------------------------------------------------------------------

# GWSM Estimation Function
gwsm_estimator <- function(sample_A, Theta_AB_tilde, pi_A, Y_B) {
  N_A <- nrow(Theta_AB_tilde)
  t_A <- sparseVector(i = sample_A, x = 1, length = N_A)
  weight_vec <- t_A / pi_A
  W <- as.numeric(crossprod(Theta_AB_tilde, weight_vec))
  Y_B_hat <- sum(W * Y_B)
  return(Y_B_hat)
}

# Run simulation
# Initialize storage vectors
simulations <- seq(10, n_simulations, by = 10)
mean_estimates <- numeric(length(simulations))

# Run simulation and store running means
cat("Step 6: Running simulation...\n")
estimates <- numeric(n_simulations)
sampled_A_size <- numeric(n_simulations)
sampled_B_size <- numeric(n_simulations)

for (sim in 1:n_simulations) {
  if (sim %% 100 == 0) cat("Progress:", sim, "\n")
  # Draw sample from U_A using SRSWOR
  sample_A <- sample(N_A, n_A, replace = FALSE)
  sampled_A_size[sim] <- length(sample_A)
  
  linked_B <- which(colSums(Theta_AB[sample_A, , drop = FALSE]) > 0)
  sampled_B_size[sim] <- length(linked_B)
  
  # Calculate GWSM estimate
  estimates[sim] <- gwsm_estimator(sample_A, Theta_AB_tilde, pi_A, Y_B)
  # Store running mean every 10th simulation
  if (sim %% 10 == 0) {
    mean_estimates[sim/10] <- mean(estimates[1:sim])
  }
}

# Plot results
data_plot <- data.frame(
  Simulation = simulations,
  Mean_Estimate = mean_estimates
)

plt <- ggplot(data_plot, aes(x = Simulation, y = Mean_Estimate)) +
  geom_line(color = "blue", linetype = "dashed", size = 1) +
  geom_hline(yintercept = true_total_Y_B, color = "green", linetype = "solid", size = 1) +
  labs(title = "Sim 1: Population v+m | Mean Estimate of Y_B over Simulations",
       x = "Simulation Number",
       y = "Mean Estimate of Y_B") +
  theme_minimal() +
  theme(text = element_text(size = 14))

print(plt)
ggsave(file.path("output", "plt1.png")) 

# Analyze results
cat("\nRESULTS - GWSM Unbiasedness Verification\n")
cat("Empirical ==========================================\n")

mean_estimate <- mean(estimates)
std_estimate <- sd(estimates)
var_estimate <- var(estimates)
bias <- mean_estimate - true_total_Y_B

cat("True total Y_B:          ", round(true_total_Y_B, 2), "\n")
cat("Mean of estimates:       ", round(mean_estimate, 2), "\n")
cat("Var of estimates:       ", round(var_estimate, 2), "\n")
cat("Standard deviation:      ", round(std_estimate, 2), "\n")
cat("Bias:                    ", round(bias, 2), "\n")
cat("Average sample size of A: ", mean(sampled_A_size))
cat("Average sample size of B: ", mean(sampled_B_size))

cat("Theoretical ========================================\n")
cat("\nTheoretical verification:\n")
cat("According to equation (3.1): E[W] = Θ̃_AB' * 1_A\n")
theoretical_E_W <- t(Theta_AB_tilde) %*% ones_A
cat("E[W] equals 1_B:", all(abs(theoretical_E_W - ones_B) < 1e-10), "\n")