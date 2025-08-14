#-------------------------------------------------------------------------------
# TODO:
# - switch from simple random sampling to Poisson sampling
# - implement and assess different scenarios based on section 6 
# (optimal standardized link matrices depends on sampling)
# - make sample sizes consistent with actual data
# - is there a way to show consistency? i.e. show that the var goes to 0?
# ------------------------------------------------------------------------------

# GWSM Unbiasedness Simulation in R
# Based on Deville & Lavallée (2006) - Property 3.1
# Demonstrates that the Generalized Weight Share Method produces unbiased estimates
rm(list = ls())

# Load required libraries
set.seed(2025)
library(Matrix)
library(ggplot2)

# Simulation parameters
n_simulations <- 1000
N_A <- 3200  # Population size for frame population U_A
N_B <- 2400  # Population size for target population U_B
n_A <- 800   # Sample size from U_A

cat("GWSM Unbiasedness Simulation\n")
cat("============================\n")
cat("Population sizes: N_A =", N_A, ", N_B =", N_B, "\n")
cat("Sample size: n_A =", n_A, "\n")
cat("Number of simulations:", n_simulations, "\n\n")

# Step 1: Create link matrix Theta_AB
# Each unit in U_A can be linked to 1-3 units in U_B
Theta_AB <- Matrix(0, nrow = N_A, ncol = N_B, sparse = TRUE)

for (j in 1:N_A) {
  n_links <- sample(1:3, 1)
  linked_units <- sample(N_B, n_links)
  Theta_AB[j, linked_units] <- 1
}
# Print dimensions of Theta_AB
cat("Dimensions of Theta_AB:", dim(Theta_AB), "\n")

# Count number of cells with value 0
num_zeros <- sum(Theta_AB == 0)
cat("Number of 0 cells:", num_zeros, "\n")

# Ensure every unit in U_B is linked to at least one unit in U_A
for (i in 1:N_B) {
  # If col sum is empty, i is not linked to any j. In this case
  # randomly select j and set (i,j) to one
  if (sum(Theta_AB[, i]) == 0) {
    j <- sample(N_A, 1)
    Theta_AB[j, i] <- 1
  }
}

# Count number of cells with value 0
num_zeros <- sum(Theta_AB == 0)
cat("Number of 0 cells:", num_zeros, "\n")

cat("Total links:", sum(Theta_AB > 0), "\n")
cat("All units in U_B linked:", all(colSums(Theta_AB) > 0), "\n\n")

# Step 2: Create standardized link matrix Theta_AB_tilde
# Following equation in paper: Theta_AB_tilde[j,i] = Theta_AB[j,i] / theta_i^+
theta_i_plus <- colSums(Theta_AB)
Theta_AB_tilde <- Theta_AB
for (i in 1:N_B) {
  if (theta_i_plus[i] > 0) {
    Theta_AB_tilde[, i] <- Theta_AB_tilde[, i] / theta_i_plus[i]
  }
}

cat("Step 2: Created standardized link matrix\n")
cat("Column sums all equal to 1:", all(abs(colSums(Theta_AB_tilde) - 1) < 1e-10), "\n")

# Verify Result 1: Theta_AB_tilde' %*% 1_A = 1_B
ones_A <- rep(1, N_A)
ones_B <- rep(1, N_B)
result1_check <- t(Theta_AB_tilde) %*% ones_A
cat("Result 1 verification: all(Theta_AB_tilde' * 1_A = 1_B):",
    all(abs(result1_check - ones_B) < 1e-10), "\n\n")

# Step 3: Create variable of interest Y for target population U_B
p_y1 <- 0.95
# Generate binary vector: 1 with probability p, 0 otherwise
Y_B <- rbinom(N_B, size = 1, prob = p_y1)
# Check proportion of ones
prop_one <- mean(Y_B == 1)
cat("Proportion of ones:", prop_one, "\n")
# Sum (true total)
true_total_Y_B <- sum(Y_B)
cat("True total Y_B:", true_total_Y_B, "\n")

cat("Step 3: Created variable of interest Y_B\n")
cat("True total Y_B =", round(true_total_Y_B, 2), "\n")
cat("Average Y_i =", round(mean(Y_B), 2), "\n\n")

# Step 4: Define sampling design (Simple Random Sampling Without Replacement)
pi_A <- rep(n_A / N_A, N_A)  # Inclusion probabilities

cat("Step 4: Sampling design defined\n")
cat("Inclusion probability π_j =", round(pi_A[1], 3), "\n\n")

# Step 5: GWSM Estimation Function
gwsm_estimator <- function(sample_A, Theta_AB_tilde, pi_A, Y_B) {
  N_A <- nrow(Theta_AB_tilde)
  t_A <- sparseVector(i = sample_A, x = 1, length = N_A)
  weight_vec <- t_A / pi_A
  W <- as.numeric(crossprod(Theta_AB_tilde, weight_vec))
  Y_B_hat <- sum(W * Y_B)
  return(Y_B_hat)
}

# Step 6: Run simulation
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

ggplot(data_plot, aes(x = Simulation, y = Mean_Estimate)) +
  geom_line(color = "blue", linetype = "dashed", size = 1) +
  geom_hline(yintercept = true_total_Y_B, color = "green", linetype = "solid", size = 1) +
  labs(title = "Mean Estimate of Y_B over Simulations",
       x = "Simulation Number",
       y = "Mean Estimate of Y_B") +
  theme_minimal() +
  theme(text = element_text(size = 14))

# Step 7: Analyze results
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