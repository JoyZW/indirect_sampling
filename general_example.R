# GWSM Unbiasedness Simulation in R
# Based on Deville & Lavallée (2006) - Property 3.1
# Demonstrates that the Generalized Weight Share Method produces unbiased estimates

# Load required libraries
set.seed(42)

# Simulation parameters
n_simulations <- 1000
N_A <- 100  # Population size for frame population U_A
N_B <- 150  # Population size for target population U_B
n_A <- 30   # Sample size from U_A

cat("GWSM Unbiasedness Simulation\n")
cat("============================\n")
cat("Population sizes: N_A =", N_A, ", N_B =", N_B, "\n")
cat("Sample size: n_A =", n_A, "\n")
cat("Number of simulations:", n_simulations, "\n\n")

# Step 1: Create link matrix Theta_AB
# Each unit in U_A can be linked to 1-3 units in U_B
Theta_AB <- matrix(0, nrow = N_A, ncol = N_B)

for (j in 1:N_A) {
  n_links <- sample(1:3, 1)  # Random number of links (1-3)
  linked_units <- sample(N_B, n_links, replace = FALSE)
  Theta_AB[j, linked_units] <- 1
}

# Ensure every unit in U_B is linked to at least one unit in U_A
for (i in 1:N_B) {
  if (sum(Theta_AB[, i]) == 0) {
    j <- sample(N_A, 1)
    Theta_AB[j, i] <- 1
  }
}

cat("Step 1: Created link matrix Theta_AB\n")
cat("Total links:", sum(Theta_AB > 0), "\n")
cat("All units in U_B linked:", all(colSums(Theta_AB) > 0), "\n\n")

# Step 2: Create standardized link matrix Theta_AB_tilde
# Following equation in paper: Theta_AB_tilde[j,i] = Theta_AB[j,i] / theta_i^+
theta_i_plus <- colSums(Theta_AB)  # Column sums
Theta_AB_tilde <- Theta_AB
for (i in 1:N_B) {
  if (theta_i_plus[i] > 0) {
    Theta_AB_tilde[, i] <- Theta_AB[, i] / theta_i_plus[i]
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
set.seed(42)
Y_B <- pmax(rnorm(N_B, mean = 50, sd = 15), 0)  # Ensure non-negative
true_total_Y_B <- sum(Y_B)

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
  N_B <- ncol(Theta_AB_tilde)
  
  # Create indicator vector t_A
  t_A <- rep(0, N_A)
  t_A[sample_A] <- 1
  
  # Calculate weights: W[i] = sum over j of (Theta_AB_tilde[j,i] * t_A[j] / pi_A[j])
  W <- rep(0, N_B)
  for (i in 1:N_B) {
    for (j in 1:N_A) {
      W[i] <- W[i] + Theta_AB_tilde[j, i] * t_A[j] / pi_A[j]
    }
  }
  
  # Calculate estimate: Y_B_hat = W' * Y_B
  Y_B_hat <- sum(W * Y_B)
  
  return(Y_B_hat)
}

# Step 6: Run simulation
cat("Step 6: Running simulation...\n")
estimates <- numeric(n_simulations)

for (sim in 1:n_simulations) {
  if (sim %% 100 == 0) cat("Progress:", sim, "\n")
  
  # Draw sample from U_A using SRSWOR
  sample_A <- sample(N_A, n_A, replace = FALSE)
  
  # Calculate GWSM estimate
  estimates[sim] <- gwsm_estimator(sample_A, Theta_AB_tilde, pi_A, Y_B)
}

# Step 7: Analyze results
cat("\nRESULTS - GWSM Unbiasedness Verification\n")
cat("==========================================\n")

mean_estimate <- mean(estimates)
std_estimate <- sd(estimates)
bias <- mean_estimate - true_total_Y_B
relative_bias <- bias / true_total_Y_B * 100

cat("True total Y_B:          ", round(true_total_Y_B, 2), "\n")
cat("Mean of estimates:       ", round(mean_estimate, 2), "\n")
cat("Standard deviation:      ", round(std_estimate, 2), "\n")
cat("Bias:                    ", round(bias, 2), "\n")
cat("Relative bias:           ", round(relative_bias, 4), "%\n")

# 95% confidence interval for the mean
ci_margin <- 1.96 * std_estimate / sqrt(n_simulations)
cat("95% CI for mean:         [", round(mean_estimate - ci_margin, 2), 
    ", ", round(mean_estimate + ci_margin, 2), "]\n")

# Statistical test for unbiasedness
t_test <- t.test(estimates, mu = true_total_Y_B)
cat("\nStatistical test for unbiasedness (H0: E[Ŷ_B] = Y_B):\n")
cat("t-statistic:             ", round(t_test$statistic, 4), "\n")
cat("p-value:                 ", round(t_test$p.value, 6), "\n")
cat("Reject H0 at α=0.05:     ", t_test$p.value < 0.05, "\n")

if (abs(relative_bias) < 0.5) {
  cat("\n✓ CONCLUSION: The GWSM estimator is UNBIASED\n")
  cat("  (Relative bias < 0.5%:", round(abs(relative_bias), 4), "%)\n")
} else {
  cat("\n✗ WARNING: Possible bias detected (Relative bias:", 
      round(relative_bias, 4), "%)\n")
}

cat("\nThis simulation demonstrates Property 3.1 from Deville & Lavallée (2006):\n")
cat("'The GWSM provides unbiased estimates if and only if Θ̃_AB is a standardized link matrix.'\n")

# Display some diagnostic information
cat("\nDiagnostic Information:\n")
cat("Range of estimates: [", round(min(estimates), 2), ", ", 
    round(max(estimates), 2), "]\n")
cat("Coefficient of variation:", round(std_estimate/mean_estimate * 100, 2), "%\n")

# Theoretical verification that E[W] = 1_B
cat("\nTheoretical verification:\n")
cat("According to equation (3.1): E[W] = Θ̃_AB' * 1_A\n")
theoretical_E_W <- t(Theta_AB_tilde) %*% ones_A
cat("E[W] equals 1_B:", all(abs(theoretical_E_W - ones_B) < 1e-10), "\n")
cat("This confirms unbiasedness: E[Ŷ_B] = E[W' * Y_B] = 1_B' * Y_B = Y_B\n")
