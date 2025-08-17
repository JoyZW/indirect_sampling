#-------------------------------------------------------------------------------
# TODO:
# - allow for. minor, moderate and severe differences between empirical v+m and population v+m
# - switch from simple random sampling to Poisson sampling
# - implement and assess different scenarios based on section 6 
# (optimal standardized link matrices depends on sampling)
# - make sample sizes consistent with actual data
# - make vcm and vplusm sample values consistent with actual data
# - equation 2 of Prof Li-Chun's note / uncalibrated GWSM estimator using v+m to assess bias

# ------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set up
# Population ratio of consumer to merchants is approximately the same as observed diary data
# Sample size of consumers is approximately the same as diary data
# 95% of merchants have Y = 1, 5% of Y = 0
# ------------------------------------------------------------------------------

# GWSM Unbiasedness Simulation in R
# Based on Deville & Lavallée (2006) - Property 3.1
# Demonstrates that the Generalized Weight Share Method produces unbiased estimates
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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

# Create link matrix Theta_AB
# Each unit in U_A can be linked to vplusm_min-vplusm_max units in U_B
Theta_AB <- Matrix(0, nrow = N_A, ncol = N_B, sparse = TRUE)
vplusm_min <- 10
vplusm_max <- 30
for (j in 1:N_A) {
  n_links <- sample(vplusm_min:vplusm_max, 1)
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

# Create standardized link matrix Theta_AB_tilde
# Following equation in paper: Theta_AB_tilde[j,i] = Theta_AB[j,i] / theta_i^+
theta_i_plus <- colSums(Theta_AB)
Theta_AB_tilde <- Theta_AB
for (i in 1:N_B) {
  if (theta_i_plus[i] > 0) {
    Theta_AB_tilde[, i] <- Theta_AB_tilde[, i] / theta_i_plus[i]
  }
}

cat("Column sums all equal to 1:", all(abs(colSums(Theta_AB_tilde) - 1) < 1e-10), "\n")

# Verify Result 1: Theta_AB_tilde' %*% 1_A = 1_B
ones_A <- rep(1, N_A)
ones_B <- rep(1, N_B)
result1_check <- t(Theta_AB_tilde) %*% ones_A
cat("Result 1 verification: all(Theta_AB_tilde' * 1_A = 1_B):",
    all(abs(result1_check - ones_B) < 1e-10), "\n\n")

# Create variable of interest Y for target population U_B
p_y1 <- 0.95
# Generate binary vector: 1 with probability p, 0 otherwise
Y_B <- rbinom(N_B, size = 1, prob = p_y1)
# Check proportion of ones
prop_one <- mean(Y_B == 1)
cat("Proportion of ones:", prop_one, "\n")
# Sum (true total)
true_total_Y_B <- sum(Y_B)
cat("True total Y_B:", true_total_Y_B, "\n")

# cat("Step 3: Created variable of interest Y_B\n")
# cat("True total Y_B =", round(true_total_Y_B, 2), "\n")
# cat("Average Y_i =", round(mean(Y_B), 2), "\n\n")

# Define sampling design (Simple Random Sampling Without Replacement)
pi_A <- rep(n_A / N_A, N_A)  # Inclusion probabilities

cat("Sampling design defined\n")
cat("Inclusion probability π_j =", round(pi_A[1], 3), "\n\n")

#-------------------------------------------------------------------------------
# Simulations
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------ Sim 1
source("sim1.R")

# ------------------------------------------------------------------------ Sim 2
source("sim1.R")