# Required packages
library(dplyr)
library(lubridate)
library(zoo)
library(DPmix)    # For DPM
library(kernlab)  # For GP regression
library(igraph)   # For network sampling
library(MASS)     # For mvrnorm
library(ks)       # For K-S test

# Load FRED data
gdp <- read.csv("GDPC1.csv")
emp <- read.csv("PAYEMS.csv")
inv <- read.csv("GPDI.csv")

# Convert dates
gdp$observation_date <- as.Date(gdp$observation_date)
emp$observation_date <- as.Date(emp$observation_date)
inv$observation_date <- as.Date(inv$observation_date)

# Filter to 1947-2024
gdp <- gdp %>% filter(observation_date >= as.Date("1947-01-01") & 
                      observation_date <= as.Date("2024-12-31"))
emp <- emp %>% filter(observation_date >= as.Date("1947-01-01") & 
                      observation_date <= as.Date("2024-12-31"))
inv <- inv %>% filter(observation_date >= as.Date("1947-01-01") & 
                      observation_date <= as.Date("2024-12-31"))

# Aggregate employment to quarterly
emp$quarter <- as.yearqtr(emp$observation_date)
emp_quarterly <- emp %>%
  group_by(quarter) %>%
  summarise(PAYEMS = mean(PAYEMS, na.rm = TRUE)) %>%
  mutate(observation_date = as.Date(quarter)) %>%
  select(observation_date, PAYEMS)

# Merge data (trim to T=100)
data <- gdp %>%
  inner_join(emp_quarterly, by = "observation_date") %>%
  inner_join(inv, by = "observation_date") %>%
  select(observation_date, GDPC1, PAYEMS, GPDI) %>%
  slice_head(n = 100)

# Standardize data
data$GDPC1 <- scale(data$GDPC1)
data$PAYEMS <- scale(data$PAYEMS)
data$GPDI <- scale(data$GPDI)

# Prepare Y array (N=1000, T=100, d=3)
set.seed(123)
N <- 1000
T <- 100
d <- 3
Y <- array(NA, dim = c(N, T, d))
for (i in 1:N) {
  Y[i,,1] <- data$GDPC1[1:T] + rnorm(T, 0, 0.1)
  Y[i,,2] <- data$PAYEMS[1:T] + rnorm(T, 0, 0.1)
  Y[i,,3] <- data$GPDI[1:T] + rnorm(T, 0, 0.1)
}

# Initialize adjacency matrix (5% connectivity)
A_init <- array(0, dim = c(N, N, T))
for (t in 1:T) {
  A_init[,,t] <- matrix(rbinom(N * N, 1, 0.05), N, N)
  diag(A_init[,,t]) <- 0
}

# Gibbs sampling for ANBEM
anbem_gibbs <- function(Y, A_init, alpha, G0, n_iter = 1000, sigma2 = 0.1) {
  N <- dim(Y)[1]
  T <- dim(Y)[2]
  d <- dim(Y)[3]
  A <- A_init
  theta <- matrix(0, N, d)
  G <- list(clusters = rep(1, N), weights = 1, means = matrix(0, 1, d))
  
  # Initialize theta
  for (i in 1:N) theta[i,] <- rnorm(d, G0$mean, G0$sd)
  
  start_time <- Sys.time()
  
  for (iter in 1:n_iter) {
    # Step 1: Update theta_i
    for (i in 1:N) {
      neighbors <- which(A[i,,1] == 1)
      Y_neighbors <- if (length(neighbors) > 0) Y[neighbors,,] else array(0, c(1,T,d))
      psi_t <- gp_regression(Y_neighbors, A[i,,1])
      mu_i <- psi_t(Y[i,,]) + theta[i,]
      theta[i,] <- mvrnorm(1, mu_i, sigma2 * diag(d))
    }
    
    # Step 2: Update G with DPmix
    dpm_fit <- dp_mix(theta, alpha = alpha, base_mean = G0$mean, base_sd = G0$sd)
    G <- list(clusters = dpm_fit$clusters, weights = dpm_fit$weights, means = dpm_fit$means)
    
    # Step 3: Update A_t
    for (t in 2:T) {
      deg <- colSums(A[,,t-1])
      probs <- deg + 0.1 * rowMeans(Y[,,t-1])
      probs <- probs / sum(probs)
      A[,,t] <- matrix(rbinom(N * N, 1, probs), N, N)
      diag(A[,,t]) <- 0
    }
    
    # Step 4: Update psi_t
    psi_t_list <- list()
    for (t in 1:T) {
      psi_t_list[[t]] <- gausspr(Y[,,t], kernel = "rbfdot", kpar = list(sigma = 1))
    }
  }
  
  comp_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  return(list(theta = theta, G = G, A = A, psi_t = psi_t_list, comp_time = comp_time))
}

# Gaussian Process regression
gp_regression <- function(Y_neighbors, A_i) {
  if (dim(Y_neighbors)[1] == 0) return(function(x) rep(0, dim(x)[3]))
  model <- gausspr(Y_neighbors, kernel = "rbfdot", kpar = list(sigma = 1))
  return(function(x) predict(model, x))
}

# Example usage
G0 <- list(mean = rep(0, d), sd = rep(1, d))
alpha <- 1
result <- anbem_gibbs(Y, A_init, alpha, G0, n_iter = 100)

# Compute MSE and K-S distance
Y_pred <- array(NA, dim = dim(Y))
for (i in 1:N) {
  for (t in 2:T) {
    neighbors <- which(result$A[i,,t-1] == 1)
    Y_pred[i,t,] <- mean(Y[neighbors,t-1,]) + rnorm(d, 0, 0.1)
  }
}
mse <- mean((Y[,,] - Y_pred[,,])^2, na.rm = TRUE)
ks_dist <- ks.test(result$theta[,1], rnorm(1000, 0, 1))$statistic

# Robustness
mse_break <- mean((Y[,80:T,] - Y_pred[,80:T,])^2, na.rm = TRUE)

# Output
cat("MSE:", mse, "\n")
cat("K-S Distance:", ks_dist, "\n")
cat("Computational Time (minutes):", result$comp_time, "\n")
cat("Robustness MSE:", mse_break, "\n")
\end{verbatim}