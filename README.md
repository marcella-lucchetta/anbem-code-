# anbem-code-
Supplementary materials for the manuscript ‘ANBEM: Flexible Econometric Modeling Through Dynamic Networks’ submitted to JASA
# ANBEM: Flexible Econometric Modeling Through Dynamic Networks

## Supplementary Materials

This repository contains the supplementary materials for the manuscript *"ANBEM: Flexible Econometric Modeling Through Dynamic Networks"* by Marcella Lucchetta, submitted to the *Journal of the American Statistical Association* (JASA) on June 4, 2025.

### Contents
- **`anbem_estimation.R`**: The R code implementing the Gibbs sampling algorithm for ANBEM estimation, as detailed in Appendix A of the manuscript. This code includes data preprocessing, model estimation, and performance evaluation using synthetic and real-world U.S. FRED data (1947Q1–2024Q4).

### Overview
The `anbem_estimation.R` script provides a complete implementation of the Adaptive Network-Based Econometric Model (ANBEM). It performs the following tasks:
- Loads and preprocesses U.S. FRED data (GDPC1, PAYEMS, GPDI) from 1947Q1 to 2024Q4.
- Generates synthetic data for \(N = 1000\) agents over \(T = 100\) periods.
- Implements Gibbs sampling to estimate latent parameters \(\theta_i\), the mixing distribution \(G\), adjacency matrices \(A_t\), and the Gaussian Process kernel \(\psi_t\).
- Computes performance metrics, including mean squared error (MSE) and Kolmogorov-Smirnov (K-S) distances, for synthetic and real data.

### Requirements
To run the code, install the following R packages:
- `dplyr`, `lubridate`, `zoo`: For data manipulation and date handling.
- `DPmix`: For Dirichlet Process Mixture modeling.
- `kernlab`: For Gaussian Process regression.
- `igraph`: For network sampling.
- `MASS`: For multivariate normal sampling.
- `ks`: For K-S tests.

Install them using:
```R
install.packages(c("dplyr", "lubridate", "zoo", "MASS", "ks"))
# Install DPmix, kernlab, and igraph separately if needed
