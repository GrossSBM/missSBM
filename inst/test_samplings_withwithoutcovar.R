rm(list=ls())
library(missSBM)
library(aricode)

set.seed(1111)

### A SBM model : ###
N <- 200
Q <- 3
alpha <- rep(1,Q)/Q                     # mixture parameter
pi <- diag(.15, Q) + .01                 # connectivity matrix
directed <- FALSE

covarParam  <- c(-5, 5)
M <- length(covarParam)
X <- t(rmultinom(N, 1, c(1/3,1/3)))
covariates <- array(dim = c(N, N, M))
for (i in 1:N)
  for (j in 1:N)
    covariates[i,j,] <- -abs(X[i, ] - X[j, ])

### Draw a undirected SBM model with covariates
mySBM_covar <- simulateSBM(N, alpha, pi, directed, covariates, covarParam)

### Draw a undirected SBM model without covariates
mySBM_nocov <- simulateSBM(N, alpha, pi, directed)

plot(mySBM_covar)
plot(mySBM_nocov)
