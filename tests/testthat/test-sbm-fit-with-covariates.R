context("test-sbm-fit-with-covariates")

library(aricode)
library(blockmodels)

set.seed(178303)
### A SBM model : ###
N <- 300
Q <- 3
alpha <- rep(1,Q)/Q                     # mixture parameter
pi <- diag(.45,Q) + .05                 # connectivity matrix
directed <- FALSE

### Draw a SBM model (Bernoulli, undirected) with covariates
M <- 10
covariates <- matrix(rnorm(N*M,mean = 0, sd = 1), N, M)
covarParam  <- rnorm(M,0,1)

mySBM <- simulateSBM(N, alpha, pi, directed, covariates, covarParam)
A <- mySBM$adjacencyMatrix


test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
