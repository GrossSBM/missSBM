context("test missSBM-fit with covariate")

## ========================================================================
## With covariates

set.seed(178303)
### A SBM model : ###
N <- 300
Q <- 3
alpha <- rep(1,Q)/Q                     # mixture parameter
pi <- diag(.45,Q) + .05                 # connectivity matrix
directed <- FALSE

### Draw a SBM model (Bernoulli, undirected) with covariates
M <- 5
covarMatrix <- matrix(rnorm(N*M,mean = 0, sd = 1), N, M)
covarParam  <- rnorm(M,-.1, .1)

mySBM <- simulateSBM(N, alpha, pi, directed, covarMatrix, covarParam)

test_that("missSBM with covariates and dyad sampling works", {

  ## sampled the network
  sampledNet <- sampleNetwork(mySBM$adjMatrix, "dyad", covarParam, covarMatrix = covarMatrix)

  ## Perform inference
  missSBM <- missSBM:::missingSBM_fit$new(sampledNet, Q, "dyad", covarMatrix = covarMatrix)
  out <- missSBM$doVEM(control)

  ## Sanity check
  expect_is(missSBM, "missingSBM_fit")
  expect_is(missSBM$fittedSBM, "SBM_fit_covariates")
  expect_is(missSBM$fittedSampling, "dyadSampling_fit_covariates")
  expect_is(missSBM$sampledNetwork, "sampledNetwork")
  expect_equal(out, missSBM$monitoring)

  ## Consistency
  tol <- 1e-1

  ## Optimization success
  expect_gte(diff(range(out$objective, na.rm = TRUE)), 0)
  ## SBM: parameters estimation
  expect_lt(sum((missSBM$fittedSBM$connectParam - mySBM$connectParam)^2)/(Q*(Q + 1)/2), tol)
  ## sampling design: parameters estimation
  expect_lt(sum((missSBM$fittedSampling$parameters - covarParam)^2)/M, tol)

  ## clustering
  tol <- .8
  expect_gt(ARI(missSBM$fittedSBM$memberships, mySBM$memberships), tol)

})

test_that("miss SBM with covariates and node sampling works", {

  ## sampled the network
  sampledNet <- sampleNetwork(mySBM$adjMatrix, "node", covarParam, covarMatrix = covarMatrix)

  ## Perform inference
  missSBM <- missSBM:::missingSBM_fit$new(sampledNet, Q, "node", covarMatrix = covarMatrix)
  out <- missSBM$doVEM(control)

  ## Sanity check
  expect_is(missSBM, "missingSBM_fit")
  expect_is(missSBM$fittedSBM, "SBM_fit_covariates")
  expect_is(missSBM$fittedSampling, "nodeSampling_fit_covariates")
  expect_is(missSBM$sampledNetwork, "sampledNetwork")

  ## Consistency
  tol <- 1e-2
  ## Optimization success
  expect_gte(diff(range(out$objective, na.rm = TRUE)), 0)
  ## SBM: parameters estimation
  expect_lt(sum((missSBM$fittedSBM$connectParam - mySBM$connectParam)^2)/(Q*(Q + 1)/2), tol)
  ## sampling design: parameters estimation
  expect_lt(sum((missSBM$fittedSampling$parameters - covarParam)^2)/M, tol)
  ## clustering
  tol <- .8
  expect_gt(ARI(missSBM$fittedSBM$memberships, mySBM$memberships), tol)

})
