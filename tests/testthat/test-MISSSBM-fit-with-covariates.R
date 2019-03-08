context("test missSBM-fit with covariates")

source("utils_test.R")

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

sbm <- simulateSBM(N, alpha, pi, directed, covarMatrix, covarParam)

## control parameter for the VEM
control <- list(threshold = 1e-4, maxIter = 200, fixPointIter = 5, trace = FALSE)

## Consistency
tol_truth <- 1e-2
tol_ARI   <- .9

test_that("missSBM with covariates and dyad sampling works", {

  ## sampled the network
  sampledNet <- sampleNetwork(sbm$adjMatrix, "dyad", covarParam, covarMatrix = covarMatrix)

  ## Perform inference
  missSBM <- missSBM:::missingSBM_fit$new(sampledNet, Q, "dyad", covarMatrix = covarMatrix)
  out <- missSBM$doVEM(control)

  ## Sanity check
  expect_is(missSBM, "missingSBM_fit")
  expect_is(missSBM$fittedSBM, "SBM_fit_covariates")
  expect_is(missSBM$fittedSampling, "dyadSampling_fit_covariates")
  expect_is(missSBM$sampledNetwork, "sampledNetwork")
  expect_equal(out, missSBM$monitoring)

  ## Optimization success
  expect_gte(diff(range(out$objective, na.rm = TRUE)), 0)

  ## SBM: parameters estimation
  expect_lt(error(missSBM$fittedSBM$mixtureParam, sbm$mixtureParam, sort = TRUE), tol_truth)

### FIXME: estimation of gamma is poor...
  expect_lt(error(missSBM$fittedSBM$connectProb, sbm$connectProb), tol_truth*10)

  ## sampling design: parameters estimation
  expect_lt(error(missSBM$fittedSampling$parameters, sbm$covarParam), tol_truth)

  ## clustering
  expect_gt(ARI(missSBM$fittedSBM$memberships, sbm$memberships), tol_ARI)

})

test_that("miss SBM with covariates and node sampling works", {

  ## sampled the network
  sampledNet <- sampleNetwork(sbm$adjMatrix, "node", covarParam, covarMatrix = covarMatrix)

  ## Perform inference
  missSBM <- missSBM:::missingSBM_fit$new(sampledNet, Q, "node", covarMatrix = covarMatrix)
  out <- missSBM$doVEM(control)

  ## Sanity check
  expect_is(missSBM, "missingSBM_fit")
  expect_is(missSBM$fittedSBM, "SBM_fit_covariates")
  expect_is(missSBM$fittedSampling, "nodeSampling_fit_covariates")
  expect_is(missSBM$sampledNetwork, "sampledNetwork")

  ## Optimization success
  expect_gte(diff(range(out$objective, na.rm = TRUE)), 0)

  ## SBM: parameters estimation
  expect_lt(error(missSBM$fittedSBM$mixtureParam, sbm$mixtureParam, sort = TRUE), tol_truth)

### FIXME: estimation of gamma is poor ...
  expect_lt(error(missSBM$fittedSBM$connectProb, sbm$connectProb), tol_truth*10)

  ## sampling design: parameters estimation
  expect_lt(error(missSBM$fittedSampling$parameters, sbm$covarParam), tol_truth)

  ## clustering
  expect_gt(ARI(missSBM$fittedSBM$memberships, sbm$memberships), tol_ARI)

})
