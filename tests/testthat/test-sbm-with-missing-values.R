context("test-sbm-with-missing-values")

library(aricode)

set.seed(1890718)
### A SBM model : ###
N <- 400
Q <- 5
alpha <- rep(1,Q)/Q       # mixture parameter
pi <- diag(.45,Q) + .05   # connectivity matrix
directed <- FALSE         # if the network is directed or not

### Draw a SBM model
mySBM <- simulateSBM(N, alpha, pi, directed) # simulation of ad Bernoulli non-directed SBM
A <- mySBM$adjMatrix             # the adjacency matrix

## control parameter for the VEM
control <- list(threshold = 1e-4, maxIter = 200, fixPointIter = 5, trace = FALSE)

test_that("miss SBM with dyad sampling works", {

  ## sampled the network
  psi <- 0.3
  sampledNet <- samplingSBM(A, "dyad", psi)

  ## Perform inference
  missSBM <- missSBM:::missingSBM_fit$new(sampledNet, Q, "dyad")
  out <- missSBM$doVEM(control)

  ## Sanity check
  expect_is(missSBM, "missingSBM_fit")
  expect_is(missSBM$fittedSBM, "SBM_fit_nocovariate")
  expect_is(missSBM$fittedSampling, "dyadSampling_fit")
  expect_is(missSBM$sampledNetwork, "sampledNetwork")
  expect_equal(out, missSBM$monitoring)

  ## Consistency
  tol <- 1e-2

  ## Optimization success
  expect_gt(diff(range(out$objective, na.rm = TRUE)), 0)
  ## SBM: parameters estimation
  expect_lt(sum((missSBM$fittedSBM$connectParam - mySBM$connectParam)^2)/(Q*(Q + 1)/2), tol)
  ## sampling design: parameters estimation
  expect_lt(abs(missSBM$fittedSampling$parameters - psi), tol)
  ## clustering
  tol <- .8
  expect_gt(ARI(missSBM$fittedSBM$memberships, mySBM$memberships), tol)

})

test_that("miss SBM with node sampling works", {

  ## sampled the network
  psi <- 0.2
  sampledNet <- samplingSBM(A, "node", psi)

  ## Perform inference
  missSBM <- missSBM:::missingSBM_fit$new(sampledNet, Q, "node")
  out <- missSBM$doVEM(control)

  ## Sanity check
  expect_is(missSBM, "missingSBM_fit")
  expect_is(missSBM$fittedSBM, "SBM_fit_nocovariate")
  expect_is(missSBM$fittedSampling, "nodeSampling_fit")
  expect_is(missSBM$sampledNetwork, "sampledNetwork")

  ## Consistency
  tol <- 1e-2
  ## Optimization success
  expect_gt(diff(range(out$objective, na.rm = TRUE)), 0)
  ## SBM: parameters estimation
  expect_lt(sum((missSBM$fittedSBM$connectParam - mySBM$connectParam)^2)/(Q*(Q + 1)/2), tol)
  ## sampling design: parameters estimation
  expect_lt(abs(missSBM$fittedSampling$parameters - psi)^2, tol)
  ## clustering
  tol <- .8
  expect_gt(ARI(missSBM$fittedSBM$memberships, mySBM$memberships), tol)

})

test_that("miss SBM with double-standard sampling works", {

  psi <- c(.3, .6)
  sampledNet <- samplingSBM(A, "double-standard", psi)

  ## Perform inference
  missSBM <- missSBM:::missingSBM_fit$new(sampledNet, Q, "double-standard")
  out <- missSBM$doVEM(control)

  ## Sanity check
  expect_is(missSBM, "missingSBM_fit")
  expect_is(missSBM$fittedSBM, "SBM_fit_nocovariate")
  expect_is(missSBM$fittedSampling, "doubleStandardSampling_fit")
  expect_is(missSBM$sampledNetwork, "sampledNetwork")

  ## Consistency
  tol <- 1e-2
  ## Optimization success
  expect_gt(diff(range(out$objective, na.rm = TRUE)), 0)
  ## SBM: parameters estimation
  expect_lt(sum((missSBM$fittedSBM$connectParam - mySBM$connectParam)^2)/(Q*(Q + 1)/2), tol)
  ## sampling design: parameters estimation
  expect_lt(sum((missSBM$fittedSampling$parameters - psi)^2/2), tol)
  ## clustering
  tol <- .9
  expect_gt(ARI(missSBM$fittedSBM$memberships, mySBM$memberships), tol)

})

test_that("miss SBM with block sampling works", {

  psi <- c(.1, .3, .2, .5, .7)
  sampledNet <- samplingSBM(A, "block-node", psi, mySBM$memberships)
  ## Perform inference
  missSBM <- missSBM:::missingSBM_fit$new(sampledNet, Q, "block-node")
  out <- missSBM$doVEM(control)

  ## Sanity check
  expect_is(missSBM, "missingSBM_fit")
  expect_is(missSBM$fittedSBM, "SBM_fit_nocovariate")
  expect_is(missSBM$fittedSampling, "blockSampling_fit")
  expect_is(missSBM$sampledNetwork, "sampledNetwork")

  ## Consistency
  tol <- 1e-2
  ## Optimization success
  expect_gt(diff(range(out$objective, na.rm = TRUE)), 0)
  ## SBM: parameters estimation
  expect_lt(sum((missSBM$fittedSBM$connectParam - mySBM$connectParam)^2)/(Q*(Q + 1)/2), tol)
  ## sampling design: parameters estimation
  expect_lt(sum((sort(missSBM$fittedSampling$parameters) - sort(psi))^2/Q), tol)
  ## clustering
  tol <- .9
  expect_gt(ARI(missSBM$fittedSBM$memberships, mySBM$memberships), tol)

})

test_that("miss SBM with degree sampling works", {

  psi <- c(-5, .1)
  sampledNet <- samplingSBM(A, "degree", psi)
  ## Perform inference
  missSBM <- missSBM:::missingSBM_fit$new(sampledNet, Q, "degree")
  out <- missSBM$doVEM(control)

  ## Sanity check
  expect_is(missSBM, "missingSBM_fit")
  expect_is(missSBM$fittedSBM, "SBM_fit_nocovariate")
  expect_is(missSBM$fittedSampling, "degreeSampling_fit")
  expect_is(missSBM$sampledNetwork, "sampledNetwork")

  ## Consistency
  tol <- 1e-2
  ## Optimization success
  expect_gte(diff(range(out$objective, na.rm = TRUE)), 0)
  ## SBM: parameters estimation
  ## FIXME: expect_lt(sum((missSBM$fittedSBM$connectParam - mySBM$connectParam)^2)/(Q*(Q + 1)/2), tol)
  ## sampling design: parameters estimation
  ## FIXME: this does work!!! expect_lt(sum((sort(missSBM$fittedSampling$parameters) - sort(psi))^2/Q), tol)
  ## clustering
  tol <- .9
  expect_gt(ARI(missSBM$fittedSBM$memberships, mySBM$memberships), tol)

})

