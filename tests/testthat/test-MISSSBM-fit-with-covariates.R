context("test missSBM-fit with covariates")

library(aricode)
source("utils_test.R")

## ========================================================================
## A SBM model with covariates

set.seed(787)
N <- 200
Q <- 3
alpha <- rep(1,Q)/Q                     # mixture parameter
pi <- diag(.45, Q, Q) + .05                 # connectivity matrix
gamma <- missSBM:::logit(pi)
directed <- FALSE

### Draw a SBM model (Bernoulli, undirected) with covariates
M <- 2
covariates_node <- replicate(M, rnorm(N,mean = 0, sd = 1), simplify = FALSE)
covarMatrix <- simplify2array(covariates_node)
covarArray  <- missSBM:::getCovarArray(covarMatrix, missSBM:::l1_similarity)
covariates_dyad <- lapply(seq(dim(covarArray)[3]), function(x) covarArray[ , , x])
covarParam  <- rnorm(M, -1, 1)

## control parameter for the VEM
control <- list(threshold = 1e-4, maxIter = 200, fixPointIter = 5, trace = TRUE)

## Consistency
tol_truth <- 1e-2
tol_ARI   <- .9

test_that("missSBM with covariates and dyad sampling works", {

  sbm <- missSBM::simulate(N, alpha, gamma, directed, covariates_dyad, covarParam)

  ## ACCOUNT FOR COVARIATES IN THE SAMPLING

  ## sampled the network
  sampledNet <- missSBM::sample(sbm$adjacencyMatrix, "covar-dyad", covarParam, covariates = covariates_dyad)

  ## Perform inference
  missSBM <- missSBM:::missSBM_fit$new(sampledNet, Q, "covar-dyad", "spectral", TRUE)
  out <- missSBM$doVEM(control)

  ## Sanity check
  expect_is(missSBM, "missSBM_fit")
  expect_is(missSBM$fittedSBM, "SBM_fit_covariates")
  expect_is(missSBM$fittedSampling, "covarDyadSampling_fit")
  expect_is(missSBM$sampledNetwork, "sampledNetwork")
  expect_equal(out, missSBM$monitoring)

  ## Optimization success
  expect_gte(diff(range(out$objective, na.rm = TRUE)), 0)

  ## SBM: parameters estimation
  expect_lt(error(missSBM$fittedSBM$mixtureParam, sbm$mixtureParam, sort = TRUE), tol_truth)

  expect_lt(error(logistic(missSBM$fittedSBM$connectParam), pi), tol_truth*10)

  ## sampling design: parameters estimation
  expect_lt(error(missSBM$fittedSBM$covarParam, sbm$covarParam), tol_truth)

  ## clustering
  expect_gt(aricode::ARI(missSBM$fittedSBM$memberships, sbm$memberships), tol_ARI)

  ## DO NOT ACCOUNT FOR COVARIATES IN THE SAMPLING (JUST IN THE SBM)

  ## sampled the network
  sampledNet <- missSBM::sample(sbm$adjacencyMatrix, "dyad", 0.9, covariates = covariates_dyad)

  ## Perform inference
  missSBM <- missSBM:::missSBM_fit$new(sampledNet, Q, "dyad", "spectral", TRUE)
  out <- missSBM$doVEM(control)

  ## Sanity check
  expect_is(missSBM, "missSBM_fit")
  expect_is(missSBM$fittedSBM, "SBM_fit_covariates")
  expect_is(missSBM$fittedSampling, "dyadSampling_fit")
  expect_is(missSBM$sampledNetwork, "sampledNetwork")
  expect_equal(out, missSBM$monitoring)

  ## Optimization success
  expect_gte(diff(range(out$objective, na.rm = TRUE)), 0)

  ## SBM: parameters estimation
  expect_lt(error(missSBM$fittedSBM$mixtureParam, sbm$mixtureParam, sort = TRUE), tol_truth)

  expect_lt(error(logistic(missSBM$fittedSBM$connectParam), pi), tol_truth*10)

  ## sampling design: parameters estimation
  expect_lt(error(missSBM$fittedSBM$covarParam, sbm$covarParam), tol_truth)

  ## clustering
  expect_gt(aricode::ARI(missSBM$fittedSBM$memberships, sbm$memberships), tol_ARI)

})

test_that("miss SBM with covariates and node sampling works", {

  sbm <- missSBM::simulate(N, alpha, gamma, directed, covariates_dyad, covarParam)

  ## sampled the network
  sampledNet <- missSBM::sample(sbm$adjacencyMatrix, "covar-node", covarParam, covariates = covariates_node)

  ## Perform inference
  missSBM <- missSBM:::missSBM_fit$new(sampledNet, Q, "covar-node", clusterInit = "spectral", TRUE)
  out <- missSBM$doVEM(control)

  ## Sanity check
  expect_is(missSBM, "missSBM_fit")
  expect_is(missSBM$fittedSBM, "SBM_fit_covariates")
  expect_is(missSBM$fittedSampling, "covarNodeSampling_fit")
  expect_is(missSBM$sampledNetwork, "sampledNetwork")

  ## Optimization success
  expect_gte(diff(range(out$objective, na.rm = TRUE)), 0)

  ## SBM: parameters estimation
  expect_lt(error(missSBM$fittedSBM$mixtureParam, sbm$mixtureParam, sort = TRUE), tol_truth)

  expect_lt(error(logistic(missSBM$fittedSBM$connectParam), pi), tol_truth)

  ## sampling design: parameters estimation
  expect_lt(error(missSBM$fittedSampling$parameters, sbm$covarParam), 10 * tol_truth)

  ## clustering
  expect_gt(aricode::ARI(missSBM$fittedSBM$memberships, sbm$memberships), tol_ARI)

  ## do not account for covariate in the sampling (just in SBM)

  ## sampled the network
  sampledNet <- missSBM::sample(sbm$adjacencyMatrix, "node", 0.9, covariates = covariates_node)

  ## Perform inference
  missSBM <- missSBM:::missSBM_fit$new(sampledNet, Q, "node", clusterInit = "spectral", TRUE)
  out <- missSBM$doVEM(control)

  ## Sanity check
  expect_is(missSBM, "missSBM_fit")
  expect_is(missSBM$fittedSBM, "SBM_fit_covariates")
  expect_is(missSBM$fittedSampling, "nodeSampling_fit")
  expect_is(missSBM$sampledNetwork, "sampledNetwork")

  ## Optimization success
  expect_gte(diff(range(out$objective, na.rm = TRUE)), 0)

  ## SBM: parameters estimation
  expect_lt(error(missSBM$fittedSBM$mixtureParam, sbm$mixtureParam, sort = TRUE), tol_truth)

  expect_lt(error(logistic(missSBM$fittedSBM$connectParam), pi), tol_truth)

  ## sampling design: parameters estimation
  expect_lt(error(missSBM$fittedSampling$parameters, 0.9), 10 * tol_truth)

  ## clustering
  expect_gt(aricode::ARI(missSBM$fittedSBM$memberships, sbm$memberships), tol_ARI)

})
