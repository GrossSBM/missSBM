context("test missSBM-fit with covariates")

library(aricode)
source("utils_test.R")

## ========================================================================
## A SBM model with covariates

set.seed(1827)
N <- 100
Q <- 2
pi <- rep(1,Q)/Q                            # mixture parameter
theta <- list(mean = diag(.45, Q, Q) + .05) # connectivity matrix

### Draw a SBM model (Bernoulli, undirected) with covariates
M <- 1
covariates_node <- replicate(M, rnorm(N,mean = 0, sd = 1), simplify = FALSE)
covarMatrix <- simplify2array(covariates_node)
covarArray  <- missSBM:::getCovarArray(covarMatrix, missSBM:::l1_similarity)
covariates_dyad <- lapply(seq(dim(covarArray)[3]), function(x) covarArray[ , , x])
covarParam  <- rnorm(M, -1, 1)

## control parameter for the VEM
control <- list(threshold = 1e-4, maxIter = 50, fixPointIter = 2, trace = TRUE)

## Consistency
tol_truth <- 1e-2
tol_ARI   <- .9

test_that("missSBM with covariates and dyad sampling works", {

  sbm <- sbm::sampleSimpleSBM(N, pi, theta, covariates = covariates_dyad, covariatesParam = covarParam)

  ## ACCOUNT FOR COVARIATES IN THE SAMPLING

  ## sampled the network
  adjMatrix <- missSBM::observeNetwork(sbm$networkData, "covar-dyad", covarParam, covariates = covariates_dyad)

  ## Prepare network data for estimation with missing data
  partlyObservedNet <- missSBM:::partlyObservedNetwork$new(adjMatrix, covariates_dyad, missSBM:::l1_similarity)
  cl <- partlyObservedNet$clustering(Q)[[Q]]

  ## Perform inference
  missSBM <- missSBM:::missSBM_fit$new(partlyObservedNet, "covar-dyad", cl, TRUE)
  out <- missSBM$doVEM(control)

  ## Sanity check
  expect_is(missSBM, "missSBM_fit")
  expect_is(missSBM$fittedSBM, "SimpleSBM_fit_missSBM")
  expect_is(missSBM$fittedSampling, "covarDyadSampling_fit")
  expect_equal(out, missSBM$monitoring)

  ## Optimization success
  expect_gte(diff(range(out$objective, na.rm = TRUE)), 0)

  ## SBM: parameters estimation
  expect_lt(error(missSBM$fittedSBM$blockProp, sbm$blockProp, sort = TRUE), tol_truth)

  expect_lt(error(missSBM$fittedSBM$connectParam$mean, theta$mean), tol_truth*10)

  ## sampling design: parameters estimation
  expect_lt(error(missSBM$fittedSBM$covarParam, sbm$covarParam), 0.25)

  ## clustering
  expect_gt(aricode::ARI(missSBM$fittedSBM$memberships, sbm$memberships), tol_ARI)

  ## DO NOT ACCOUNT FOR COVARIATES IN THE SAMPLING (JUST IN THE SBM)

  ## sampled the network
  adjMatrix <- missSBM::observeNetwork(sbm$networkData, "dyad", 0.9)

  ## Prepare network data for estimation with missing data
  partlyObservedNet <- missSBM:::partlyObservedNetwork$new(adjMatrix, covariates_dyad, missSBM:::l1_similarity)
  cl <- partlyObservedNet$clustering(Q)[[1]]

  ## Perform inference
  missSBM <- missSBM:::missSBM_fit$new(partlyObservedNet, "dyad", cl, TRUE)
  out <- missSBM$doVEM(control)

  ## Sanity check
  expect_is(missSBM, "missSBM_fit")
  expect_is(missSBM$fittedSBM, "SimpleSBM_fit_missSBM")
  expect_is(missSBM$fittedSampling, "dyadSampling_fit")
  expect_equal(out, missSBM$monitoring)

  ## Optimization success
  expect_gte(diff(range(out$objective, na.rm = TRUE)), 0)

  ## SBM: parameters estimation
  expect_lt(error(missSBM$fittedSBM$blockProp, sbm$blockProp, sort = TRUE), tol_truth)

  expect_lt(error(missSBM$fittedSBM$connectParam$mean, theta$mean), tol_truth)

  ## sampling design: parameters estimation
  expect_lt(error(missSBM$fittedSBM$covarParam, sbm$covarParam), 0.25)

  ## clustering
  expect_gt(aricode::ARI(missSBM$fittedSBM$memberships, sbm$memberships), tol_ARI)

})

test_that("miss SBM with covariates and node sampling works", {

  sbm <- sbm::sampleSimpleSBM(N, pi, theta, covariates = covariates_dyad, covariatesParam = covarParam)

  ## sampled the network
  adjMatrix <- missSBM::observeNetwork(sbm$networkData, "covar-node", covarParam, covariates = covariates_node, similarity = missSBM:::l1_similarity)

  ## Prepare network data for estimation with missing data
  partlyObservedNet <- missSBM:::partlyObservedNetwork$new(adjMatrix, covariates_node, missSBM:::l1_similarity)
  cl <- partlyObservedNet$clustering(Q)[[1]]

  ## Perform inference
  missSBM <- missSBM:::missSBM_fit$new(partlyObservedNet, "covar-node", cl, TRUE)
  out <- missSBM$doVEM(control)

  ## Sanity check
  expect_is(missSBM, "missSBM_fit")
  expect_is(missSBM$fittedSBM, "SimpleSBM_fit_missSBM")
  expect_is(missSBM$fittedSampling, "covarNodeSampling_fit")

  ## Optimization success
  expect_gte(diff(range(out$objective, na.rm = TRUE)), 0)

  ## SBM: parameters estimation
  expect_lt(error(missSBM$fittedSBM$blockProp, sbm$blockProp, sort = TRUE), tol_truth)

  expect_lt(error(missSBM$fittedSBM$connectParam$mean, theta$mean), tol_truth)

  ## sampling design: parameters estimation
  expect_lt(error(missSBM$fittedSBM$covarParam, sbm$covarParam), .5)

  ## clustering
  expect_gt(aricode::ARI(missSBM$fittedSBM$memberships, sbm$memberships), tol_ARI)

  ## do not account for covariate in the SBM (just in the sampling)
  sbm <- sbm::sampleSimpleSBM(N, pi, theta)

  ## sampled the network
  adjMatrix <- missSBM::observeNetwork(sbm$networkData, "node", 0.9, covariates = covariates_node, similarity = missSBM:::l1_similarity)

  ## Prepare network data for estimation with missing data
  partlyObservedNet <- missSBM:::partlyObservedNetwork$new(adjMatrix, covariates_node, missSBM:::l1_similarity)
  cl <- partlyObservedNet$clustering(Q)[[1]]

  ## Perform inference
  missSBM <- missSBM:::missSBM_fit$new(partlyObservedNet, "node", cl, FALSE)
  out <- missSBM$doVEM(control)

  ## Sanity check
  expect_is(missSBM, "missSBM_fit")
  expect_is(missSBM$fittedSBM, "SimpleSBM_fit_missSBM")
  expect_is(missSBM$fittedSampling, "nodeSampling_fit")

  ## Optimization success
  expect_gte(diff(range(out$objective, na.rm = TRUE)), 0)

  ## SBM: parameters estimation
  expect_lt(error(missSBM$fittedSBM$blockProp, sbm$blockProp, sort = TRUE), tol_truth)
  expect_lt(error(missSBM$fittedSBM$connectParam$mean, theta$mean), tol_truth)

  ## sampling design: parameters estimation
  expect_lt(error(missSBM$fittedSampling$parameters, 0.9), 10 * tol_truth)

  ## clustering
  expect_gt(aricode::ARI(missSBM$fittedSBM$memberships, sbm$memberships), tol_ARI)

})
