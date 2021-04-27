context("test missSBM-fit with covariates")

N_cov <- 150
Q <- 2
M <- 1
source("utils_test.R", local = TRUE)

## control parameter for the VEM
control <- list(threshold = 1e-4, maxIter = 50, fixPointIter = 5, trace = TRUE)

## Consistency
tol_truth <- .2
tol_ARI   <- .7

test_that("missSBM with covariates and dyad sampling works", {

  sampler_undirected_cov$rNetwork(store = TRUE)

  ## ACCOUNT FOR COVARIATES IN THE SAMPLING

  ## sampled the network
  adjMatrix <- missSBM::observeNetwork(sampler_undirected_cov$networkData, "covar-dyad", covarParam, covariates = covarList_undirected)

  ## Prepare network data for estimation with missing data
  partlyObservedNet <- missSBM:::partlyObservedNetwork$new(adjMatrix, covarList_undirected)
  cl <- partlyObservedNet$clustering(Q)[[1]]

  ## Perform inference
  missSBM <- missSBM:::missSBM_fit$new(partlyObservedNet, "covar-dyad", cl, TRUE)
  out <- missSBM$doVEM(control)

  ## Sanity check
  expect_is(missSBM, "missSBM_fit")
  expect_is(missSBM$fittedSBM, "SimpleSBM_fit_withCov")
  expect_is(missSBM$fittedSampling, "covarDyadSampling_fit")
  expect_equal(out, missSBM$monitoring)

  ## Optimization success
  expect_gte(diff(range(out$elbo, na.rm = TRUE)), 0)

  ## SBM: parameters estimation
  expect_lt(error(missSBM$fittedSBM$blockProp, sampler_undirected_cov$blockProp, sort = TRUE), tol_truth)

  expect_lt(error(missSBM$fittedSBM$connectParam$mean, sampler_undirected_cov$connectParam$mean), tol_truth)

  ## sampling design: parameters estimation
  expect_lt(error(missSBM$fittedSBM$covarParam, sampler_undirected_cov$covarParam), 0.25)

  ## clustering
  expect_gt(aricode::ARI(missSBM$fittedSBM$memberships, sampler_undirected_cov$memberships), tol_ARI)

  ## DO NOT ACCOUNT FOR COVARIATES IN THE SAMPLING (JUST IN THE SBM)

  ## sampled the network
  adjMatrix <- missSBM::observeNetwork(sampler_undirected_cov$networkData, "dyad", 0.9)

  ## Prepare network data for estimation with missing data
  partlyObservedNet <- missSBM:::partlyObservedNetwork$new(adjMatrix, covarList_undirected)
  cl <- partlyObservedNet$clustering(Q)[[1]]

  ## Perform inference
  missSBM <- missSBM:::missSBM_fit$new(partlyObservedNet, "dyad", cl, TRUE)
  out <- missSBM$doVEM(control)

  ## Sanity check
  expect_is(missSBM, "missSBM_fit")
  expect_is(missSBM$fittedSBM, "SimpleSBM_fit_withCov")
  expect_is(missSBM$fittedSampling, "dyadSampling_fit")
  expect_equal(out, missSBM$monitoring)

  ## Optimization success
  expect_gte(diff(range(out$elbo, na.rm = TRUE)), 0)

  ## SBM: parameters estimation
  expect_lt(error(missSBM$fittedSBM$blockProp, sampler_undirected_cov$blockProp, sort = TRUE), tol_truth)

  expect_lt(error(missSBM$fittedSBM$connectParam$mean, sampler_undirected_cov$connectParam$mean), tol_truth)

  ## sampling design: parameters estimation
  expect_lt(error(missSBM$fittedSBM$covarParam, sampler_undirected_cov$covarParam), 0.25)

  ## clustering
  expect_gt(aricode::ARI(missSBM$fittedSBM$memberships, sampler_undirected_cov$memberships), tol_ARI)

})

test_that("miss SBM with covariates and node sampling works", {

  sampler_undirected_cov_node$rNetwork(store = TRUE)

  ## sampled the network
  intercept  <- .5
  adjMatrix <- missSBM::observeNetwork(sampler_undirected_cov_node$networkData, "covar-node", covarParam, covariates = covarList_node, intercept = intercept, similarity = missSBM:::l1_similarity)

  ## Prepare network data for estimation with missing data
  partlyObservedNet <- missSBM:::partlyObservedNetwork$new(adjMatrix, covarList_node, missSBM:::l1_similarity)
  cl <- partlyObservedNet$clustering(Q)[[1]]

  ## Perform inference
  missSBM <- missSBM:::missSBM_fit$new(partlyObservedNet, "covar-node", cl, TRUE)
  out <- missSBM$doVEM(control)

  ## Sanity check
  expect_is(missSBM, "missSBM_fit")
  expect_is(missSBM$fittedSBM, "SimpleSBM_fit_withCov")
  expect_is(missSBM$fittedSampling, "covarNodeSampling_fit")

  ## Optimization success
  expect_gte(diff(range(out$elbo, na.rm = TRUE)), 0)

  ## SBM: parameters estimation
  expect_lt(error(missSBM$fittedSBM$blockProp, sampler_undirected_cov_node$blockProp, sort = TRUE), tol_truth)

  expect_lt(error(missSBM$fittedSBM$connectParam$mean, sampler_undirected_cov_node$connectParam$mean), tol_truth)

  ## sampling design: parameters estimation
  expect_lt(error(missSBM$fittedSBM$covarParam, sampler_undirected_cov_node$covarParam), tol_truth * 3)

  ## clustering
  expect_gt(aricode::ARI(missSBM$fittedSBM$memberships, sampler_undirected_cov_node$memberships), tol_ARI)

  ## sampled the network
  psi <- 0.9
  adjMatrix <- missSBM::observeNetwork(sampler_undirected_cov_node$networkData, "node", psi)

  ## Prepare network data for estimation with missing data
  partlyObservedNet <- missSBM:::partlyObservedNetwork$new(adjMatrix, covarList_node, missSBM:::l1_similarity)
  cl <- partlyObservedNet$clustering(Q)[[1]]

  ## Perform inference
  missSBM <- missSBM:::missSBM_fit$new(partlyObservedNet, "node", cl, TRUE)
  out <- missSBM$doVEM(control)

  ## Sanity check
  expect_is(missSBM, "missSBM_fit")
  expect_is(missSBM$fittedSBM, "SimpleSBM_fit_withCov")
  expect_is(missSBM$fittedSampling, "nodeSampling_fit")

  ## Optimization success
  expect_gte(diff(range(out$elbo, na.rm = TRUE)), 0)

  ## SBM: parameters estimation
  expect_lt(error(missSBM$fittedSBM$blockProp, sampler_undirected_cov_node$blockProp, sort = TRUE), tol_truth)
  expect_lt(error(missSBM$fittedSBM$connectParam$mean, sampler_undirected_cov_node$connectParam$mean), tol_truth)

  ## sampling design: parameters estimation
  expect_lt(error(missSBM$fittedSampling$parameters, psi), tol_truth)

  ## clustering
  expect_gt(aricode::ARI(missSBM$fittedSBM$memberships, sampler_undirected_cov_node$memberships), tol_ARI)

})
