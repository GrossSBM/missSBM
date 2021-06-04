context("test network sampling fit (Class networkSampling_fit and chidren)")

N_cov <- 300
M <- 4
source("utils_test.R", local = TRUE)

test_that("Parameter estimation in dyad-centered sampling with covariates", {

  sampler_undirected_cov$rNetwork(store = TRUE)

  ### Draw a SBM model (Bernoulli, undirected) with covariates
  intercept  <- .5

  adjMatrix  <- missSBM::observeNetwork(sampler_undirected_cov$networkData, "covar-dyad", covarParam, covariates = covarList_undirected, intercept = intercept)
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix, covarList_undirected)

  fittedSampling <- missSBM:::covarDyadSampling_fit$new(net)
  expect_is(fittedSampling, "covarDyadSampling_fit")

  tolerance <- .1
  expect_equal(fittedSampling$df, 1 + length(covarParam))
  expect_equal(fittedSampling$penalty, log(net$nbDyads) * (1 + length(covarParam)) )
  expect_lt(fittedSampling$vExpec, 0)
})

test_that("Parameter estimation in dyad-centered sampling with covariates but ignoring them", {

  psi <- 0.8
  adjMatrix  <- missSBM::observeNetwork(sampler_undirected_cov$networkData, "dyad", psi, covariates = covarList_undirected)
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix, covarList_undirected)

  fittedSampling <- missSBM:::dyadSampling_fit$new(net)
  expect_is(fittedSampling, "dyadSampling_fit")

  tolerance <- .1
  expect_lt(error(fittedSampling$parameters, psi), tolerance)
  expect_equal(fittedSampling$df, 1)
  expect_equal(fittedSampling$penalty, log(net$nbDyads))
  expect_lt(fittedSampling$vExpec, 0)
})

test_that("Parameter estimation in node-centered sampling with covariates", {

  sampler_undirected_cov_node$rNetwork(store = TRUE)

  ### Draw a SBM model (Bernoulli, undirected) with covariates
  intercept  <- .5

  adjMatrix <- missSBM::observeNetwork(sampler_undirected_cov_node$networkData, "covar-node", covarParam, covariates = covarList_node, intercept = intercept, similarity = missSBM:::l1_similarity)
  partlyObservedNet <- missSBM:::partlyObservedNetwork$new(adjMatrix, covarList_node, missSBM:::l1_similarity)

  fittedSampling <- missSBM:::covarNodeSampling_fit$new(partlyObservedNet)
  expect_is(fittedSampling, "covarNodeSampling_fit")

  tolerance <- .3
  expect_equal(fittedSampling$df, 1 + length(covarParam))
  expect_equal(fittedSampling$penalty, log(N_cov) * (1 + length(covarParam)))
  expect_lt(fittedSampling$vExpec, 0)
})

test_that("Parameter estimation in node-centered sampling with covariates but ignoring them", {

  psi <- 0.9
  adjMatrix <- missSBM::observeNetwork(sampler_undirected_cov_node$networkData, "node", psi)
  partlyObservedNet <- missSBM:::partlyObservedNetwork$new(adjMatrix, covarList_node, missSBM:::l1_similarity)

  fittedSampling <- missSBM:::nodeSampling_fit$new(partlyObservedNet)
  expect_is(fittedSampling, "nodeSampling_fit")

  tolerance <- .1
  expect_lt(error(fittedSampling$parameters, psi), tolerance)
  expect_equal(fittedSampling$df, 1)
  expect_equal(fittedSampling$penalty, log(N_cov))
  expect_lt(fittedSampling$vExpec, 0)
})
