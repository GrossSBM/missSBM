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

#### There is a problem here !
  tolerance <- .1
  expect_lt(error(fittedSampling$parameters, c(intercept, covarParam)), tolerance)
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

  tolerance <- 1e-2
  expect_lt(error(fittedSampling$parameters, psi), tolerance)
  expect_equal(fittedSampling$df, 1)
  expect_equal(fittedSampling$penalty, log(net$nbDyads))
  expect_lt(fittedSampling$vExpec, 0)
})

# test_that("Parameter estimation in node-centered sampling with covariates", {
#
#   M <- 10
#   covariates_node <- replicate(M, rnorm(N ,mean = 0, sd = 1), simplify = FALSE)
#   covarArray <- missSBM:::getCovarArray(simplify2array(covariates_node), missSBM:::l1_similarity)
#   covariates <- lapply(1:M, function(m) covarArray[,,m])
#   covarParam <- rnorm(M, 0, 1)
#   intercept  <- .5
#
#   sbm <- sbm::sampleSimpleSBM(N, pi, theta, covariates = covariates, covariatesParam = covarParam)
#   adjMatrix <- missSBM::observeNetwork(sbm$networkData, "covar-node", covarParam, covariates = covariates_node, intercept = intercept, similarity = missSBM:::l1_similarity)
#   partlyObservedNet <- missSBM:::partlyObservedNetwork$new(adjMatrix, covariates_node, missSBM:::l1_similarity)
#
#   fittedSampling <- missSBM:::covarNodeSampling_fit$new(partlyObservedNet, simplify2array(covariates_node))
#   expect_is(fittedSampling, "covarNodeSampling_fit")
#   expect_true(all(fittedSampling$prob_obs > 0, fittedSampling$prob_obs < 1))
#
#   tolerance <- .1
#   expect_lt(error(fittedSampling$parameters, c(intercept, covarParam)), tolerance)
#   expect_equal(fittedSampling$df, 1 + length(covarParam))
#   expect_equal(fittedSampling$penalty, log(N) * (1 + length(covarParam)))
#   expect_lt(fittedSampling$vExpec, 0)
# })
#
# test_that("Parameter estimation in node-centered sampling with covariates but ignoring them", {
#
#   M <- 10
#   covariates_node <- replicate(M, rnorm(N ,mean = 0, sd = 1), simplify = FALSE)
#   covarArray <- missSBM:::getCovarArray(simplify2array(covariates_node), missSBM:::l1_similarity)
#   covariates <- lapply(1:M, function(m) covarArray[,,m])
#   covarParam  <- rnorm(M, 0, 1)
#
#   sbm <- sbm::sampleSimpleSBM(N, pi, theta, covariates = covariates, covariatesParam = covarParam)
#
#   adjMatrix <- missSBM::observeNetwork(sbm$networkData, "node", .9, covariates = covariates_node)
#   partlyObservedNet <- missSBM:::partlyObservedNetwork$new(adjMatrix, covariates_node, missSBM:::l1_similarity)
#
#   fittedSampling <- missSBM:::nodeSampling_fit$new(partlyObservedNet, simplify2array(covariates_node))
#   expect_is(fittedSampling, "nodeSampling_fit")
#   expect_true(all(fittedSampling$prob_obs > 0, fittedSampling$prob_obs < 1))
#
#   tolerance <- .01
#   expect_lt(sum((fittedSampling$parameters - .9)^2), tolerance)
#   expect_equal(fittedSampling$df, 1)
#   expect_equal(fittedSampling$penalty, log(N))
#   expect_lt(fittedSampling$vExpec, 0)
# })
