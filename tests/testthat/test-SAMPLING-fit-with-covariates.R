context("test network sampling fit (Class networkSampling_fit and chidren)")

set.seed(1234)
### A SBM model : ###
N <- 200
Q <- 3
pi <- rep(1, Q)/Q                        # block proportion parameter
theta <- list(mean = diag(.45,Q) + .05 ) # connectivity matrix

test_that("Parameter estimation in dyad-centered sampling with covariates", {

  ### Draw a SBM model (Bernoulli, undirected) with covariates
  M <- 5
  covariates <- replicate(M, matrix(rnorm(N * N ,mean = 0, sd = 1), N, N), simplify = FALSE)
  covarParam  <- rnorm(M, 0, 1)
  intercept  <- .5

  sbm <- sbm::sampleSimpleSBM(N, pi, theta, covariates = covariates, covariatesParam = covarParam)

  adjMatrix  <- missSBM::observeNetwork(sbm$networkData, "covar-dyad", covarParam, covariates = covariates, intercept = intercept)
  partlyObservedNet <- missSBM:::partlyObservedNetwork$new(adjMatrix, covariates, missSBM:::l1_similarity)

  fittedSampling <- missSBM:::covarDyadSampling_fit$new(partlyObservedNet, sbm$covarArray)
  expect_is(fittedSampling, "covarDyadSampling_fit")
  expect_true(all(fittedSampling$prob_obs > 0, fittedSampling$prob_obs < 1))

  tolerance <- .5
  expect_lt(sum((fittedSampling$parameters - c(intercept, covarParam))^2)/(M + 1), tolerance)
  expect_equal(fittedSampling$df, 1 + length(covarParam))
  expect_equal(fittedSampling$penalty, log(N * (N - 1)/2) * (1 + length(covarParam)) )
  expect_lt(fittedSampling$vExpec, 0)
})

test_that("Parameter estimation in dyad-centered sampling with covariates but ignoring them", {

  ### Draw a SBM model (Bernoulli, undirected) with covariates
  M <- 5
  covariates <- replicate(M, matrix(rnorm(N * N ,mean = 0, sd = 1), N, N), simplify = FALSE)
  covarParam  <- rnorm(M, 0, 1)
  sbm <- sbm::sampleSimpleSBM(N, pi, theta, covariates = covariates, covariatesParam = covarParam)

  adjMatrix <- missSBM::observeNetwork(sbm$networkData, "dyad", parameters = .9, covariates = covariates)
  partlyObservedNet <- missSBM:::partlyObservedNetwork$new(adjMatrix, covariates, missSBM:::l1_similarity)

  fittedSampling <- missSBM:::dyadSampling_fit$new(partlyObservedNet, sbm$covarArray)
  expect_is(fittedSampling, "dyadSampling_fit")
  expect_true(all(fittedSampling$prob_obs > 0, fittedSampling$prob_obs < 1))

  tolerance <- 1e-2
  expect_lt(sum((fittedSampling$parameters - .9)^2/(M + 1)), tolerance)
  expect_equal(fittedSampling$df, 1)
  expect_equal(fittedSampling$penalty, log(N * (N - 1)/2))
  expect_lt(fittedSampling$vExpec, 0)
})

test_that("Parameter estimation in node-centered sampling with covariates", {

  M <- 10
  covariates_node <- replicate(M, rnorm(N ,mean = 0, sd = 1), simplify = FALSE)
  covarArray <- missSBM:::getCovarArray(simplify2array(covariates_node), missSBM:::l1_similarity)
  covariates <- lapply(1:M, function(m) covarArray[,,m])
  covarParam <- rnorm(M, 0, 1)
  intercept  <- .5

  sbm <- sbm::sampleSimpleSBM(N, pi, theta, covariates = covariates, covariatesParam = covarParam)
  adjMatrix <- missSBM::observeNetwork(sbm$networkData, "covar-node", covarParam, covariates = covariates_node, intercept = intercept, similarity = missSBM:::l1_similarity)
  partlyObservedNet <- missSBM:::partlyObservedNetwork$new(adjMatrix, covariates_node, missSBM:::l1_similarity)

  fittedSampling <- missSBM:::covarNodeSampling_fit$new(partlyObservedNet, simplify2array(covariates_node))
  expect_is(fittedSampling, "covarNodeSampling_fit")
  expect_true(all(fittedSampling$prob_obs > 0, fittedSampling$prob_obs < 1))

  tolerance <- .25
  expect_lt(error(fittedSampling$parameters,c(intercept, covarParam)), tolerance)
  expect_equal(fittedSampling$df, 1 + length(covarParam))
  expect_equal(fittedSampling$penalty, log(N) * (1 + length(covarParam)))
  expect_lt(fittedSampling$vExpec, 0)
})

test_that("Parameter estimation in node-centered sampling with covariates but ignoring them", {

  M <- 10
  covariates_node <- replicate(M, rnorm(N ,mean = 0, sd = 1), simplify = FALSE)
  covarArray <- missSBM:::getCovarArray(simplify2array(covariates_node), missSBM:::l1_similarity)
  covariates <- lapply(1:M, function(m) covarArray[,,m])
  covarParam  <- rnorm(M, 0, 1)

  sbm <- sbm::sampleSimpleSBM(N, pi, theta, covariates = covariates, covariatesParam = covarParam)

  adjMatrix <- missSBM::observeNetwork(sbm$networkData, "node", .9, covariates = covariates_node)
  partlyObservedNet <- missSBM:::partlyObservedNetwork$new(adjMatrix, covariates_node, missSBM:::l1_similarity)

  fittedSampling <- missSBM:::nodeSampling_fit$new(partlyObservedNet, simplify2array(covariates_node))
  expect_is(fittedSampling, "nodeSampling_fit")
  expect_true(all(fittedSampling$prob_obs > 0, fittedSampling$prob_obs < 1))

  tolerance <- .2
  expect_lt(sum((fittedSampling$parameters - .9)^2), tolerance)
  expect_equal(fittedSampling$df, 1)
  expect_equal(fittedSampling$penalty, log(N))
  expect_lt(fittedSampling$vExpec, 0)
})

test_that("Parameter estimation in degree sampling", {

  ### Draw a SBM model (Bernoulli, undirected) with covariates
  sbm <- sbm::sampleSimpleSBM(N, pi, theta)

  psi <- c(-.5,0.01)
  adjMatrix <- missSBM::observeNetwork(sbm$networkData, "degree", psi)
  partlyObservedNet <- missSBM:::partlyObservedNetwork$new(adjMatrix)

  Z0 <- missSBM:::clustering_indicator(sbm$memberships)
  fittedSampling <- suppressWarnings(missSBM:::degreeSampling_fit$new(partlyObservedNet, Z0, sbm$connectParam$mean))

  # tolerance <- 1 ## not expected good after one iterate
  # expect_lt(sum((fittedSampling$parameters - psi)^2), tolerance)
  expect_equal(fittedSampling$df, length(psi))
  expect_equal(fittedSampling$penalty, log(N) * length(psi))
  expect_lt(fittedSampling$vExpec, 0)
})

