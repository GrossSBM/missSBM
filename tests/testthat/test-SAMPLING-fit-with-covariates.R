context("test network sampling fit (Class networkSampling_fit and chidren)")

set.seed(178303)
### A SBM model : ###
N <- 200
Q <- 3
alpha <- rep(1, Q)/Q                       # mixture parameter
pi <- diag(.45, Q, Q) + .05                   # connectivity matrix
directed <- FALSE                         # if the network is directed or not

test_that("Parameter estimation in dyad-centered sampling with covariates", {

  ### Draw a SBM model (Bernoulli, undirected) with covariates
  M <- 5
  covariates <- replicate(M, matrix(rnorm(N * N ,mean = 0, sd = 1), N, N), simplify = FALSE)
  covarParam  <- rnorm(M, 0, 1)
  intercept  <- .5

  sbm <- missSBM::simulate(N, alpha, log(pi/(1 - pi)), directed, covariates, covarParam)

  adjMatrix  <- missSBM::sample(sbm$adjacencyMatrix, "covar-dyad", covarParam, covariates = covariates, intercept = intercept)
  covar <- missSBM:::format_covariates(covariates, missSBM:::l1_similarity)
  sampledNet <- missSBM:::sampledNetwork$new(adjMatrix, covar$Matrix, covar$Array)

  fittedSampling <- missSBM:::covarDyadSampling_fit$new(sampledNet, sbm$covarArray)
  expect_is(fittedSampling, "covarDyadSampling_fit")
  expect_true(all(fittedSampling$prob_obs > 0, fittedSampling$prob_obs < 1))

  tolerance <- .5
  # expect_lt(sum((fittedSampling$parameters - c(intercept, covarParam))^2)/(M + 1), tolerance)
  expect_equal(fittedSampling$df, 1 + length(covarParam))
  expect_equal(fittedSampling$penalty, log(N * (N - 1)/2) * (1 + length(covarParam)) )
  expect_lt(fittedSampling$vExpec, 0)
})

test_that("Parameter estimation in dyad-centered sampling with covariates but ignoring them", {

  ### Draw a SBM model (Bernoulli, undirected) with covariates
  M <- 5
  covariates <- replicate(M, matrix(rnorm(N * N ,mean = 0, sd = 1), N, N), simplify = FALSE)
  covarParam  <- rnorm(M, 0, 1)
  sbm <- missSBM::simulate(N, alpha, log(pi/(1 - pi)), directed, covariates, covarParam)

  adjMatrix <- missSBM::sample(sbm$adjacencyMatrix, "dyad", parameters = .9, covariates = covariates)
  covar <- missSBM:::format_covariates(covariates, missSBM:::l1_similarity)
  sampledNet <- missSBM:::sampledNetwork$new(adjMatrix, covar$Matrix, covar$Array)

  fittedSampling <- missSBM:::dyadSampling_fit$new(sampledNet, sbm$covarArray)
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

  sbm <- missSBM::simulate(N, alpha, log(pi/(1 - pi)), directed, covariates, covarParam)
  adjMatrix <- missSBM::sample(sbm$adjacencyMatrix, "covar-node", covarParam, covariates = covariates_node, intercept = intercept)
  covar <- missSBM:::format_covariates(covariates_node, missSBM:::l1_similarity)
  sampledNet <- missSBM:::sampledNetwork$new(adjMatrix, covar$Matrix, covar$Array)

  fittedSampling <- missSBM:::covarNodeSampling_fit$new(sampledNet, simplify2array(covariates_node))
  expect_is(fittedSampling, "covarNodeSampling_fit")
  expect_true(all(fittedSampling$prob_obs > 0, fittedSampling$prob_obs < 1))

  tolerance <- .05
  expect_lt(sum((fittedSampling$parameters - c(intercept, covarParam))^2)/(M + 1), tolerance)
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

  sbm <- missSBM::simulate(N, alpha, log(pi/(1 - pi)), directed, covariates, covarParam)

  adjMatrix <- missSBM::sample(sbm$adjacencyMatrix, "node", .9, covariates = covariates_node)
  covar <- missSBM:::format_covariates(covariates_node, missSBM:::l1_similarity)
  sampledNet <- missSBM:::sampledNetwork$new(adjMatrix, covar$Matrix, covar$Array)

  fittedSampling <- missSBM:::nodeSampling_fit$new(sampledNet, simplify2array(covariates_node))
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
  sbm <- missSBM::simulate(N, alpha, pi, directed)

  psi <- c(-.5,0.01)
  adjMatrix <- missSBM::sample(sbm$adjacencyMatrix, "degree", psi)
  sampledNet <- missSBM:::sampledNetwork$new(adjMatrix)

  Z0 <- missSBM:::clustering_indicator(sbm$memberships)
  fittedSampling <- missSBM:::degreeSampling_fit$new(sampledNet, Z0, sbm$connectParam)

  # tolerance <- 1 ## not expected good after one iterate
  # expect_lt(sum((fittedSampling$parameters - psi)^2), tolerance)
  expect_equal(fittedSampling$df, length(psi))
  expect_equal(fittedSampling$penalty, log(N) * length(psi))
  expect_lt(fittedSampling$vExpec, 0)
})

