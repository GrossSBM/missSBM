context("test network samplers (Class networkSampling_fit and chidren)")

set.seed(178303)
### A SBM model : ###
N <- 400
Q <- 5
alpha <- rep(1,Q)/Q                       # mixture parameter
pi <- diag(.45,Q) + .05                   # connectivity matrix
directed <- FALSE                         # if the network is directed or not

### Draw a SBM undirected model
mySBM <- simulateSBM(N, alpha, pi, directed)
A <- mySBM$adjMatrix

### Draw a SBM model (Bernoulli, undirected) with covariates
M <- 10
covarMatrix <- matrix(rnorm(N*M,mean = 0, sd = 1), N, M)
covarParam  <- rnorm(M,0,1)
mySBM_cov <- simulateSBM(N, alpha, pi, directed, covarMatrix, covarParam)
A_cov <- mySBM$adjMatrix

test_that("Parameter estimation in dyad-centered sampling", {
  psi <- 0.1
  sampledNet <- sampleNetwork(A, "dyad", psi)
  fittedSampling <- missSBM:::dyadSampling_fit$new(sampledNet)
  expect_is(fittedSampling, c("dyadSampling_fit", "networkSamplingDyads_fit", "networkSampling", "R6"))

  tolerance <- 1e-2
  expect_lt(abs(fittedSampling$parameters - psi), tolerance)

  expect_equal(fittedSampling$df, length(psi))
  expect_equal(fittedSampling$penalty, log(N * (N - 1)/2) * length(psi))
  expect_lt(fittedSampling$vExpec, 0)
})

test_that("Parameter estimation in dyad-centered sampling with covariates", {

  sampledNet <- sampleNetwork(A_cov, "dyad", covarParam, covarMatrix = covarMatrix)

  fittedSampling <- missSBM:::dyadSampling_fit_covariates$new(sampledNet, mySBM_cov$covarArray)
  expect_is(fittedSampling, "dyadSampling_fit_covariates")
  expect_true(all(fittedSampling$prob_obs > 0, fittedSampling$prob_obs < 1))

  tolerance <- 1e-2
  expect_lt(sum((fittedSampling$parameters - covarParam)^2), tolerance)
  expect_equal(fittedSampling$df, length(covarParam))
  expect_equal(fittedSampling$penalty, log(N * (N - 1)/2) * length(covarParam))
  expect_lt(fittedSampling$vExpec, 0)
})

test_that("Parameter estimation in node-centered sampling", {
  psi <- 0.1
  sampledNet <- sampleNetwork(A, "node", psi)
  fittedSampling <- missSBM:::nodeSampling_fit$new(sampledNet)
  expect_is(fittedSampling, c("nodeSampling_fit", "networkSamplingNodes_fit", "networkSampling", "R6"))

  tolerance <- 1e-1
  expect_lt(abs(fittedSampling$parameters - psi), tolerance)
  expect_equal(fittedSampling$df, length(psi))
  expect_equal(fittedSampling$penalty, log(N) * length(psi))
  expect_lt(fittedSampling$vExpec, 0)
})

test_that("Parameter estimation in node-centered sampling with covariates", {

  sampledNet <- sampleNetwork(A_cov, "node", covarParam, covarMatrix = covarMatrix)

  fittedSampling <- missSBM:::nodeSampling_fit_covariates$new(sampledNet, covarMatrix)
  expect_is(fittedSampling, "nodeSampling_fit_covariates")
  expect_true(all(fittedSampling$prob_obs > 0, fittedSampling$prob_obs < 1))

  tolerance <- .2
  expect_lt(sum((fittedSampling$parameters - covarParam)^2)/length(covarParam), tolerance)
  expect_equal(fittedSampling$df, length(covarParam))
  expect_equal(fittedSampling$penalty, log(N) * length(covarParam))
  expect_lt(fittedSampling$vExpec, 0)
})

test_that("Parameter estimation in double-standard sampling", {
  psi <- c(0.3, 0.6)
  sampledNet <- sampleNetwork(A, "double-standard", psi)
  fittedSampling <- missSBM:::doubleStandardSampling_fit$new(sampledNet)
  expect_is(fittedSampling, c("doubleStandardSampling_fit", "networkSamplingDyads_fit", "networkSampling", "R6"))

  tolerance <- .5 ## not expected good after one iterate
  expect_lt(sum((fittedSampling$parameters - psi)^2), tolerance)
  expect_equal(fittedSampling$df, length(psi))
  expect_equal(fittedSampling$penalty, log(N * (N - 1)/2) * length(psi))
  expect_lt(fittedSampling$vExpec, 0)
})

test_that("Parameter estimation in block sampling", {
  psi <- c(.1, .2, .3, .5, .7)
  sampledNet <- sampleNetwork(A, "block-node", psi, mySBM$memberships)
  Z0 <- missSBM:::clustering_indicator(mySBM$memberships)
  fittedSampling <- missSBM:::blockSampling_fit$new(sampledNet, Z0)
  expect_is(fittedSampling, c("blockSampling_fit", "networkSamplingNodes_fit", "networkSampling", "R6"))

  tolerance <- .5 ## not expected good after one iterate
  expect_lt(sqrt(sum((fittedSampling$parameters - psi)^2)), tolerance)
  expect_equal(fittedSampling$df, length(psi))
  expect_equal(fittedSampling$penalty, log(N) * length(psi))
  expect_lt(fittedSampling$vExpec, 0)
})

test_that("Parameter estimation in degree sampling", {
  psi <- c(-.5,0.01)
  sampledNet <- sampleNetwork(A,"degree", psi)
  Z0 <- missSBM:::clustering_indicator(mySBM$memberships)
  fittedSampling <- missSBM:::degreeSampling_fit$new(sampledNet, Z0, mySBM$connectParam)

  # tolerance <- 1 ## not expected good after one iterate
  # expect_lt(sum((fittedSampling$parameters - psi)^2), tolerance)
  expect_equal(fittedSampling$df, length(psi))
  expect_equal(fittedSampling$penalty, log(N) * length(psi))
  expect_lt(fittedSampling$vExpec, 0)
})

