context("test-inference-samplings-with-covariates")

set.seed(178303)
### A SBM model : ###
N <- 300
Q <- 3
alpha <- rep(1,Q)/Q                     # mixture parameter
pi <- diag(.45,Q) + .05                 # connectivity matrix
directed <- FALSE

### Draw a SBM model (Bernoulli, undirected) with covariates
M <- 10
covarMatrix <- matrix(rnorm(N*M,mean = 0, sd = 1), N, M)
covarParam  <- rnorm(M,0,1)

mySBM <- simulateSBM(N, alpha, pi, directed, covarMatrix, covarParam)
A <- mySBM$adjMatrix

test_that("Parameter estimation in dyad-centered sampling", {

  sampledNet <- sampleNetwork(A, "dyad", covarParam, covarMatrix = covarMatrix)

  fittedSampling <- missSBM:::dyadSampling_fit_covariates$new(sampledNet, mySBM$covarArray)
  expect_is(fittedSampling, "dyadSampling_fit_covariates")
  expect_true(all(fittedSampling$prob_obs > 0, fittedSampling$prob_obs < 1))

  tolerance <- 1e-2
  expect_lt(sum((fittedSampling$parameters - covarParam)^2), tolerance)
  expect_equal(fittedSampling$df, length(covarParam))
  expect_equal(fittedSampling$penalty, log(N * (N - 1)/2) * length(covarParam))
  expect_lt(fittedSampling$vExpec, 0)
})

test_that("Parameter estimation in node-centered sampling", {

  sampledNet <- sampleNetwork(A, "node", covarParam, covarMatrix = covarMatrix)

  fittedSampling <- missSBM:::nodeSampling_fit_covariates$new(sampledNet, mySBM$covarMatrix)
  expect_is(fittedSampling, "nodeSampling_fit_covariates")
  expect_true(all(fittedSampling$prob_obs > 0, fittedSampling$prob_obs < 1))

  tolerance <- 1e-1
  expect_lt(sum((fittedSampling$parameters - covarParam)^2)/length(covarParam), tolerance)
  expect_equal(fittedSampling$df, length(covarParam))
  expect_equal(fittedSampling$penalty, log(N) * length(covarParam))
  expect_lt(fittedSampling$vExpec, 0)
})
