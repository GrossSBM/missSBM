context("test-inference-samplings")

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

test_that("Parameter estimation in dyad-centered sampling", {
  psi <- 0.1
  sampledNet <- samplingSBM(A, "dyad", psi)
  fittedSampling <- missSBM:::dyadSampling_fit$new(sampledNet)
  expect_is(fittedSampling, c("dyadSampling_fit", "networkSamplingDyads_fit", "networkSampling", "R6"))

  tolerance <- 1e-2
  expect_lt(abs(fittedSampling$parameters - psi), tolerance)

  expect_equal(fittedSampling$df, length(psi))
  expect_equal(fittedSampling$penalty, log(N * (N - 1)/2) * length(psi))
  expect_lt(fittedSampling$vExpec, 0)
})

test_that("Parameter estimation in node-centered sampling", {
  psi <- 0.1
  sampledNet <- samplingSBM(A, "node", psi)
  fittedSampling <- missSBM:::nodeSampling_fit$new(sampledNet)
  expect_is(fittedSampling, c("nodeSampling_fit", "networkSamplingNodes_fit", "networkSampling", "R6"))

  tolerance <- 1e-2
  expect_lt(abs(fittedSampling$parameters - psi), tolerance)
  expect_equal(fittedSampling$df, length(psi))
  expect_equal(fittedSampling$penalty, log(N) * length(psi))
  expect_lt(fittedSampling$vExpec, 0)
})


test_that("Parameter estimation in double-standard sampling", {
  psi <- c(0.3, 0.6)
  sampledNet <- samplingSBM(A, "double_standard", psi)
  fittedSampling <- missSBM:::doubleStandardSampling_fit$new(sampledNet)
  expect_is(fittedSampling, c("doubleStandardSampling_fit", "networkSamplingDyads_fit", "networkSampling", "R6"))

  tolerance <- .5 ## not expected good after one iterate
  expect_lt(sqrt(sum((fittedSampling$parameters - psi)^2)), tolerance)
  expect_equal(fittedSampling$df, length(psi))
  expect_equal(fittedSampling$penalty, log(N * (N - 1)/2) * length(psi))
  expect_lt(fittedSampling$vExpec, 0)
})

test_that("Parameter estimation in block sampling", {
  psi <- c(.1, .2, .3, .5, .7)
  sampledNet <- samplingSBM(A, "block", psi, NULL, mySBM$memberships)
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
  sampledNet <- samplingSBM(A,"degree", psi)
  Z0 <- missSBM:::clustering_indicator(mySBM$memberships)
  expect_warning(fittedSampling <- missSBM:::degreeSampling_fit$new(sampledNet, Z0, mySBM$connectParam))

  tolerance <- .5 ## not expected good after one iterate
  expect_lt(sqrt(sum((fittedSampling$parameters - psi)^2)), tolerance)
  expect_equal(fittedSampling$df, length(psi))
  expect_equal(fittedSampling$penalty, log(N) * length(psi))
  expect_lt(fittedSampling$vExpec, 0)
})

