context("test-network-data-samplings")

set.seed(178303)
### A SBM model : ###
N <- 500
Q <- 3
alpha <- rep(1,Q)/Q                     # mixture parameter
pi <- diag(.45,Q) + .05                 # connectivity matrix
directed <- FALSE

### Draw a SBM model (Bernoulli, undirected)
mySBM <- simulateSBM(N, alpha, pi, directed)
A <- mySBM$adjMatrix

### Draw a SBM model (Bernoulli, undirected) with covariates
M <- 10
covarMatrix <- matrix(rnorm(N*M,mean = 0, sd = 1), N, M)
covarParam  <- rnorm(M,0,1)
mySBM_cov <- simulateSBM(N, alpha, pi, directed, covarMatrix, covarParam)
A_cov <- mySBM_cov$adjMatrix

test_that("Consistency of dyad-centered sampling", {

  ## testing the formatting of the output
  dyad  <- samplingSBM(A, "dyad", .1)
  expect_is(dyad, "sampledNetwork", "R6")
  expect_equal(dim(dyad$adjMatrix), dim(A))
  ## expect error if psi is negative
  expect_error(samplingSBM(A, "dyad", -.1))

  # With dyad sampling, psi is the probability of sampling a dyad
  # The samplign rate is very well controlled
  for (psi in c(.05, .1, .25, .5)) {
    dyad <- samplingSBM(A, "dyad", psi)
    expect_lt(abs(dyad$samplingRate - psi), psi/10)
  }

  ## with covariates
  psi <- runif(M, -5, 5)
  dyad <- samplingSBM(A_cov, "dyad", psi, covarMatrix = covarMatrix)
  expect_is(dyad, "sampledNetwork", "R6")
  expect_equal(dim(dyad$adjMatrix), dim(A_cov))

})

test_that("Consistency of node-centered network sampling", {

  node  <- samplingSBM(A, "node", .1)
  expect_is(node, "sampledNetwork", "R6")
  expect_equal(dim(node$adjMatrix), dim(A))
  expect_error(samplingSBM(A, "node", -.1))

  # With node sampling, psi is the probability of sampling a node
  # The expected samplign rate is psi * (2-psi)
  for (psi in c(.05, .1, .25, .5)) {
    node <- samplingSBM(A, "node", psi)
    expect_lt(abs(node$samplingRate - psi * (2 - psi)), .1)
  }

  ## with covariates
  psi <- runif(M, -5, 5)
  node <- samplingSBM(A, "node", psi, covarMatrix = covarMatrix)
  expect_is(node, "sampledNetwork", "R6")
  expect_equal(dim(node$adjMatrix), dim(A))

})

test_that("Consistency of block-node network sampling", {

  block <- samplingSBM(A, "block-node", c(.1, .2, .7), clusters = mySBM$memberships)
  expect_is(block, "sampledNetwork", "R6")
  expect_equal(dim(block$adjMatrix), dim(A))
  ## error if psi is not of size Q
  expect_error(samplingSBM(A, "block-node", c(.1, .2), clusters = mySBM$memberships))
  ## error if no clustering is given
  expect_error(samplingSBM(A, "block-node", c(.1, .2, .7)))

})

test_that("Consistency of block-node network sampling", {

  block <- samplingSBM(A, "block-dyad", mySBM$connectParam, clusters = mySBM$memberships)
  expect_is(block, "sampledNetwork", "R6")
  expect_equal(dim(block$adjMatrix), dim(A))
  ## error if psi is not of size Q x Q
  expect_error(samplingSBM(A, "block-dyad", c(.1, .2, .7), clusters = mySBM$memberships))
  ## error if psi is not probabilities
  expect_error(samplingSBM(A, "block-dyad", -mySBM$connectParam, clusters = mySBM$memberships))
  ## error if no clustering is given
  expect_error(samplingSBM(A, "block-dyad", mySBM$connectParam))
})

test_that("Consistency of double-standard sampling", {

  double_standard <- samplingSBM(A,"double-standard", c(0.1, 0.5))
  expect_is(double_standard, "sampledNetwork", "R6")
  expect_equal(dim(double_standard$adjMatrix), dim(A))
  expect_error(samplingSBM(A, "double-standard", c(-0.1, 0.5)))
  expect_error(samplingSBM(A, "double-standard", c(0.1, -0.5)))
  expect_error(samplingSBM(A, "double-standard", c(-0.1, -0.5)))
})


test_that("Consistency of degree network sampling", {

  degree <- samplingSBM(A,"degree", c(0.01,0.01))
  expect_is(degree, "sampledNetwork", "R6")
  expect_equal(dim(degree$adjMatrix), dim(A))

})

