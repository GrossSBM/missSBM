context("testing network samplers (top-level function missSBM::sample)")

set.seed(178303)
### A SBM model : ###
N <- 500
Q <- 3
alpha <- rep(1, Q)/Q                     # mixture parameter
pi <- diag(.45, Q, Q) + .05                 # connectivity matrix
gamma <- missSBM:::logit(pi)
directed <- FALSE

### Draw a SBM model (Bernoulli, undirected)
sbm <- missSBM::simulate(N, alpha, pi, directed)
A <- sbm$adjMatrix

### Draw a SBM model (Bernoulli, undirected) with covariates
M <- 10
covariates  <- replicate(M, rnorm(N,mean = 0, sd = 1), simplify = FALSE)
covarMatrix <- simplify2array(covariates)
covarArray  <- missSBM:::getCovarArray(covarMatrix, missSBM:::l1_similarity)
covarParam  <- rnorm(M, 0, 1)
sbm_cov <- missSBM::simulate(N, alpha, gamma, directed, covariates, covarParam)
A_cov <- sbm_cov$adjMatrix

test_that("Consistency of dyad-centered sampling", {

  ## testing the formatting of the output
  dyad  <- missSBM::sample(A, "dyad", .1)
  expect_is(dyad, "sampledNetwork", "R6")
  expect_lte(dyad$samplingRate, 1)
  expect_gte(dyad$samplingRate, 0)
  expect_gte(dyad$nNodes, N)
  expect_gte(dyad$nDyads, N * (N-1)/2)
  expect_equal(dyad$is_directed, directed)
  expect_equal(dim(dyad$adjMatrix), dim(A))

  ## expect error if psi is negative
  expect_error(missSBM::sample(A, "dyad", -.1))

  # With dyad sampling, psi is the probability of sampling a dyad
  # The samplign rate is very well controlled
  for (psi in c(.05, .1, .25, .5)) {
    dyad <- missSBM::sample(A, "dyad", psi)
    expect_lt(abs(dyad$samplingRate - psi), psi/10)
  }

  ## with covariates
  psi <- runif(M, -5, 5)
  dyad <- missSBM::sample(A_cov, "dyad", psi, covariates = covariates)
  expect_is(dyad, "sampledNetwork", "R6")
  expect_equal(dim(dyad$adjMatrix), dim(A_cov))

})

test_that("Consistency of node-centered network sampling", {

  node  <- missSBM::sample(A, "node", .1)
  expect_is(node, "sampledNetwork", "R6")
  expect_lte(node$samplingRate, 1)
  expect_gte(node$samplingRate, 0)
  expect_gte(node$nNodes, N)
  expect_gte(node$nDyads, N * (N-1)/2)
  expect_equal(node$is_directed, directed)
  expect_equal(dim(node$adjMatrix), dim(A))
  expect_error(missSBM::sample(A, "node", -.1))

  # With node sampling, psi is the probability of sampling a node
  # The expected samplign rate is psi * (2-psi)
  for (psi in c(.05, .1, .25, .5)) {
    node <- missSBM::sample(A, "node", psi)
    expect_lt(abs(node$samplingRate - psi * (2 - psi)), .1)
  }

  ## with covariates
  psi <- runif(M, -5, 5)
  node <- missSBM::sample(A, "node", psi, covariates = covarMatrix)
  expect_is(node, "sampledNetwork", "R6")
  expect_equal(dim(node$adjMatrix), dim(A))

})

test_that("Consistency of block-node network sampling", {

  block <- missSBM::sample(A, "block-node", c(.1, .2, .7), clusters = sbm$memberships)
  expect_is(block, "sampledNetwork", "R6")
  expect_lte(block$samplingRate, 1)
  expect_gte(block$samplingRate, 0)
  expect_gte(block$nNodes, N)
  expect_gte(block$nDyads, N * (N-1)/2)
  expect_equal(block$is_directed, directed)
  expect_equal(dim(block$adjMatrix), dim(A))
  ## error if psi is not of size Q
  expect_error(missSBM::sample(A, "block-node", c(.1, .2), clusters = sbm$memberships))
  ## error if no clustering is given
  expect_error(missSBM::sample(A, "block-node", c(.1, .2, .7)))

})

test_that("Consistency of block-node network sampling", {

  block <- missSBM::sample(A, "block-dyad", sbm$connectParam, clusters = sbm$memberships)
  expect_is(block, "sampledNetwork", "R6")
  expect_lte(block$samplingRate, 1)
  expect_gte(block$samplingRate, 0)
  expect_gte(block$nNodes, N)
  expect_gte(block$nDyads, N * (N-1)/2)
  expect_equal(block$is_directed, directed)
  expect_equal(dim(block$adjMatrix), dim(A))
  ## error if psi is not of size Q x Q
  expect_error(missSBM::sample(A, "block-dyad", c(.1, .2, .7), clusters = sbm$memberships))
  ## error if psi is not probabilities
  expect_error(missSBM::sample(A, "block-dyad", -sbm$connectParam, clusters = sbm$memberships))
  ## error if no clustering is given
  expect_error(missSBM::sample(A, "block-dyad", sbm$connectParam))
})

test_that("Consistency of double-standard sampling", {

  double_standard <- missSBM::sample(A,"double-standard", c(0.1, 0.5))
  expect_is(double_standard, "sampledNetwork", "R6")
  expect_lte(double_standard$samplingRate, 1)
  expect_gte(double_standard$samplingRate, 0)
  expect_gte(double_standard$nNodes, N)
  expect_gte(double_standard$nDyads, N * (N-1)/2)
  expect_equal(double_standard$is_directed, directed)
  expect_equal(dim(double_standard$adjMatrix), dim(A))
  expect_error(missSBM::sample(A, "double-standard", c(-0.1, 0.5)))
  expect_error(missSBM::sample(A, "double-standard", c(0.1, -0.5)))
  expect_error(missSBM::sample(A, "double-standard", c(-0.1, -0.5)))
})


test_that("Consistency of degree network sampling", {

  degree <- missSBM::sample(A,"degree", c(0.01,0.01))
  expect_is(degree, "sampledNetwork", "R6")
  expect_lte(degree$samplingRate, 1)
  expect_gte(degree$samplingRate, 0)
  expect_gte(degree$nNodes, N)
  expect_gte(degree$nDyads, N * (N-1)/2)
  expect_equal(degree$is_directed, directed)
  expect_equal(dim(degree$adjMatrix), dim(A))

})

