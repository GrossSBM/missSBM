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

### Draw a SBM model (Bernoulli, undirected) with covariates
M <- 10
covariates_node <- replicate(M, rnorm(N,mean = 0, sd = 1), simplify = FALSE)
covarMatrix <- simplify2array(covariates_node)
covarArray  <- missSBM:::getCovarArray(covarMatrix, missSBM:::l1_similarity)
covariates_dyad <- lapply(seq(dim(covarArray)[3]), function(x) covarArray[ , , x])

covarParam  <- rnorm(M, 0, 1)
sbm_cov_dyad <- missSBM::simulate(N, alpha, gamma, directed, covariates_dyad, covarParam)
sbm_cov_node <- missSBM::simulate(N, alpha, gamma, directed, covariates_dyad, covarParam)

test_that("Consistency of dyad-centered sampling", {

  ## testing the formatting of the output
  dyad  <- missSBM::sample(sbm$adjMatrix, "dyad", .1)
  expect_is(dyad, "sampledNetwork", "R6")
  expect_lte(dyad$samplingRate, 1)
  expect_gte(dyad$samplingRate, 0)
  expect_gte(dyad$nNodes, N)
  expect_gte(dyad$nDyads, N * (N - 1)/2)
  expect_equal(dyad$is_directed, directed)
  expect_equal(dim(dyad$adjMatrix), dim(sbm$adjMatrix))

  ## expect error if psi is negative
  expect_error(missSBM::sample(sbm$adjMatrix, "dyad", -.1))

  # With dyad sampling, psi is the probability of sampling a dyad
  # The samplign rate is very well controlled
  for (psi in c(.05, .1, .25, .5)) {
    dyad <- missSBM::sample(sbm$adjMatrix, "dyad", psi)
    expect_lt(abs(dyad$samplingRate - psi), psi/10)
  }

  ## with covariates
  psi <- runif(M, -5, 5)
  dyad <- missSBM::sample(sbm_cov_dyad$adjMatrix, "dyad", psi, covariates = covariates_dyad)
  expect_is(dyad, "sampledNetwork", "R6")
  expect_equal(dim(dyad$adjMatrix), dim(sbm_cov_dyad$adjMatrix))

})

test_that("Consistency of node-centered network sampling", {

  node  <- missSBM::sample(sbm$adjMatrix, "node", .1)
  expect_is(node, "sampledNetwork", "R6")
  expect_lte(node$samplingRate, 1)
  expect_gte(node$samplingRate, 0)
  expect_gte(node$nNodes, N)
  expect_gte(node$nDyads, N * (N-1)/2)
  expect_equal(node$is_directed, directed)
  expect_equal(dim(node$adjMatrix), dim(sbm$adjMatrix))
  expect_error(missSBM::sample(sbm$adjMatrix, "node", -.1))

  # With node sampling, psi is the probability of sampling a node
  # The expected samplign rate is psi * (2-psi)
  for (psi in c(.05, .1, .25, .5)) {
    node <- missSBM::sample(sbm$adjMatrix, "node", psi)
    expect_lt(abs(node$samplingRate - psi * (2 - psi)), .1)
  }

  ## with covariates
  psi <- runif(M, -5, 5)
  node <- missSBM::sample(sbm_cov_node$adjMatrix, "node", psi, covariates = covariates_node)
  expect_is(node, "sampledNetwork", "R6")
  expect_equal(dim(node$adjMatrix), dim(sbm_cov_node$adjMatrix))

})

test_that("Consistency of block-node network sampling", {

  block <- missSBM::sample(sbm$adjMatrix, "block-node", c(.1, .2, .7), clusters = sbm$memberships)
  expect_is(block, "sampledNetwork", "R6")
  expect_lte(block$samplingRate, 1)
  expect_gte(block$samplingRate, 0)
  expect_gte(block$nNodes, N)
  expect_gte(block$nDyads, N * (N - 1)/2)
  expect_equal(block$is_directed, directed)
  expect_equal(dim(block$adjMatrix), dim(sbm$adjMatrix))
  ## error if psi is not of size Q
  expect_error(missSBM::sample(sbm$adjMatrix, "block-node", c(.1, .2), clusters = sbm$memberships))
  ## error if no clustering is given
  expect_error(missSBM::sample(sbm$adjMatrix, "block-node", c(.1, .2, .7)))

})

test_that("Consistency of block-node network sampling", {

  block <- missSBM::sample(sbm$adjMatrix, "block-dyad", sbm$connectParam, clusters = sbm$memberships)
  expect_is(block, "sampledNetwork", "R6")
  expect_lte(block$samplingRate, 1)
  expect_gte(block$samplingRate, 0)
  expect_gte(block$nNodes, N)
  expect_gte(block$nDyads, N * (N - 1)/2)
  expect_equal(block$is_directed, directed)
  expect_equal(dim(block$adjMatrix), dim(sbm$adjMatrix))
  ## error if psi is not of size Q x Q
  expect_error(missSBM::sample(sbm$adjMatrix, "block-dyad", c(.1, .2, .7), clusters = sbm$memberships))
  ## error if psi is not probabilities
  expect_error(missSBM::sample(sbm$adjMatrix, "block-dyad", -sbm$connectParam, clusters = sbm$memberships))
  ## error if no clustering is given
  expect_error(missSBM::sample(sbm$adjMatrix, "block-dyad", sbm$connectParam))
})

test_that("Consistency of double-standard sampling", {

  double_standard <- missSBM::sample(sbm$adjMatrix,"double-standard", c(0.1, 0.5))
  expect_is(double_standard, "sampledNetwork", "R6")
  expect_lte(double_standard$samplingRate, 1)
  expect_gte(double_standard$samplingRate, 0)
  expect_gte(double_standard$nNodes, N)
  expect_gte(double_standard$nDyads, N * (N - 1)/2)
  expect_equal(double_standard$is_directed, directed)
  expect_equal(dim(double_standard$adjMatrix), dim(sbm$adjMatrix))
  expect_error(missSBM::sample(sbm$adjMatrix, "double-standard", c(-0.1, 0.5)))
  expect_error(missSBM::sample(sbm$adjMatrix, "double-standard", c(0.1, -0.5)))
  expect_error(missSBM::sample(sbm$adjMatrix, "double-standard", c(-0.1, -0.5)))
})


test_that("Consistency of degree network sampling", {

  degree <- missSBM::sample(sbm$adjMatrix,"degree", c(0.01,0.01))
  expect_is(degree, "sampledNetwork", "R6")
  expect_lte(degree$samplingRate, 1)
  expect_gte(degree$samplingRate, 0)
  expect_gte(degree$nNodes, N)
  expect_gte(degree$nDyads, N * (N-1)/2)
  expect_equal(degree$is_directed, directed)
  expect_equal(dim(degree$adjMatrix), dim(sbm$adjMatrix))

})

