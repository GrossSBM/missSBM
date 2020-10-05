context("testing network samplers (top-level function missSBM::sample)")

set.seed(178304)
### A SBM model : ###
N <- 100
Q <- 3
pi <- rep(1, Q)/Q           # block proportion
theta <- list(mean = diag(.45,Q) + .05 ) # connectivity matrix

### Draw a SBM model (Bernoulli, undirected)
sbm <- sbm::sampleSimpleSBM(N, pi, theta)

### Draw a SBM model (Bernoulli, undirected) with covariates
M <- 2
covariates_node <- replicate(M, rnorm(N,mean = 0, sd = 1), simplify = FALSE)
covarMatrix <- simplify2array(covariates_node)
covarArray  <- missSBM:::getCovarArray(covarMatrix, missSBM:::l1_similarity)
covariates_dyad <- missSBM:::array2list(covarArray)

covarParam  <- rnorm(M, 0, 1)
sbm_cov_dyad <- sbm::sampleSimpleSBM(N, pi, theta, covariates = covariates_dyad, covariatesParam = covarParam)
sbm_cov_node <- sbm::sampleSimpleSBM(N, pi, theta, covariates = covariates_dyad, covariatesParam = covarParam)

test_that("Consistency of dyad-centered sampling", {

  ## testing the formatting of the output
  adjMatrix  <- missSBM::sample(sbm$netMatrix, "dyad", .1)
  dyad <- missSBM:::sampledNetwork$new(adjMatrix)

  expect_is(dyad, "sampledNetwork", "R6")
  expect_lte(dyad$samplingRate, 1)
  expect_gte(dyad$samplingRate, 0)
  expect_gte(dyad$nbNodes, N)
  expect_gte(dyad$nbDyads, N * (N - 1)/2)
  expect_equal(dyad$is_directed, FALSE)
  expect_equal(dim(dyad$netMatrix), dim(sbm$netMatrix))

  ## expect error if psi is negative
  expect_error(missSBM::sample(sbm$netMatrix, "dyad", -.1))

  # With dyad sampling, psi is the probability of sampling a dyad
  # The samplign rate is very well controlled
  for (psi in c(.1, .25, .4)) {
    adjMatrix <-missSBM::sample(sbm$netMatrix, "dyad", psi)
    dyad  <- missSBM:::sampledNetwork$new(adjMatrix)
    expect_lt(abs(dyad$samplingRate - psi), psi/10)
  }

  ## with covariates
  psi <- runif(M, -5, 5)
  adjMatrix <- missSBM::sample(sbm_cov_dyad$netMatrix, "covar-dyad", psi, covariates = covariates_dyad)
  ## Prepare network data for estimation with missing data
  dyad <- missSBM:::sampledNetwork$new(adjMatrix, covariates_dyad, missSBM:::l1_similarity)
  expect_is(dyad, "sampledNetwork", "R6")
  expect_equal(dim(dyad$netMatrix), dim(sbm_cov_dyad$netMatrix))

})

test_that("Consistency of node-centered network sampling", {

  adjMatrix <- missSBM::sample(sbm$netMatrix, "node", .1)
  node <- missSBM:::sampledNetwork$new(adjMatrix)

  expect_is(node, "sampledNetwork", "R6")
  expect_lte(node$samplingRate, 1)
  expect_gte(node$samplingRate, 0)
  expect_gte(node$nbNodes, N)
  expect_gte(node$nbDyads, N * (N-1)/2)
  expect_equal(node$is_directed, FALSE)
  expect_equal(dim(node$netMatrix), dim(sbm$netMatrix))
  expect_error(missSBM::sample(sbm$adjacency, "node", -.1))

  # With node sampling, psi is the probability of sampling a node
  # The expected samplign rate is psi * (2-psi)
  for (psi in c(.05, .1, .25, .5)) {
    adjMatrix <- missSBM::sample(sbm$netMatrix, "node", psi)
    node <- missSBM:::sampledNetwork$new(adjMatrix)
    expect_lt(abs(node$samplingRate - psi * (2 - psi)), .1)
  }

  ## with covariates
  psi <- runif(M, -5, 5)
  adjMatrix <- missSBM::sample(sbm_cov_node$netMatrix, "covar-node", psi, covariates = covariates_node)
  node <- missSBM:::sampledNetwork$new(adjMatrix, covariates_node, missSBM:::l1_similarity)
  expect_is(node, "sampledNetwork", "R6")
  expect_equal(dim(node$netMatrix), dim(sbm_cov_node$netMatrix))

})

test_that("Consistency of block-node network sampling", {

  adjMatrix <- missSBM::sample(sbm$netMatrix, "block-node", c(.1, .2, .7), clusters = sbm$memberships)
  block <- missSBM:::sampledNetwork$new(adjMatrix)
  expect_is(block, "sampledNetwork", "R6")
  expect_lte(block$samplingRate, 1)
  expect_gte(block$samplingRate, 0)
  expect_gte(block$nbNodes, N)
  expect_gte(block$nbDyads, N * (N - 1)/2)
  expect_equal(block$is_directed, FALSE)
  expect_equal(dim(block$netMatrix), dim(sbm$netMatrix))
  ## error if psi is not of size Q
  expect_error(missSBM::sample(sbm$netMatrix, "block-node", c(.1, .2), clusters = sbm$memberships))
  ## error if no clustering is given
  expect_error(missSBM::sample(sbm$adjacency, "block-node", c(.1, .2, .7)))

})

test_that("Consistency of block-node network sampling", {

  adjMatrix <- missSBM::sample(sbm$netMatrix, "block-dyad", sbm$connectParam$mean, clusters = sbm$memberships)
  block <- missSBM:::sampledNetwork$new(adjMatrix)
  expect_is(block, "sampledNetwork", "R6")
  expect_lte(block$samplingRate, 1)
  expect_gte(block$samplingRate, 0)
  expect_gte(block$nbNodes, N)
  expect_gte(block$nbDyads, N * (N - 1)/2)
  expect_equal(block$is_directed, FALSE)
  expect_equal(dim(block$netMatrix), dim(sbm$netMatrix))
  ## error if psi is not of size Q x Q
  expect_error(missSBM::sample(sbm$netMatrix, "block-dyad", c(.1, .2, .7), clusters = sbm$memberships))
  ## error if psi is not probabilities
  expect_error(missSBM::sample(sbm$netMatrix, "block-dyad", -sbm$connectParam, clusters = sbm$memberships))
  ## error if no clustering is given
  expect_error(missSBM::sample(sbm$netMatrix, "block-dyad", sbm$connectParam))
})

test_that("Consistency of double-standard sampling", {

  adjMatrix <- missSBM::sample(sbm$netMatrix,"double-standard", c(0.1, 0.5))
  double_standard <- missSBM:::sampledNetwork$new(adjMatrix)

  expect_is(double_standard, "sampledNetwork", "R6")
  expect_lte(double_standard$samplingRate, 1)
  expect_gte(double_standard$samplingRate, 0)
  expect_gte(double_standard$nbNodes, N)
  expect_gte(double_standard$nbDyads, N * (N - 1)/2)
  expect_equal(double_standard$is_directed, FALSE)
  expect_equal(dim(double_standard$netMatrix), dim(sbm$netMatrix))
  expect_error(missSBM::sample(sbm$adjacency, "double-standard", c(-0.1, 0.5)))
  expect_error(missSBM::sample(sbm$adjacency, "double-standard", c(0.1, -0.5)))
  expect_error(missSBM::sample(sbm$adjacency, "double-standard", c(-0.1, -0.5)))
})


test_that("Consistency of degree network sampling", {

  adjMatrix <- missSBM::sample(sbm$netMatrix,"degree", c(0.01,0.01))
  degree <- missSBM:::sampledNetwork$new(adjMatrix)
  expect_is(degree, "sampledNetwork", "R6")
  expect_lte(degree$samplingRate, 1)
  expect_gte(degree$samplingRate, 0)
  expect_gte(degree$nbNodes, N)
  expect_gte(degree$nbDyads, N * (N-1)/2)
  expect_equal(degree$is_directed, FALSE)
  expect_equal(dim(degree$netMatrix), dim(sbm$netMatrix))

})

