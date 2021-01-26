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
  adjMatrix  <- missSBM::observeNetwork(sbm$networkData, "dyad", .1)
  dyad <- missSBM:::partlyObservedNetwork$new(adjMatrix)

  expect_is(dyad, "partlyObservedNetwork", "R6")
  expect_lte(dyad$samplingRate, 1)
  expect_gte(dyad$samplingRate, 0)
  expect_gte(dyad$nbNodes, N)
  expect_gte(dyad$nbDyads, N * (N - 1)/2)
  expect_equal(dyad$is_directed, FALSE)
  expect_equal(dim(dyad$networkData), dim(sbm$networkData))

  ## expect error if psi is negative
  expect_error(missSBM::observeNetwork(sbm$networkData, "dyad", -.1))

  # With dyad sampling, psi is the probability of sampling a dyad
  # The samplign rate is very well controlled
  for (psi in c(.1, .25, .4)) {
    adjMatrix <-missSBM::observeNetwork(sbm$networkData, "dyad", psi)
    dyad  <- missSBM:::partlyObservedNetwork$new(adjMatrix)
    expect_lt(abs(dyad$samplingRate - psi), psi/10)
  }

  ## with covariates
  psi <- runif(M, -5, 5)
  adjMatrix <- missSBM::observeNetwork(sbm_cov_dyad$networkData, "covar-dyad", psi, covariates = covariates_dyad)
  ## Prepare network data for estimation with missing data
  dyad <- missSBM:::partlyObservedNetwork$new(adjMatrix, covariates_dyad, missSBM:::l1_similarity)
  expect_is(dyad, "partlyObservedNetwork", "R6")
  expect_equal(dim(dyad$networkData), dim(sbm_cov_dyad$networkData))

})

test_that("Consistency of node-centered network sampling", {

  adjMatrix <- missSBM::observeNetwork(sbm$networkData, "node", .1)
  node <- missSBM:::partlyObservedNetwork$new(adjMatrix)

  expect_is(node, "partlyObservedNetwork", "R6")
  expect_lte(node$samplingRate, 1)
  expect_gte(node$samplingRate, 0)
  expect_gte(node$nbNodes, N)
  expect_gte(node$nbDyads, N * (N-1)/2)
  expect_equal(node$is_directed, FALSE)
  expect_equal(dim(node$networkData), dim(sbm$networkData))
  expect_error(missSBM::observeNetwork(sbm$adjacency, "node", -.1))

  # With node sampling, psi is the probability of sampling a node
  # The expected samplign rate is psi * (2-psi)
  for (psi in c(.05, .1, .25, .5)) {
    adjMatrix <- missSBM::observeNetwork(sbm$networkData, "node", psi)
    node <- missSBM:::partlyObservedNetwork$new(adjMatrix)
    expect_lt(abs(node$samplingRate - psi * (2 - psi)), .1)
  }

  ## with covariates
  psi <- runif(M, -5, 5)
  adjMatrix <- missSBM::observeNetwork(sbm_cov_node$networkData, "covar-node", psi, covariates = covariates_node)
  node <- missSBM:::partlyObservedNetwork$new(adjMatrix, covariates_node, missSBM:::l1_similarity)
  expect_is(node, "partlyObservedNetwork", "R6")
  expect_equal(dim(node$networkData), dim(sbm_cov_node$networkData))

})

test_that("Consistency of block-node network sampling", {

  adjMatrix <- missSBM::observeNetwork(sbm$networkData, "block-node", c(.1, .2, .7), clusters = sbm$memberships)
  block <- missSBM:::partlyObservedNetwork$new(adjMatrix)
  expect_is(block, "partlyObservedNetwork", "R6")
  expect_lte(block$samplingRate, 1)
  expect_gte(block$samplingRate, 0)
  expect_gte(block$nbNodes, N)
  expect_gte(block$nbDyads, N * (N - 1)/2)
  expect_equal(block$is_directed, FALSE)
  expect_equal(dim(block$networkData), dim(sbm$networkData))
  ## error if psi is not of size Q
  expect_error(missSBM::observeNetwork(sbm$networkData, "block-node", c(.1, .2), clusters = sbm$memberships))
  ## error if no clustering is given
  expect_error(missSBM::observeNetwork(sbm$adjacency, "block-node", c(.1, .2, .7)))

})

test_that("Consistency of block-node network sampling", {

  adjMatrix <- missSBM::observeNetwork(sbm$networkData, "block-dyad", sbm$connectParam$mean, clusters = sbm$memberships)
  block <- missSBM:::partlyObservedNetwork$new(adjMatrix)
  expect_is(block, "partlyObservedNetwork", "R6")
  expect_lte(block$samplingRate, 1)
  expect_gte(block$samplingRate, 0)
  expect_gte(block$nbNodes, N)
  expect_gte(block$nbDyads, N * (N - 1)/2)
  expect_equal(block$is_directed, FALSE)
  expect_equal(dim(block$networkData), dim(sbm$networkData))
  ## error if psi is not of size Q x Q
  expect_error(missSBM::observeNetwork(sbm$networkData, "block-dyad", c(.1, .2, .7), clusters = sbm$memberships))
  ## error if psi is not probabilities
  expect_error(missSBM::observeNetwork(sbm$networkData, "block-dyad", -sbm$connectParam, clusters = sbm$memberships))
  ## error if no clustering is given
  expect_error(missSBM::observeNetwork(sbm$networkData, "block-dyad", sbm$connectParam))
})

test_that("Consistency of double-standard sampling", {

  adjMatrix <- missSBM::observeNetwork(sbm$networkData,"double-standard", c(0.1, 0.5))
  double_standard <- missSBM:::partlyObservedNetwork$new(adjMatrix)

  expect_is(double_standard, "partlyObservedNetwork", "R6")
  expect_lte(double_standard$samplingRate, 1)
  expect_gte(double_standard$samplingRate, 0)
  expect_gte(double_standard$nbNodes, N)
  expect_gte(double_standard$nbDyads, N * (N - 1)/2)
  expect_equal(double_standard$is_directed, FALSE)
  expect_equal(dim(double_standard$networkData), dim(sbm$networkData))
  expect_error(missSBM::observeNetwork(sbm$adjacency, "double-standard", c(-0.1, 0.5)))
  expect_error(missSBM::observeNetwork(sbm$adjacency, "double-standard", c(0.1, -0.5)))
  expect_error(missSBM::observeNetwork(sbm$adjacency, "double-standard", c(-0.1, -0.5)))
})


test_that("Consistency of degree network sampling", {

  adjMatrix <- missSBM::observeNetwork(sbm$networkData,"degree", c(0.01,0.01))
  degree <- missSBM:::partlyObservedNetwork$new(adjMatrix)
  expect_is(degree, "partlyObservedNetwork", "R6")
  expect_lte(degree$samplingRate, 1)
  expect_gte(degree$samplingRate, 0)
  expect_gte(degree$nbNodes, N)
  expect_gte(degree$nbDyads, N * (N-1)/2)
  expect_equal(degree$is_directed, FALSE)
  expect_equal(dim(degree$networkData), dim(sbm$networkData))

})

