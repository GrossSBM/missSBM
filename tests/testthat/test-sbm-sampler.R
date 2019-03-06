context("test-sbm-sampler")

set.seed(178303)

test_that("SBM sampler without covariates", {

  ### A SBM model : ###
  N <- 400
  Q <- 5
  alpha <- rep(1,Q)/Q                     # mixture parameter
  pi <- diag(.45,Q) + .05                 # connectivity matrix
  directed <- FALSE

  mySBM <- missSBM:::SBM_sampler$new(directed, N, alpha, pi)
  expect_null(mySBM$adjacencyMatrix)
  expect_null(mySBM$blocks)
  expect_error(mySBM$memberships)
  expect_equal(mySBM$connectParam, pi)
  expect_equal(mySBM$mixtureParam, alpha)
  expect_null(mySBM$covarParam)
  expect_null(mySBM$covariates)
  expect_false(mySBM$has_covariates)
  expect_equal(mySBM$direction, "undirected")

  mySBM$rBlocks()
  expect_equal(length(mySBM$memberships), N)
  expect_equal(dim(mySBM$blocks), c(N, Q))
  expect_equal(length(unique(mySBM$memberships)), Q)

  mySBM$rAdjMatrix()
  expect_equal(dim(mySBM$adjacencyMatrix), c(N, N))
  expect_true(isSymmetric(mySBM$adjacencyMatrix))
})

test_that("SBM sampler with covariates", {

  set.seed(178303)
  ### A SBM model : ###
  N <- 300
  Q <- 3
  alpha <- rep(1,Q)/Q                     # mixture parameter
  pi <- diag(.45,Q) + .05                 # connectivity matrix
  directed <- FALSE

  ### Draw a SBM model (Bernoulli, undirected) with covariates
  M <- 4
  covariates <- matrix(rnorm(N*M,mean = 0, sd = 1), N, M)
  covarParam  <- rnorm(M,0,1)

  mySBM <- missSBM:::SBM_sampler$new(directed, N, alpha, pi, covariates, covarParam)
  expect_null(mySBM$adjacencyMatrix)
  expect_null(mySBM$blocks)
  expect_error(mySBM$memberships)
  expect_equal(mySBM$connectParam, pi)
  expect_equal(mySBM$mixtureParam, alpha)
  expect_true(mySBM$has_covariates)
  expect_equal(mySBM$covarParam, covarParam)
  expect_equal(dim(mySBM$covariates), c(N, N, M))
  expect_equal(mySBM$direction, "undirected")

  mySBM$rBlocks()
  expect_equal(length(mySBM$memberships), N)
  expect_equal(dim(mySBM$blocks), c(N, Q))
  expect_equal(length(unique(mySBM$memberships)), Q)

  mySBM$rAdjMatrix()
  expect_equal(dim(mySBM$adjacencyMatrix), c(N, N))
  expect_true(isSymmetric(mySBM$adjacencyMatrix))

})
