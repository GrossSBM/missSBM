context("test-sbm-sampler")


test_that("SBM sampler without covariates", {

  set.seed(178303)
  ### A SBM model : ###
  N <- 400
  Q <- 5
  alpha <- rep(1,Q)/Q                     # mixture parameter
  pi <- diag(.45,Q) + .05                 # connectivity matrix
  directed <- FALSE

  mySBM <- missSBM:::SBM_sampler$new(directed, N, alpha, pi)
  expect_null(mySBM$adjMatrix)
  expect_null(mySBM$blocks)
  expect_error(mySBM$memberships)
  expect_equal(mySBM$connectParam, pi)
  expect_equal(mySBM$mixtureParam, alpha)
  expect_null(mySBM$covarParam)
  expect_null(mySBM$covarArray)
  expect_false(mySBM$hasCovariates)
  expect_equal(mySBM$direction, "undirected")

  mySBM$rBlocks()
  expect_equal(length(mySBM$memberships), N)
  expect_equal(dim(mySBM$blocks), c(N, Q))
  expect_equal(length(unique(mySBM$memberships)), Q)

  mySBM$rAdjMatrix()
  expect_equal(dim(mySBM$adjMatrix), c(N, N))
  expect_true(isSymmetric(mySBM$adjMatrix))
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
  covarMatrix <- matrix(rnorm(N*M,mean = 0, sd = 1), N, M)
  covarParam  <- rnorm(M,0,1)
  covarArray <- missSBM:::getCovarArray(covarMatrix, missSBM:::l1_similarity)
  mySBM <- missSBM:::SBM_sampler$new(directed, N, alpha, pi, covarParam, covarArray)
  expect_null(mySBM$adjMatrix)
  expect_null(mySBM$blocks)
  expect_error(mySBM$memberships)
  expect_equal(mySBM$connectParam, pi)
  expect_equal(mySBM$mixtureParam, alpha)
  expect_true(mySBM$hasCovariates)
  expect_equal(mySBM$covarParam, covarParam)
  expect_equal(dim(mySBM$covarArray) , c(N, N, M))
  expect_equal(mySBM$direction, "undirected")

  mySBM$rBlocks()
  expect_equal(length(mySBM$memberships), N)
  expect_equal(dim(mySBM$blocks), c(N, Q))
  expect_equal(length(unique(mySBM$memberships)), Q)

  mySBM$rAdjMatrix()
  expect_equal(dim(mySBM$adjMatrix), c(N, N))
  expect_true(isSymmetric(mySBM$adjMatrix))

})
