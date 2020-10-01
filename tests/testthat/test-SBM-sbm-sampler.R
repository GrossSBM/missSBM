context("testing SBM samplers (class SBM_sampler)")


test_that("SBM sampler without covariates", {

  set.seed(178303)
  ### A SBM model : ###
  N <- 100
  Q <- 3
  pi <- rep(1,Q)/Q                     # mixture parameter
  theta <- diag(.45,Q) + .05                 # connectivity matrix
  directed <- FALSE

  mySBM <- missSBM:::SBM_sampler$new(directed, N, pi, theta)
  expect_null(mySBM$adjacencyMatrix)
  expect_null(mySBM$indMemberships)
  expect_error(mySBM$memberships)
  expect_equal(mySBM$connectParam, theta)
  expect_equal(mySBM$blockProp, pi)
  expect_equal(mySBM$covarParam, numeric(0))
  expect_equal(mySBM$covarList, list())
  expect_false(mySBM$nbCovariates > 0)
  expect_equal(mySBM$directed, FALSE)

  mySBM$rMemberships()
  expect_equal(length(mySBM$memberships), N)
  expect_equal(dim(mySBM$indMemberships), c(N, Q))
  expect_equal(length(unique(mySBM$memberships)), Q)

  mySBM$rAdjacency()
  expect_equal(dim(mySBM$adjacencyMatrix), c(N, N))
  expect_true(isSymmetric(mySBM$adjacencyMatrix))
})

test_that("SBM sampler with covariates", {

  set.seed(178303)
  ### A SBM model : ###
  N <- 100
  Q <- 3
  pi <- rep(1,Q)/Q                     # mixture parameter
  theta <- diag(.45,Q) + .05                 # connectivity matrix
  gamma <- missSBM:::.logit(theta)
  directed <- FALSE

  ### Draw a SBM model (Bernoulli, undirected) with covariates
  M <- 2
  covarMatrix <- matrix(rnorm(N*M,mean = 0, sd = 1), N, M)
  covarParam  <- rnorm(M,0,1)
  covarArray <- missSBM:::getCovarArray(covarMatrix, missSBM:::l1_similarity)
  mySBM <- missSBM:::SBM_sampler$new(directed, N, pi, gamma, covarParam, missSBM:::array2list(covarArray))
  expect_null(mySBM$adjacencyMatrix)
  expect_null(mySBM$indMemberships)
  expect_error(mySBM$memberships)
  expect_equal(missSBM:::.logistic(mySBM$connectParam), theta)
  expect_equal(mySBM$blockProp, pi)
  expect_true(mySBM$nbCovariates > 0)
  expect_equal(mySBM$covarParam, covarParam)
  expect_equal(dim(mySBM$covarArray) , c(N, N, M))
  expect_equal(mySBM$directed, FALSE)

  mySBM$rMemberships()
  expect_equal(length(mySBM$memberships), N)
  expect_equal(dim(mySBM$indMemberships), c(N, Q))
  expect_equal(length(unique(mySBM$memberships)), Q)

  mySBM$rAdjacency()
  expect_equal(dim(mySBM$adjacencyMatrix), c(N, N))
  expect_true(isSymmetric(mySBM$adjacencyMatrix))

})
