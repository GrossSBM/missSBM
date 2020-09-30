context("test sbm fit without covariate (class SBM_fit_nocovariate)")

library(aricode)
library(blockmodels)

set.seed(178303)
### A SBM model : ###
N <- 100
Q <- 3
alpha <- rep(1, Q)/Q                     # mixture parameter
pi <- diag(.45, Q, Q) + .05                 # connectivity matrix
directed <- FALSE

### Draw a undirected SBM model
mySBM <- missSBM::simulate(N, alpha, pi, directed)
A <- mySBM$adjacencyMatrix
cl_rand <- base::sample(mySBM$memberships)
cl_spec <- missSBM:::init_clustering(A, Q, NULL, "spectral")
cl_hier <- missSBM:::init_clustering(A, Q, NULL, "hierarchical")
cl_kmns <- missSBM:::init_clustering(A, Q, NULL, "kmeans")

test_that("Creation of a SBM_fit_nocovariate", {

  mySBM_fit <- missSBM:::SBM_fit_nocovariate$new(A, cl_rand)
  expect_is(mySBM_fit, c("SBM_fit_nocovariate", "SBM_fit", "SBM", "R6"))
  expect_equal(mySBM_fit$memberships, cl_rand)
  expect_equal(mySBM_fit$df_connectParams, Q * (Q + 1)/2)
  expect_equal(mySBM_fit$df_covarParams, 0)
  expect_equal(mySBM_fit$df_blockProps, Q - 1)
  expect_equal(mySBM_fit$probMemberships, missSBM:::clustering_indicator(cl_rand))
  expect_equal(dim(mySBM_fit$connectParam), dim(mySBM$connectParam))
  expect_equal(length(mySBM_fit$blockProp), length(mySBM$blockProp))
  expect_equal(mySBM_fit$direction, "undirected")

})

test_that("Consistency of VEM of a SBM_fit_nocovariate when the number of block is given", {

  tol <- 2e-3

  ## testing all initialization
  mySBM_fit_hier <- missSBM:::SBM_fit_nocovariate$new(A, cl_hier)
  mySBM_fit_spec <- missSBM:::SBM_fit_nocovariate$new(A, cl_spec)
  mySBM_fit_kmns <- missSBM:::SBM_fit_nocovariate$new(A, cl_kmns)

  out_hier <- mySBM_fit_hier$doVEM(trace = FALSE, threshold = tol, maxIter = 50, fixPointIter = 3)
  out_spec <- mySBM_fit_spec$doVEM(trace = FALSE, threshold = tol, maxIter = 50, fixPointIter = 3)
  out_kmns <- mySBM_fit_kmns$doVEM(trace = FALSE, threshold = tol, maxIter = 50, fixPointIter = 3)

  BM <- blockmodels::BM_bernoulli("SBM_sym", A, verbosity = 0, explore_max = Q, plotting = "", ncores = 1)
  BM$estimate()
  pi_BM <- BM$model_parameters[[Q]]$pi

  ## checking estimation consistency
  expect_lt(sum((mySBM_fit_spec$connectParam - mySBM$connectParam)^2), tol)
  expect_lt(sum((mySBM_fit_hier$connectParam - mySBM$connectParam)^2), tol)
  expect_lt(sum((mySBM_fit_kmns$connectParam - mySBM$connectParam)^2), tol)

  ## checking estimation consistency with block model
  expect_lt(sum((mySBM_fit_spec$connectParam - pi_BM)^2), tol)
  expect_lt(sum((mySBM_fit_hier$connectParam - pi_BM)^2), tol)
  expect_lt(sum((mySBM_fit_kmns$connectParam - pi_BM)^2), tol)

  ## checking consistency of the clustering
  expect_lt(1 - ARI(mySBM_fit_hier$memberships, mySBM$memberships), tol)
  expect_lt(1 - ARI(mySBM_fit_spec$memberships, mySBM$memberships), tol)
  expect_lt(1 - ARI(mySBM_fit_kmns$memberships, mySBM$memberships), tol)

  expect_equal(mySBM_fit_hier$loglik,
               mySBM_fit_spec$loglik,
               mySBM_fit_kmns$loglik)
  expect_gt(mySBM_fit_hier$loglik - .5 * mySBM_fit_hier$penalty,  BM$ICL[[Q]])

})


test_that("Consistency of VEM of a SBM_fit_nocovariate on a series of values for nbBlocks", {

  BM <- blockmodels::BM_bernoulli("SBM_sym", A, verbosity = 0, explore_min = 5, explore_max = 5, plotting = "", ncores = 1)
  BM$estimate()

  vBlocks <- 1:5
  models <- lapply(vBlocks, function(nbBlocks) {
    cl0 <- missSBM:::init_clustering(mySBM$adjacencyMatrix, nbBlocks, NULL, "hierarchical")
    myFit <- missSBM:::SBM_fit_nocovariate$new(mySBM$adjacencyMatrix, cl0)
    myFit$doVEM()
    myFit
  })

  ICLs  <- sapply(models, function(model) model$ICL)
  bestICL <- models[[which.min(ICLs)]]

  expect_equal(which.min(ICLs), which.max(BM$ICL))

  expect_lt(sum((-.5 * ICLs - BM$ICL)/BM$ICL^2), 5e-6)
})

