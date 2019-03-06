context("test-sbm-fit-without-covariate")

library(aricode)
library(blockmodels)

set.seed(178303)
### A SBM model : ###
N <- 400
Q <- 5
alpha <- rep(1,Q)/Q                     # mixture parameter
pi <- diag(.45,Q) + .05                 # connectivity matrix
directed <- FALSE

### Draw a undirected SBM model
mySBM <- simulateSBM(N, alpha, pi, directed)
A <- mySBM$adjMatrix
cl_rand <- sample(mySBM$memberships)
cl_spec <- missSBM:::init_clustering(A, Q, NULL, "spectral")
cl_hier <- missSBM:::init_clustering(A, Q, NULL,"hierarchical")
cl_kmns <- missSBM:::init_clustering(A, Q, NULL,"kmeans")

test_that("Creation of a SBM_fit_nocovariate", {

  mySBM_fit <- missSBM:::SBM_fit_nocovariate$new(A, cl_rand)
  expect_is(mySBM_fit, c("SBM_fit_nocovariate", "SBM_fit", "SBM", "R6"))
  expect_equal(mySBM_fit$memberships, cl_rand)
  expect_equal(mySBM_fit$df_connectParams, Q * (Q + 1)/2)
  expect_equal(mySBM_fit$df_covarParams, 0)
  expect_equal(mySBM_fit$df_mixtureParams, Q - 1)
  expect_equal(mySBM_fit$blocks, missSBM:::clustering_indicator(cl_rand))
  expect_equal(dim(mySBM_fit$connectParam), dim(mySBM$connectParam))
  expect_equal(length(mySBM_fit$mixtureParam), length(mySBM$mixtureParam))
  expect_equal(mySBM_fit$direction, "undirected")

})

test_that("Consistency of VEM of a SBM_fit_nocovariate when the number of block is given", {

  tol <- 1e-3

  ## testing all initialization
  mySBM_fit_hier <- missSBM:::SBM_fit_nocovariate$new(A, cl_hier)
  mySBM_fit_spec <- missSBM:::SBM_fit_nocovariate$new(A, cl_spec)
  mySBM_fit_kmns <- missSBM:::SBM_fit_nocovariate$new(A, cl_kmns)

  out_hier <- mySBM_fit_hier$doVEM(A, trace = FALSE, threshold = tol, maxIter = 50, fixPointIter = 3)
  out_spec <- mySBM_fit_spec$doVEM(A, trace = FALSE, threshold = tol, maxIter = 50, fixPointIter = 3)
  out_kmns <- mySBM_fit_kmns$doVEM(A, trace = FALSE, threshold = tol, maxIter = 50, fixPointIter = 3)

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

  expect_equal(mySBM_fit_hier$vBound(A),
               mySBM_fit_spec$vBound(A),
               mySBM_fit_kmns$vBound(A))
  expect_gt(mySBM_fit_hier$vBound(A) - .5 * mySBM_fit_hier$penalty,  BM$ICL[[Q]])

})


test_that("Consistency of VEM of a SBM_fit_nocovariate on a series of values for nBlocks", {

  BM <- blockmodels::BM_bernoulli("SBM_sym", A, verbosity = 0, explore_min = 8, explore_max = 8, plotting = "", ncores = 1)
  BM$estimate()

  vBlocks <- 1:8
  models <- lapply(vBlocks, function(nBlocks) {
    cl0 <- missSBM:::init_clustering(mySBM$adjMatrix, nBlocks, NULL, "hierarchical")
    myFit <- missSBM:::SBM_fit_nocovariate$new(mySBM$adjMatrix, cl0)
    myFit$doVEM(mySBM$adjMatrix)
    myFit
  })

  vICLs  <- sapply(models, function(model) model$vICL(mySBM$adjMatrix))
  bestICL <- models[[which.min(vICLs)]]

  expect_equal(which.min(vICLs), which.max(BM$ICL))

  expect_lt(sum((-.5 * vICLs - BM$ICL)/BM$ICL^2), 1e-6)
})

