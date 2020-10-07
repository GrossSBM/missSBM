context("test sbm fit without covariate (class SimpleSBM_fit_missSBM)")

library(aricode)
library(blockmodels)

set.seed(178303)
### A SBM model : ###
N <- 100
Q <- 3
pi <- rep(1, Q)/Q                     # mixture parameter
theta <- list(mean = diag(.45, Q, Q) + .05)                 # connectivity matrix
directed <- FALSE

### Draw a undirected SBM model
mySBM <- sbm::sampleSimpleSBM(N, pi, theta)
A <- mySBM$netMatrix
### UGLY FIX - FIXME
diag(A) <- 0
cl_rand <- base::sample(mySBM$memberships)
cl_spec <- missSBM:::init_clustering(A, Q, NULL, "spectral")
cl_hier <- missSBM:::init_clustering(A, Q, NULL, "hierarchical")
cl_kmns <- missSBM:::init_clustering(A, Q, NULL, "kmeans")

test_that("Creation of a SimpleSBM_fit_missSBM", {

  mySBM_fit <- missSBM:::SimpleSBM_fit_missSBM$new(A, cl_rand)
  expect_is(mySBM_fit, c("SimpleSBM_fit_missSBM", "SimpleSBM_fit", "SBM", "R6"))
  expect_equal(mySBM_fit$memberships, cl_rand)
  expect_equal(mySBM_fit$nbConnectParam, Q * (Q + 1)/2)
  expect_equal(mySBM_fit$nbCovariates, 0)
  expect_equal(mySBM_fit$probMemberships, missSBM:::clustering_indicator(cl_rand))
  expect_equal(dim(mySBM_fit$connectParam$mean), dim(mySBM$connectParam$mean))
  expect_equal(length(mySBM_fit$blockProp), length(mySBM$blockProp))
  expect_equal(mySBM_fit$directed, FALSE)

})

test_that("Consistency of VEM of a SimpleSBM_fit_missSBM when the number of block is given", {

  tol <- 2e-3

  ## testing all initialization
  mySBM_fit_hier <- missSBM:::SimpleSBM_fit_missSBM$new(A, cl_hier)
  mySBM_fit_spec <- missSBM:::SimpleSBM_fit_missSBM$new(A, cl_spec)
  mySBM_fit_kmns <- missSBM:::SimpleSBM_fit_missSBM$new(A, cl_kmns)

  out_hier <- mySBM_fit_hier$doVEM(trace = FALSE, threshold = tol, maxIter = 50, fixPointIter = 3)
  out_spec <- mySBM_fit_spec$doVEM(trace = FALSE, threshold = tol, maxIter = 50, fixPointIter = 3)
  out_kmns <- mySBM_fit_kmns$doVEM(trace = FALSE, threshold = tol, maxIter = 50, fixPointIter = 3)

  BM <- blockmodels::BM_bernoulli("SBM_sym", A, verbosity = 0, explore_max = Q, plotting = "", ncores = 1)
  BM$estimate()
  theta_BM <- BM$model_parameters[[Q]]$pi

  ## checking estimation consistency
  expect_lt(sum((mySBM_fit_spec$connectParam$mean - mySBM$connectParam$mean)^2), tol)
  expect_lt(sum((mySBM_fit_hier$connectParam$mean - mySBM$connectParam$mean)^2), tol)
  expect_lt(sum((mySBM_fit_kmns$connectParam$mean - mySBM$connectParam$mean)^2), tol)

  ## checking estimation consistency with block model
  expect_lt(sum((mySBM_fit_spec$connectParam$mean - theta_BM)^2), tol)
  expect_lt(sum((mySBM_fit_hier$connectParam$mean - theta_BM)^2), tol)
  expect_lt(sum((mySBM_fit_kmns$connectParam$mean - theta_BM)^2), tol)

  ## checking consistency of the clustering
  expect_lt(1 - ARI(mySBM_fit_hier$memberships, mySBM$memberships), tol)
  expect_lt(1 - ARI(mySBM_fit_spec$memberships, mySBM$memberships), tol)
  expect_lt(1 - ARI(mySBM_fit_kmns$memberships, mySBM$memberships), tol)

  expect_equal(mySBM_fit_hier$loglik,
               mySBM_fit_spec$loglik,
               mySBM_fit_kmns$loglik)
  expect_gt(mySBM_fit_hier$loglik - .5 * mySBM_fit_hier$penalty,  BM$ICL[[Q]])

})


test_that("Consistency of VEM of a SimpleSBM_fit_missSBM on a series of values for nbBlocks", {

  BM <- blockmodels::BM_bernoulli("SBM_sym", A, verbosity = 0, explore_min = 5, explore_max = 5, plotting = "", ncores = 1)
  BM$estimate()

  vBlocks <- 1:5
  models <- lapply(vBlocks, function(nbBlocks) {
    cl0 <- missSBM:::init_clustering(A, nbBlocks, NULL, "hierarchical")
    myFit <- missSBM:::SimpleSBM_fit_missSBM$new(A, cl0)
    myFit$doVEM()
    myFit
  })

  ICLs  <- sapply(models, function(model) model$ICL)
  bestICL <- models[[which.min(ICLs)]]

  expect_equal(which.min(ICLs), which.max(BM$ICL))

  expect_lt(sum((-.5 * ICLs - BM$ICL)/BM$ICL^2), 5e-6)
})

