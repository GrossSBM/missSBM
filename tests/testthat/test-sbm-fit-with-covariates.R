context("test-sbm-fit-with-covariates")

library(aricode)
library(blockmodels)

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
mySBM <- simulateSBM(N, alpha, pi, directed, covarMatrix, covarParam, missSBM:::l1_similarity)
A <- mySBM$adjMatrix

## Formatting covariates for blockmodels
covariates_BM <- lapply(seq(dim(mySBM$covarArray)[3]), function(x) mySBM$covarArray[ , , x])

### Draw a undirected SBM model
cl_rand <- sample(mySBM$memberships)
cl_spec <- missSBM:::init_clustering(A, Q, mySBM$covarArray, "spectral")
cl_hier <- missSBM:::init_clustering(A, Q, mySBM$covarArray, "hierarchical")
cl_kmns <- missSBM:::init_clustering(A, Q, mySBM$covarArray, "kmeans")

test_that("Creation of a SBM_fit_covariates", {

  mySBM_fit <- missSBM:::SBM_fit_covariates$new(A, cl_rand, covarMatrix, missSBM:::l1_similarity)
  expect_is(mySBM_fit, "SBM_fit_covariates")
  expect_equal(mySBM_fit$memberships, cl_rand)
  expect_equal(mySBM_fit$df_connectParams, Q * (Q + 1)/2)
  expect_true(mySBM_fit$hasCovariates)
  expect_equal(mySBM_fit$df_covarParams, M)
  expect_equal(mySBM_fit$df_mixtureParams, Q - 1)
  expect_equal(mySBM_fit$blocks, missSBM:::clustering_indicator(cl_rand))
  expect_equal(dim(mySBM_fit$connectParam), dim(mySBM$connectParam))
  expect_equal(length(mySBM_fit$mixtureParam), length(mySBM$mixtureParam))
  expect_equal(mySBM_fit$direction, "undirected")

})

test_that("Consistency of VEM of a SBM_fit with the number of block given", {

  tol <- 1e-2

  ## testing all initialization
  mySBM_fit_hier <- missSBM:::SBM_fit_covariates$new(A, cl_hier, covarMatrix, missSBM:::l1_similarity)
  mySBM_fit_spec <- missSBM:::SBM_fit_covariates$new(A, cl_spec, covarMatrix, missSBM:::l1_similarity)
  mySBM_fit_kmns <- missSBM:::SBM_fit_covariates$new(A, cl_kmns, covarMatrix, missSBM:::l1_similarity)

  out_hier <- mySBM_fit_hier$doVEM(A, trace = FALSE, threshold = tol, maxIter = 50, fixPointIter = 3)
  out_spec <- mySBM_fit_spec$doVEM(A, trace = FALSE, threshold = tol, maxIter = 50, fixPointIter = 3)
  out_kmns <- mySBM_fit_kmns$doVEM(A, trace = FALSE, threshold = tol, maxIter = 50, fixPointIter = 3)

  ## checking estimation consistency
  expect_lt(sum((mySBM_fit_spec$connectParam - mySBM$connectParam)^2), tol)
  expect_lt(sum((mySBM_fit_hier$connectParam - mySBM$connectParam)^2), tol)
  expect_lt(sum((mySBM_fit_kmns$connectParam - mySBM$connectParam)^2), tol)

  ## checking consistency of the clustering
  expect_lt(1 - ARI(mySBM_fit_hier$memberships, mySBM$memberships), tol)
  expect_lt(1 - ARI(mySBM_fit_spec$memberships, mySBM$memberships), tol)
  expect_lt(1 - ARI(mySBM_fit_kmns$memberships, mySBM$memberships), tol)

  ## checking consistency between the ICL and with the value comptued by blockmodels
  expect_lt(abs(mySBM_fit_hier$vBound(A) - mySBM_fit_spec$vBound(A))/abs(mySBM_fit_spec$vBound(A)), tol)
  expect_lt(abs(mySBM_fit_hier$vBound(A) - mySBM_fit_kmns$vBound(A))/abs(mySBM_fit_kmns$vBound(A)), tol)
  expect_lt(abs(mySBM_fit_kmns$vBound(A) - mySBM_fit_spec$vBound(A))/abs(mySBM_fit_spec$vBound(A)), tol)

})


# test_that("Consistency of VEM of a SBM_fit_nocovariate on a series of values for nBlocks", {
#
#   BM <- blockmodels::BM_bernoulli_covariates("SBM_sym", A, covariates_BM, verbosity = 0, explore_min = 5, explore_max = 5, plotting = "", ncores = 1)
#   BM$estimate()
#
#   vBlocks <- 1:5
#   models <- lapply(vBlocks, function(nBlocks) {
#     cl0 <- missSBM:::init_clustering(mySBM$adjMatrix, nBlocks, mySBM$co, "hierarchical")
#     myFit <- missSBM:::SBM_fit_nocovariate$new(mySBM$adjMatrix, cl0)
#     myFit$doVEM(mySBM$adjMatrix)
#     myFit
#   })
#
#   vICLs  <- sapply(models, function(model) model$vICL(mySBM$adjMatrix))
#   bestICL <- models[[which.min(vICLs)]]
#
#   expect_equal(which.min(vICLs), which.max(BM$ICL))
#
#   expect_lt(sum((-.5 * vICLs - BM$ICL)/BM$ICL^2), 1e-6)
# })
