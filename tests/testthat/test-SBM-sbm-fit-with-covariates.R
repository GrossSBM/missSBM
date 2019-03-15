context("test sbm fit with covariates (class SBM_fit_covariates)")

library(aricode)
library(blockmodels)
source("utils_test.R")
## ========================================================================
## A SBM model with covariates

set.seed(178303)
N <- 300
Q <- 3
alpha <- rep(1,Q)/Q                     # mixture parameter
pi <- diag(.45,Q) + .05                 # connectivity matrix
gamma <- missSBM:::logit(pi)
directed <- FALSE

### Draw a SBM model (Bernoulli, undirected) with covariates
M <- 2
covarMatrix <- matrix(rnorm(N*M,mean = 0, sd = 1), N, M)
covarParam  <- rnorm(M, -1, 1)

sbm <- simulateSBM(N, alpha, gamma, directed, covarMatrix, covarParam)

### Draw a undirected SBM model
cl_rand <- sample(sbm$memberships)
cl_spec <- missSBM:::init_clustering(sbm$adjMatrix, Q, sbm$covarArray, "spectral")
cl_hier <- missSBM:::init_clustering(sbm$adjMatrix, Q, sbm$covarArray, "hierarchical")
cl_kmns <- missSBM:::init_clustering(sbm$adjMatrix, Q, sbm$covarArray, "kmeans")

test_that("Creation of a SBM_fit_covariates", {

  mySBM_fit <- missSBM:::SBM_fit_covariates$new(sbm$adjMatrix, cl_rand, sbm$covarArray)
  expect_is(mySBM_fit, "SBM_fit_covariates")
  expect_equal(mySBM_fit$memberships, cl_rand)
  expect_equal(mySBM_fit$df_connectParams, Q * (Q + 1)/2)
  expect_true(mySBM_fit$hasCovariates)
  expect_equal(mySBM_fit$df_covarParams, M)
  expect_equal(mySBM_fit$df_mixtureParams, Q - 1)
  expect_equal(mySBM_fit$blocks, missSBM:::clustering_indicator(cl_rand))
  expect_equal(dim(mySBM_fit$connectParam), dim(sbm$connectParam))
  expect_equal(length(mySBM_fit$mixtureParam), length(sbm$mixtureParam))
  expect_equal(mySBM_fit$direction, "undirected")

})

test_that("Consistency of VEM of a SBM_fit_covariates with the number of block given", {

  tol <- 1e-3

  ## testing just hierarchical clustering (best init)
  mySBM_fit <- missSBM:::SBM_fit_covariates$new(sbm$adjMatrix, cl_spec, sbm$covarArray)

  out <- mySBM_fit$doVEM(sbm$adjMatrix, trace = FALSE, threshold = tol, maxIter = 10, fixPointIter = 3)

  ## A bit long, but it works: we do just as good as blockmodels, sometimes better
  covariates_BM <- lapply(seq(dim(sbm$covarArray)[3]), function(x) sbm$covarArray[ , , x])
  BM <- blockmodels::BM_bernoulli_covariates("SBM_sym", sbm$adjMatrix, covariates_BM, verbosity = 0, explore_min = 3, explore_max = 3, plotting = "", ncores = 1)
  BM$estimate()

  ## similar estimation thant BM for connection parameters
  error_BM      <- error(sbm$connectParam, BM$model_parameters[[3]]$m)
  error_missSBM <- error(sbm$connectParam, mySBM_fit$connectParam)
  expect_lt(abs(error_missSBM - error_BM), 1e-2)

  ## similar estimation thant BM for regression parameters
  error_BM      <- error(sbm$covarParam, as.numeric(BM$model_parameters[[3]]$beta))
  error_missSBM <- error(sbm$covarParam, mySBM_fit$covarParam)
  expect_lt(abs(error_missSBM - error_BM), 1e-2)

  ## checking estimation consistency
  expect_lt(error(logistic(mySBM_fit$connectParam), pi), tol)

  expect_lt(error(mySBM_fit$covarParam, covarParam), tol*10)

  ## checking consistency of the clustering
  expect_lt(1 - ARI(mySBM_fit$memberships, sbm$memberships), tol)

  tol <- 1e-3

})

# ## CONSISTENCY WITH BLOCKMODELS AND ON TIMOTHÉE'S EXAMPLE
#
# referenceResults <- readRDS(system.file("extdata", "referenceResults.rds", package = "missSBM"))
#
# test_that("Consistency of VEM of a SBM_fit_covariates on a series of values for nBlocks", {
#
#   truth   <- referenceResults$true_sbm_cov
#   cl_init <- missSBM:::init_clustering(truth$adjMatrix, Q, truth$covarArray, "spectral")
#
#   ## testing just hierarchical clustering (best init)
#   mySBM_fit <- missSBM:::SBM_fit_covariates$new(truth$adjMatrix, cl_init, truth$covarArray)
#
#   out <- mySBM_fit$doVEM(truth$adjMatrix, trace = FALSE, threshold = 1e-4, maxIter = 10, fixPointIter = 3)
#
#   ## A bit long, but it works: we do just as good as blockmodels, sometimes better
#   covariates_BM <- lapply(seq(dim(truth$covarArray)[3]), function(x) truth$covarArray[ , , x])
#   BM <- blockmodels::BM_bernoulli_covariates("SBM_sym", truth$adjMatrix, covariates_BM, verbosity = 0, explore_min = 3, explore_max = 3, plotting = "", ncores = 1)
#   BM$estimate()
#
#   ## similar estimation thant BM for connection parameters
#   error_BM      <- error(logistic(sbm$connectParam), logistic(BM$model_parameters[[3]]$m))
#   error_missSBM <- error(logistic(sbm$connectParam), logistic(mySBM_fit$connectParam))
#
#   if (error_missSBM > error_BM) {
#     expect_lt(abs(error_missSBM - error_BM), 5e-2)
#   } else {
#     cat('bette than BM.')
#   }
#
#   ## similar estimation thant BM for regression parameters
#   error_BM      <- error(truth$covarParam, as.numeric(BM$model_parameters[[3]]$beta))
#   error_missSBM <- error(truth$covarParam, mySBM_fit$covarParam)
#   if (error_missSBM > error_BM) {
#     expect_lt(abs(error_missSBM - error_BM), 5e-2)
#   } else {
#     cat('bette than BM.')
#   }
#
#   ## checking consistency of the clustering
#   expect_gt(ARI(mySBM_fit$memberships, truth$memberships), 0.6)
#
# })


## CONSISTENCY WITH BLOCKMODELS

test_that("Consistency of VEM of a SBM_fit_covariates on a series of values for nBlocks", {

  ## ========================================================================
  ## A SBM model with covariates
  set.seed(178304)
  N <- 40
  Q <- 2
  alpha <- rep(1,Q)/Q                     # mixture parameter
  pi <- diag(.45,Q) + .05                 # connectivity matrix
  directed <- FALSE
  gamma <- missSBM:::logit(pi)

  ### Draw a SBM model (Bernoulli, undirected) with covariates
  M <- 2
  covarMatrix <- matrix(rnorm(N*M,mean = 0, sd = 1), N, M)
  covarParam  <- rnorm(M,0,1)
  sbm <- simulateSBM(N, alpha, gamma, directed, covarMatrix, covarParam, missSBM:::l1_similarity)

  ## Formatting covariates for blockmodels
  covariates_BM <- lapply(seq(dim(sbm$covarArray)[3]), function(x) sbm$covarArray[ , , x])

  BM <- blockmodels::BM_bernoulli_covariates("SBM_sym", sbm$adjMatrix, covariates_BM, verbosity = 0, explore_min = 4, explore_max = 4, plotting = "", ncores = 1)
  BM$estimate()

  vBlocks <- 1:4
  models <- lapply(vBlocks, function(nBlocks) {
    cl0 <- missSBM:::init_clustering(sbm$adjMatrix, nBlocks, sbm$covarArray, "hierarchical")
    myFit <- missSBM:::SBM_fit_covariates$new(sbm$adjMatrix, cl0, sbm$covarArray)
    myFit$doVEM(sbm$adjMatrix)
    myFit
  })

  vICLs  <- sapply(models, function(model) model$vICL(sbm$adjMatrix))
  bestICL <- models[[which.min(vICLs)]]

  expect_equal(which.min(vICLs), which.max(BM$ICL))

  tol_ref   <- 1e-2
  tol_truth <- 1e-2
  expect_lt(sum(((-.5 * vICLs - BM$ICL)/BM$ICL)^2), tol_ref)

  error_missSBM <- error(logistic(sbm$connectParam), logistic(bestICL$connectParam))
  error_BM      <- error(logistic(bestICL$connectParam),
                         logistic(BM$model_parameters[[2]]$m))
  error_BM_beta <- error(bestICL$covarParam, BM$model_parameters[[2]]$beta)

  expect_lt(error_missSBM, tol_truth)
  ## we are close to Blockmodels...
  expect_lt(error_BM     , tol_ref)
  expect_lt(error_BM_beta, tol_ref)

})
