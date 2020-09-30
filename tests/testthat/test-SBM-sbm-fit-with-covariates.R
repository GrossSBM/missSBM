context("test sbm fit with covariates (class SBM_fit_covariates)")

library(aricode)
library(blockmodels)
source("utils_test.R")
## ========================================================================
## A SBM model with covariates

set.seed(178303)
N <- 100
Q <- 3
alpha <- rep(1,Q)/Q                     # mixture parameter
pi <- diag(.45,Q) + .05                 # connectivity matrix
gamma <- missSBM:::logit(pi)
directed <- FALSE

### Draw a SBM model (Bernoulli, undirected) with covariates
M <- 1
covariates <- replicate(M, matrix(rnorm(N*N,mean = 0, sd = 1), N, N), simplify = FALSE)
covarParam  <- rnorm(M, 0, 1)
sbm <- missSBM::simulate(N, alpha, gamma, directed, covariates, covarParam)

### Draw a undirected SBM model
cl_rand <- base::sample(sbm$memberships)
cl_spec <- missSBM:::init_clustering(sbm$adjacencyMatrix, Q, sbm$covarArray, "spectral")
cl_hier <- missSBM:::init_clustering(sbm$adjacencyMatrix, Q, sbm$covarArray, "hierarchical")
cl_kmns <- missSBM:::init_clustering(sbm$adjacencyMatrix, Q, sbm$covarArray, "kmeans")

test_that("Creation of a SBM_fit_covariates", {

  mySBM_fit <- missSBM:::SBM_fit_covariates$new(sbm$adjacencyMatrix, cl_rand, sbm$covarArray)
  expect_is(mySBM_fit, "SBM_fit_covariates")
  expect_equal(mySBM_fit$memberships, cl_rand)
  expect_equal(mySBM_fit$df_connectParams, Q * (Q + 1)/2)
  expect_true(mySBM_fit$nbCovariates > 0)
  expect_equal(mySBM_fit$df_covarParams, M)
  expect_equal(mySBM_fit$df_blockProps, Q - 1)
  expect_equal(mySBM_fit$probMemberships, missSBM:::clustering_indicator(cl_rand))
  expect_equal(dim(mySBM_fit$connectParam), dim(sbm$connectParam))
  expect_equal(length(mySBM_fit$blockProp), length(sbm$blockProp))
  expect_equal(mySBM_fit$direction, "undirected")

})

# test_that("Consistency of VEM of a SBM_fit_covariates with the number of block given", {
#
#   tol <- 1e-3
#
#   ## testing just hierarchical clustering (best init)
#   mySBM_fit <- missSBM:::SBM_fit_covariates$new(sbm$adjMatrix, cl_spec, sbm$covarArray)
#
#   out <- mySBM_fit$doVEM(sbm$adjMatrix, trace = FALSE, threshold = tol, maxIter = 10, fixPointIter = 3)
#
#   ## A bit long, but it works: we do just as good as blockmodels, sometimes better
#   BM <- blockmodels::BM_bernoulli_covariates("SBM_sym", sbm$adjMatrix, covariates, verbosity = 0, explore_min = 3, explore_max = 3, plotting = "", ncores = 1)
#   BM$estimate()
#
#   ## similar estimation thant BM for connection parameters
#   error_BM      <- error(sbm$connectParam, BM$model_parameters[[3]]$m)
#   error_missSBM <- error(sbm$connectParam, mySBM_fit$connectParam)
#   expect_lt(abs(error_missSBM - error_BM), 1e-2)
#
#   ## similar estimation thant BM for regression parameters
#   error_BM      <- error(sbm$covarParam, as.numeric(BM$model_parameters[[3]]$beta))
#   error_missSBM <- error(sbm$covarParam, mySBM_fit$covarParam)
#   expect_lt(abs(error_missSBM - error_BM), 1e-2)
#
#   ## checking estimation consistency
#   expect_lt(error(logistic(mySBM_fit$connectParam), pi), tol)
#
#   expect_lt(error(mySBM_fit$covarParam, covarParam), tol*10)
#
#   ## checking consistency of the clustering
#   expect_lt(1 - ARI(mySBM_fit$memberships, sbm$memberships), tol)
#
#   tol <- 1e-3
#
# })

## CONSISTENCY WITH BLOCKMODELS

test_that("Consistency of VEM of a SBM_fit_covariates on a series of values for nbBlocks", {

  ## ========================================================================
  ## A SBM model with covariates
  set.seed(178304)
  N <- 40
  Q <- 2
  alpha <- rep(1, Q)/Q                     # mixture parameter
  pi <- diag(.45, Q) + .05                 # connectivity matrix
  directed <- FALSE
  gamma <- missSBM:::logit(pi)

  ### Draw a SBM model (Bernoulli, undirected) with covariates
  M <- 1
  covariates_node <- replicate(M, rnorm(N,mean = 0, sd = 1), simplify = FALSE)
  covarMatrix <- simplify2array(covariates_node)
  covarArray  <- missSBM:::getCovarArray(covarMatrix, missSBM:::l1_similarity)
  covariates_dyad <- lapply(seq(dim(covarArray)[3]), function(x) covarArray[ , , x])
  covarParam  <- rnorm(M, 0, 1)
  sbm <- missSBM::simulate(N, alpha, gamma, directed, covariates_dyad, covarParam)

  ## Formatting covariates for blockmodels
  BM <- blockmodels::BM_bernoulli_covariates("SBM_sym", sbm$adjacencyMatrix, covariates_dyad, verbosity = 0, explore_min = 3, explore_max = 3, plotting = "", ncores = 1)
  BM$estimate()

  vBlocks <- 1:3
  models <- lapply(vBlocks, function(nbBlocks) {
    cl0 <- missSBM:::init_clustering(sbm$adjacencyMatrix, nbBlocks, sbm$covarArray, "hierarchical")
    myFit <- missSBM:::SBM_fit_covariates$new(sbm$adjacencyMatrix, cl0, sbm$covarArray)
    myFit$doVEM()
    myFit
  })

  ICLs  <- sapply(models, function(model) model$ICL)
  bestICL <- models[[which.min(ICLs)]]

  expect_equal(which.min(ICLs), which.max(BM$ICL))

  tol_ref   <- 1e-2
  tol_truth <- 1e-2
  expect_lt(sum(((-.5 * ICLs - BM$ICL)/BM$ICL)^2), tol_ref)

  error_missSBM <- error(logistic(sbm$connectParam), logistic(bestICL$connectParam))
  error_BM      <- error(logistic(bestICL$connectParam),
                         logistic(BM$model_parameters[[2]]$m))
  error_BM_beta <- error(bestICL$covarParam, BM$model_parameters[[2]]$beta)

  expect_lt(error_missSBM, tol_truth)
  ## we are close to Blockmodels...
  expect_lt(error_BM     , tol_ref)
  expect_lt(error_BM_beta, tol_ref)

})
