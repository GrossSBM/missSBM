context("test sbm fit with covariates (class SBM_fit_covariates)")

library(aricode)
library(blockmodels)
source("utils_test.R")
## ========================================================================
## A SBM model with covariates

set.seed(178303)
N <- 100
Q <- 3
pi <- rep(1,Q)/Q                        # mixture parameter
theta <- list(mean = diag(.45,Q) + .05) # connectivity matrix
directed <- FALSE

### Draw a SBM model (Bernoulli, undirected) with covariates
M <- 1
covariates <- replicate(M, matrix(rnorm(N*N,mean = 0, sd = 1), N, N), simplify = FALSE)
covarParam  <- rnorm(M, 0, 1)
sbm <- sbm::sampleSimpleSBM(N, pi, theta, covariates = covariates, covariatesParam = covarParam)

### Draw a undirected SBM model
cl_rand <- base::sample(sbm$memberships)
cl_spec <- missSBM:::init_clustering(sbm$netMatrix, Q, sbm$covarArray, "spectral")
cl_hier <- missSBM:::init_clustering(sbm$netMatrix, Q, sbm$covarArray, "hierarchical")
cl_kmns <- missSBM:::init_clustering(sbm$netMatrix, Q, sbm$covarArray, "kmeans")

test_that("Creation of a SBM_fit_covariates", {

  mySBM_fit <- missSBM:::SBM_fit_covariates$new(sbm$netMatrix, cl_rand, missSBM:::array2list(sbm$covarArray))
  expect_is(mySBM_fit, "SBM_fit_covariates")
  expect_equal(mySBM_fit$memberships, cl_rand)
  expect_equal(mySBM_fit$nbConnectParam, Q * (Q + 1)/2)
  expect_equal(mySBM_fit$nbCovariates, M)
  expect_equal(mySBM_fit$probMemberships, missSBM:::clustering_indicator(cl_rand))
  expect_equal(dim(mySBM_fit$connectParam$mean), dim(sbm$connectParam$mean))
  expect_equal(length(mySBM_fit$blockProp), length(sbm$blockProp))
  expect_equal(mySBM_fit$directed, FALSE)

})

## CONSISTENCY WITH BLOCKMODELS

test_that("Consistency of VEM of a SBM_fit_covariates on a series of values for nbBlocks", {

  ## ========================================================================
  ## A SBM model with covariates
  set.seed(178314)
  N <- 40
  Q <- 2
  pi <- rep(1, Q)/Q                        # mixture parameter
  theta <- list(mean = diag(.45, Q) + .05) # connectivity matrix
  directed <- FALSE

  ### Draw a SBM model (Bernoulli, undirected) with covariates
  M <- 1
  covariates <- replicate(M, {
      cov <- matrix(rnorm(N*N,mean = 0, sd = 1), N, N)
      cov + t(cov)
    }, simplify = FALSE)

  covarParam  <- rnorm(M, 0, 1)

  covarParam  <- rnorm(M, 0, 1)
  sbm <- sbm::sampleSimpleSBM(N, pi, theta, covariates = covariates, covariatesParam = covarParam)
  adjMatrix <- sbm$netMatrix
  diag(adjMatrix) <- 0

  ## Formatting covariates for blockmodels
  BM <- blockmodels::BM_bernoulli_covariates("SBM_sym", adjMatrix, covariates, verbosity = 0, explore_min = 3, explore_max = 3, plotting = "", ncores = 1)
  BM$estimate()

  vBlocks <- 1:3
  models <- lapply(vBlocks, function(nbBlocks) {
    cl0 <- missSBM:::init_clustering(adjMatrix, nbBlocks, sbm$covarArray, "hierarchical")
    myFit <- missSBM:::SBM_fit_covariates$new(adjMatrix, cl0, covariates)
    myFit$doVEM()
    myFit
  })

  ICLs  <- sapply(models, function(model) model$ICL)
  bestICL <- models[[which.min(ICLs)]]

  expect_equal(which.min(ICLs), which.max(BM$ICL))

  tol_ref   <- 1e-2
  tol_truth <- 1e-2
  expect_lt(sum(((-.5 * ICLs - BM$ICL)/BM$ICL)^2), tol_ref)

  error_missSBM <- error(sbm$connectParam$mean, bestICL$connectParam$mean)
  error_BM      <- error(bestICL$connectParam$mean,
                         .logistic(BM$model_parameters[[2]]$m))
  error_BM_beta <- error(bestICL$covarParam, BM$model_parameters[[2]]$beta)

  expect_lt(error_missSBM, tol_truth)
  ## we are close to Blockmodels...
  expect_lt(error_BM     , tol_ref)
  expect_lt(error_BM_beta, tol_ref)

})
