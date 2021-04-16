context("test sbm fit with covariates (class simpleSBM_fit_missSBM)")

library(aricode)
library(blockmodels)
error <- function(beta1, beta2, sort = FALSE) {
  if (sort)
    err <- sum((sort(beta1) - sort(beta2))^2)/length(beta2)
  else
    err <- sum((beta1 - beta2)^2)/length(beta2)
  err
}

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
cl_rand <- base::sample.int(Q, N, replace = TRUE)

test_that("Creation of a SimpleSBM_fit_missSBM", {

  mySBM_fit <- missSBM:::SimpleSBM_fit_missSBM$new(sbm$networkData, cl_rand, missSBM:::array2list(sbm$covarArray))
  expect_is(mySBM_fit, "SimpleSBM_fit_missSBM")
  expect_equal(mySBM_fit$memberships, cl_rand)
  expect_equal(mySBM_fit$nbConnectParam, Q * (Q + 1)/2)
  expect_equal(mySBM_fit$nbCovariates, M)
  expect_equal(mySBM_fit$probMemberships, missSBM:::clustering_indicator(cl_rand))
  expect_equal(dim(mySBM_fit$connectParam$mean), dim(sbm$connectParam$mean))
  expect_equal(length(mySBM_fit$blockProp), length(sbm$blockProp))
  expect_equal(mySBM_fit$directed, FALSE)

})

## CONSISTENCY WITH BLOCKMODELS

test_that("Consistency of VEM of a SimpleSBM_fit_missSBM on a series of values for nbBlocks", {

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
  adjMatrix <- sbm$networkData
  diag(adjMatrix) <- 0

  ## Formatting covariates for blockmodels
  BM <- blockmodels::BM_bernoulli_covariates("SBM_sym", adjMatrix, covariates, verbosity = 0, explore_min = 3, explore_max = 3, plotting = "", ncores = 1)
  BM$estimate()

  myNet <- missSBM:::partlyObservedNetwork$new(sbm$networkData, covariates = covariates)
  cl_spec <- myNet$clustering(1:3)

  models <- lapply(cl_spec, function(cl0) {
    myFit <- missSBM:::SimpleSBM_fit_withCov$new(myNet, cl0, covariates)
    myFit$doVEM()
    myFit
  })

  ICLs  <- sapply(models, function(model) model$ICL)
  bestICL <- models[[which.min(ICLs)]]

  expect_equal(which.min(ICLs), which.max(BM$ICL))

  tol_ref   <- 1e-2
  tol_truth <- 1e-3
  expect_lt(sum(((-.5 * ICLs - BM$ICL)/BM$ICL)^2), tol_ref)

  error_missSBM <- error(sbm$connectParam$mean, bestICL$connectParam$mean)
  error_BM      <- error(bestICL$connectParam$mean,
                         missSBM:::.logistic(BM$model_parameters[[2]]$m))
  error_BM_beta <- error(bestICL$covarParam, BM$model_parameters[[2]]$beta)

  expect_lt(error_missSBM, tol_truth)
  ## we are close to Blockmodels...
  expect_lt(error_BM     , tol_ref)
  expect_lt(error_BM_beta, tol_ref)

})
