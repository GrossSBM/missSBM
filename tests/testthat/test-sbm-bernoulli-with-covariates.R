set.seed(1234)
library(sbm)
library(aricode)

rmse <- function(theta, theta_star) { sqrt(sum((theta - theta_star)^2)/sum(theta_star^2)) }

## Common parameters
nbNodes  <- 40
nbBlocks <- 2
blockProp <- c(.5, .5) # group proportions
covarParam <- c(-2,2)
dimLabels <- list(row = "rowLabel", col = "colLabel")
covar1 <- matrix(rnorm(nbNodes**2), nbNodes, nbNodes)
covar2 <- matrix(rnorm(nbNodes**2), nbNodes, nbNodes)
covarList_directed <- list(covar1 = covar1, covar2 = covar2)

covar1 <- covar1 + t(covar1)
covar2 <- covar2 + t(covar2)
covarList <- list(covar1 = covar1, covar2 = covar2)

test_that("SimpleSBM_fit 'Bernoulli' model, undirected, one covariate", {

  ## SIMPLE UNDIRECTED BERNOULLI SBM
  means <- diag(.4, 2) + 0.05
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM$new('bernoulli', nbNodes, FALSE, blockProp, connectParam, covarParam = covarParam[1], covarList = covarList[1])
  mySampler$rMemberships(store = TRUE)
  mySampler$rEdges(store = TRUE)

  ## Construction----------------------------------------------------------------
  mySBM_sbm <- sbm::SimpleSBM_fit$new(mySampler$networkData, 'bernoulli', FALSE, covarList = covarList[1])
  mySBM_sbm$optimize(estimOptions=list(verbosity = 0, plot = FALSE))
  mySBM_sbm$setModel(2)

  cl <- missSBM:::init_spectral(mySampler$networkData, 2)
  mySBM_missSBM <- missSBM:::SimpleSBM_fit_withCov$new(mySampler$networkData, clusterInit = cl, covarList = covarList[1])
  mySBM_missSBM$doVEM()

  ## correctness

  ## distance with blockmodels/sbm estiamtor
  expect_lt(rmse(mySBM_missSBM$connectParam$mean, mySBM_sbm$connectParam$mean), 0.1)
  expect_gt(ARI(mySBM_missSBM$memberships, mySBM_sbm$memberships), 0.95)
  expect_lt(rmse(mySBM_missSBM$loglik, mySBM_sbm$loglik), 0.01)

  ## distance to true values
  expect_lt(rmse(mySBM_missSBM$connectParam$mean, mySampler$connectParam$mean), 0.11)
  expect_gt(ARI(mySBM_missSBM$memberships, mySampler$memberships), 0.95)

})

test_that("SimpleSBM_fit 'Bernoulli' model, directed, one covariate", {

  ## SIMPLE UNDIRECTED BERNOULLI SBM
  means <- matrix(rev(c(0.1, 0.4, 0.6, 0.9)), 2,  2)
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM$new('bernoulli', nbNodes, TRUE, blockProp, connectParam, c(node = "nodeName"), covarParam[1], covarList_directed[1])
  mySampler$rMemberships(store = TRUE)
  mySampler$rEdges(store = TRUE)

  ## Construction----------------------------------------------------------------
  mySBM_sbm <- sbm::SimpleSBM_fit$new(mySampler$networkData, 'bernoulli', TRUE, covarList = covarList_directed[1])
  mySBM_sbm$optimize(estimOptions=list(verbosity = 0))
  mySBM_sbm$setModel(2)

  cl <- missSBM:::init_spectral(mySampler$networkData, 2)
  mySBM_missSBM <- missSBM:::SimpleSBM_fit_withCov$new(mySampler$networkData, clusterInit = cl, covarList = covarList_directed[1])
  mySBM_missSBM$doVEM()
  mySBM_missSBM$reorder()

  ## correctness
  ## distance with blockmodels/sbm estiamtor
  expect_lt(rmse(mySBM_missSBM$connectParam$mean, mySBM_sbm$connectParam$mean), 0.1)
  expect_gt(ARI(mySBM_missSBM$memberships, mySBM_sbm$memberships), 0.8)
  expect_lt(rmse(mySBM_missSBM$loglik, mySBM_sbm$loglik), 0.05)

  ## distance to true values
  expect_lt(rmse(mySBM_missSBM$connectParam$mean, mySampler$connectParam$mean), 0.1)
  expect_gt(ARI(mySBM_missSBM$memberships, mySampler$memberships), 0.95)

})
