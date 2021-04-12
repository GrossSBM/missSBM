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

test_that("SimpleSBM_fit 'Bernoulli' model, undirected, one covariate, dyad sampling, MAR", {

  ## SIMPLE UNDIRECTED BERNOULLI SBM
  means <- diag(.4, 2) + 0.05
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  # mySampler <- sbm::SimpleSBM$new('bernoulli', nbNodes, FALSE, blockProp, connectParam, covarParam = covarParam[1], covarList = covarList[1])
  mySampler <- sbm::SimpleSBM$new('bernoulli', nbNodes, FALSE, blockProp, connectParam, covarParam = covarParam[1], covarList = covarList[1])
  mySampler$rNetwork(store = TRUE)

  adjMatrix  <- missSBM::observeNetwork(mySampler$networkData, "dyad", 0.5)
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix)
  cl <- net$clustering(2)[[1]]

  mySBM <- missSBM:::SimpleSBM_fit_withCov$new(adjMatrix, cl, covarList[1])
  mySBM$doVEM()
  mySBM$reorder()

  mySBM_MAR <- missSBM:::SimpleSBM_fit_MAR_withCov$new(adjMatrix, cl, covarList[1])
  mySBM_MAR$doVEM()
  mySBM_MAR$reorder()

  ## correctness
  expect_lt(rmse(mySBM_MAR$connectParam$mean, mySampler$connectParam$mean), 0.075)
  expect_lt(rmse(mySBM$connectParam$mean, mySampler$connectParam$mean), 0.075)
  expect_lt(rmse(mySBM$connectParam$mean, mySBM_MAR$connectParam$mean), 1e-5)
  expect_gt(ARI(mySBM_MAR$memberships, mySampler$memberships), 0.95)
  expect_lt(mySBM_MAR$loglik, mySBM$loglik)

})

test_that("SimpleSBM_fit 'Bernoulli' model, directed, one covariate, dyad sampling, MAR", {

  ## SIMPLE UNDIRECTED BERNOULLI SBM
  means <- matrix(rev(c(0.1, 0.4, 0.6, 0.9)), 2,  2)
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM$new('bernoulli', nbNodes, TRUE, blockProp, connectParam, c(node = "nodeName"), covarParam[1], covarList_directed[1])
  mySampler$rNetwork(store = TRUE)

  adjMatrix  <- missSBM::observeNetwork(mySampler$networkData, "dyad", 0.5, covariates = covarList[1])
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix)
  cl <- net$clustering(3)[[1]]

  mySBM <- missSBM:::SimpleSBM_fit_withCov$new(adjMatrix, cl, covarList[1])
  mySBM$doVEM()
  mySBM$reorder()

  mySBM_MAR <- missSBM:::SimpleSBM_fit_MAR_withCov$new(adjMatrix, cl, covarList[1])
  mySBM_MAR$doVEM()
  mySBM_MAR$reorder()

  ## correctness
  expect_lt(rmse(mySBM_MAR$connectParam$mean, mySampler$connectParam$mean), 0.075)
  expect_lt(rmse(mySBM$connectParam$mean, mySampler$connectParam$mean), 0.075)
  expect_lt(rmse(mySBM$connectParam$mean, mySBM_MAR$connectParam$mean), 1e-5)
  expect_gt(ARI(mySBM_MAR$memberships, mySampler$memberships), 0.95)
  expect_lt(mySBM_MAR$loglik, mySBM$loglik)

})
