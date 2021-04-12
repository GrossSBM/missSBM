set.seed(1234)
library(sbm)
library(aricode)

rmse <- function(theta, theta_star) { sqrt(sum((theta - theta_star)^2)/sum(theta_star^2)) }

## Common parameters
nbNodes  <- 200
nbBlocks <- 3
blockProp <- c(.5, .25, .25) # group proportions

test_that("SimpleSBM_fit 'Bernoulli' model, undirected, no covariate, dyad sampling, MAR", {

  ## SIMPLE UNDIRECTED BERNOULLI SBM
  means <- diag(.4, 3) + 0.05
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- sbm::SimpleSBM$new('bernoulli', nbNodes, FALSE, blockProp, connectParam)
  mySampler$rNetwork(store = TRUE)

  adjMatrix  <- missSBM::observeNetwork(mySampler$networkData, "dyad", 0.5)
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix)
  cl <- net$clustering(3)[[1]]

  mySBM <- missSBM:::SimpleSBM_fit_noCov$new(adjMatrix, clusterInit = cl)
  mySBM$doVEM()
  mySBM$reorder()

  ## correctness
  expect_lt(rmse(mySBM$connectParam$mean, mySampler$connectParam$mean), 0.075)
  expect_gt(ARI(mySBM$memberships, mySampler$memberships), 0.95)

})

test_that("SimpleSBM_fit 'Bernoulli' model, undirected, no covariate, node sampling, MAR", {

  ## SIMPLE UNDIRECTED BERNOULLI SBM
  means <- diag(.4, 3) + 0.05
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- sbm::SimpleSBM$new('bernoulli', nbNodes, FALSE, blockProp, connectParam)
  mySampler$rNetwork(store = TRUE)

  adjMatrix  <- missSBM::observeNetwork(mySampler$networkData, "node", 0.5)
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix)
  cl <- net$clustering(3)[[1]]

  mySBM <- missSBM:::SimpleSBM_fit_noCov$new(adjMatrix, clusterInit = cl)
  mySBM$doVEM()
  mySBM$reorder()

  ## correctness
  expect_lt(rmse(mySBM$connectParam$mean, mySampler$connectParam$mean), 0.075)
  expect_gt(ARI(mySBM$memberships, mySampler$memberships), 0.95)

})

test_that("SimpleSBM_fit 'Bernoulli' model, directed, no covariate, dyad sampling, MAR", {

  ## SIMPLE UNDIRECTED BERNOULLI SBM
  means <- matrix(c(9:1)/10, 3,  3)
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- sbm::SimpleSBM$new('bernoulli', nbNodes, TRUE, blockProp, connectParam)
  mySampler$rNetwork(store = TRUE)

  adjMatrix  <- missSBM::observeNetwork(mySampler$networkData, "dyad", 0.5)
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix)
  cl <- net$clustering(3)[[1]]

  mySBM <- missSBM:::SimpleSBM_fit_noCov$new(adjMatrix, clusterInit = cl)
  mySBM$doVEM()
  mySBM$reorder()

  ## correctness
  expect_lt(rmse(mySBM$connectParam$mean, mySampler$connectParam$mean), 0.075)
  expect_gt(ARI(mySBM$memberships, mySampler$memberships), 0.95)

})

test_that("SimpleSBM_fit 'Bernoulli' model, directed, no covariate, node samping, MAR", {

  ## SIMPLE UNDIRECTED BERNOULLI SBM
  means <- matrix(c(9:1)/10, 3,  3)
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- sbm::SimpleSBM$new('bernoulli', nbNodes, TRUE, blockProp, connectParam)
  mySampler$rNetwork(store = TRUE)

  adjMatrix  <- missSBM::observeNetwork(mySampler$networkData, "node", 0.2)
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix)
  cl <- net$clustering(3)[[1]]

  mySBM <- missSBM:::SimpleSBM_fit_noCov$new(adjMatrix, clusterInit = cl)
  mySBM$doVEM()
  mySBM$reorder()

  ## correctness
  expect_lt(rmse(mySBM$connectParam$mean, mySampler$connectParam$mean), 0.075)
  expect_gt(ARI(mySBM$memberships, mySampler$memberships), 0.95)

})
