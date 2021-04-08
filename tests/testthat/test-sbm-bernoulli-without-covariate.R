set.seed(1234)
library(sbm)
library(aricode)

rmse <- function(theta, theta_star) { sqrt(sum((theta - theta_star)^2)/sum(theta_star^2)) }

## Common parameters
nbNodes  <- 200
nbBlocks <- 3
blockProp <- c(.5, .25, .25) # group proportions

test_that("SimpleSBM_fit 'Bernoulli' model, undirected, no covariate", {

  ## SIMPLE UNDIRECTED BERNOULLI SBM
  means <- diag(.4, 3) + 0.05
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- sbm::SimpleSBM$new('bernoulli', nbNodes, FALSE, blockProp, connectParam)
  mySampler$rMemberships(store = TRUE)
  mySampler$rEdges(store = TRUE)

  ## Construction----------------------------------------------------------------
  mySBM_sbm     <- sbm::SimpleSBM_fit$new(mySampler$networkData, 'bernoulli', FALSE)
  mySBM_sbm$optimize(estimOptions=list(verbosity = 0, plot = FALSE))
  mySBM_sbm$setModel(3)

  net <- missSBM:::partlyObservedNetwork$new(mySampler$networkData)
  cl <- net$clustering(3)[[1]]

  mySBM_missSBM <- missSBM:::SimpleSBM_fit_noCov$new(mySampler$networkData, clusterInit = cl)
  mySBM_missSBM$doVEM()

  ## correctness

  ## distance with blockmodels/sbm estiamtor
  expect_lt(rmse(mySBM_missSBM$connectParam$mean, mySBM_sbm$connectParam$mean), 0.05)
  expect_gt(ARI(mySBM_missSBM$memberships, mySBM_sbm$memberships), 0.95)
  expect_lt(rmse(mySBM_missSBM$loglik, mySBM_sbm$loglik), 0.01)

  ## distance to true values
  expect_lt(rmse(mySBM_missSBM$connectParam$mean, mySampler$connectParam$mean), 0.1)
  expect_gt(ARI(mySBM_missSBM$memberships, mySampler$memberships), 0.95)

})

test_that("SimpleSBM_fit 'Bernoulli' model, directed, no covariate", {

  ## SIMPLE UNDIRECTED BERNOULLI SBM
  means <- matrix(c(9:1)/10, 3,  3)
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- sbm::SimpleSBM$new('bernoulli', nbNodes, TRUE, blockProp, connectParam)
  mySampler$rMemberships(store = TRUE)
  mySampler$rEdges(store = TRUE)

  ## Construction----------------------------------------------------------------
  mySBM_sbm     <- sbm::SimpleSBM_fit$new(mySampler$networkData, 'bernoulli', TRUE)
  mySBM_sbm$optimize(estimOptions=list(verbosity = 0))
  mySBM_sbm$setModel(3)

  net <- missSBM:::partlyObservedNetwork$new(mySampler$networkData)
  cl <- net$clustering(3)[[1]]

  mySBM_missSBM <- missSBM:::SimpleSBM_fit_noCov$new(mySampler$networkData, clusterInit = cl)
  mySBM_missSBM$doVEM()
  mySBM_missSBM$reorder()

  ## correctness
  ## distance with blockmodels/sbm estiamtor
  expect_lt(rmse(mySBM_missSBM$connectParam$mean, mySBM_sbm$connectParam$mean), 0.1)
  expect_gt(ARI(mySBM_missSBM$memberships, mySBM_sbm$memberships), 0.95)
  expect_lt(rmse(mySBM_missSBM$loglik, mySBM_sbm$loglik), 0.01)

  ## distance to true values
  expect_lt(rmse(mySBM_missSBM$connectParam$mean, mySampler$connectParam$mean), 0.1)
  expect_gt(ARI(mySBM_missSBM$memberships, mySampler$memberships), 0.95)

})
