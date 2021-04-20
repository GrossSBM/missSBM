set.seed(1234)
library(sbm)
library(aricode)

## Common parameters
N_nocov  <- 200
Q <- 3
source("utils_test.R", local = TRUE)

test_that("SimpleSBM_fit 'Bernoulli' model, undirected, no covariate", {

  sampler_undirected_nocov$rNetwork(store = TRUE)

  ## Construction----------------------------------------------------------------
  mySBM_sbm     <- sbm::SimpleSBM_fit$new(sampler_undirected_nocov$networkData, 'bernoulli', FALSE)
  mySBM_sbm$optimize(estimOptions=list(verbosity = 0, plot = FALSE))
  mySBM_sbm$setModel(3)

  net <- missSBM:::partlyObservedNetwork$new(sampler_undirected_nocov$networkData)
  cl <- net$clustering(3)[[1]]

  mySBM_missSBM <- missSBM:::SimpleSBM_fit_noCov$new(net, clusterInit = cl)
  mySBM_missSBM$doVEM()

  ## correctness

  ## distance with blockmodels/sbm estiamtor
  expect_lt(rmse(mySBM_missSBM$connectParam$mean, mySBM_sbm$connectParam$mean), 0.05)
  expect_gt(ARI(mySBM_missSBM$memberships, mySBM_sbm$memberships), 0.95)
  expect_lt(rmse(mySBM_missSBM$loglik, mySBM_sbm$loglik), 0.01)

  ## distance to true values
  expect_lt(rmse(mySBM_missSBM$connectParam$mean, sampler_undirected_nocov$connectParam$mean), 0.05)
  expect_gt(ARI(mySBM_missSBM$memberships, sampler_undirected_nocov$memberships), 0.95)

})

test_that("SimpleSBM_fit 'Bernoulli' model, directed, no covariate", {

  sampler_directed_nocov$rNetwork(store = TRUE)

  ## Construction----------------------------------------------------------------
  mySBM_sbm     <- sbm::SimpleSBM_fit$new(sampler_directed_nocov$networkData, 'bernoulli', TRUE)
  mySBM_sbm$optimize(estimOptions=list(verbosity = 0))
  mySBM_sbm$setModel(3)

  net <- missSBM:::partlyObservedNetwork$new(sampler_directed_nocov$networkData)
  cl <- net$clustering(3)[[1]]

  mySBM_missSBM <- missSBM:::SimpleSBM_fit_noCov$new(sampler_directed_nocov$networkData, clusterInit = cl)
  mySBM_missSBM$doVEM()

  ## correctness
  ## distance with blockmodels/sbm estiamtor
  expect_lt(rmse(mySBM_missSBM$connectParam$mean, mySBM_sbm$connectParam$mean), 0.1)
  expect_gt(ARI(mySBM_missSBM$memberships, mySBM_sbm$memberships), 0.95)
  expect_lt(rmse(mySBM_missSBM$loglik, mySBM_sbm$loglik), 0.01)

  ## distance to true values
  expect_lt(rmse(mySBM_missSBM$connectParam$mean, sampler_directed_nocov$connectParam$mean), 0.1)
  expect_lt(rmse(mySBM_missSBM$covarParam, sampler_directed_cov$covarParam), 0.2)
  expect_gt(ARI(mySBM_missSBM$memberships, sampler_directed_nocov$memberships), 0.95)

})

