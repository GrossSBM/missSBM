set.seed(1234)
library(sbm)
library(aricode)

N_cov <- 80
Q <- 2
M <- 1
source("utils_test.R", local = TRUE)

test_that("SimpleSBM_fit 'Bernoulli' model, undirected, one covariate", {

  sampler_undirected_cov$rNetwork(store = TRUE)

  ## blockmodels
  mySBM_sbm <- sbm::SimpleSBM_fit$new(sampler_undirected_cov$networkData, 'bernoulli', FALSE, covarList = covarList_undirected)
  mySBM_sbm$optimize(estimOptions=list(verbosity = 0, plot = FALSE))
  mySBM_sbm$setModel(Q)

  ## missSBM
  net <- missSBM:::partlyObservedNetwork$new(sampler_undirected_cov$networkData, covariates = covarList_undirected)
  cls <- net$clustering(1:(2*Q))
  cl <- cls[[Q]]
  mySBM_missSBM <- missSBM:::SimpleSBM_fit_withCov$new(net, clusterInit = cl, covarList = covarList_undirected)
  mySBM_missSBM$doVEM()

  ## correctness

  ## distance with blockmodels/sbm estiamtor
  expect_lt(rmse(mySBM_missSBM$connectParam$mean, mySBM_sbm$connectParam$mean), 0.05)
  expect_gt(ARI(mySBM_missSBM$memberships, mySBM_sbm$memberships), 0.8)
  expect_lt(rmse(mySBM_missSBM$loglik, mySBM_sbm$loglik), 0.01)

  ## distance to true values
  expect_lt(rmse(mySBM_missSBM$connectParam$mean, sampler_undirected_cov$connectParam$mean), 0.1)
  expect_lt(rmse(mySBM_missSBM$covarParam, sampler_undirected_cov$covarParam), 0.1)
  expect_gt(ARI(mySBM_missSBM$memberships, sampler_undirected_cov$memberships), 0.85)

})

test_that("SimpleSBM_fit 'Bernoulli' model, directed, one covariate", {

  sampler_directed_cov$rNetwork(store = TRUE)

  ## Construction----------------------------------------------------------------
  mySBM_sbm <- sbm::SimpleSBM_fit$new(sampler_directed_cov$networkData, 'bernoulli', TRUE, covarList = covarList_directed)
  mySBM_sbm$optimize(estimOptions=list(verbosity = 0, plot = FALSE))
  mySBM_sbm$setModel(Q)

  ## missSBM
  net <- missSBM:::partlyObservedNetwork$new(sampler_directed_cov$networkData, covariates = covarList_directed)
  cls <- net$clustering(1:(2*Q))
  cl <- cls[[Q]]
  mySBM_missSBM <- missSBM:::SimpleSBM_fit_withCov$new(net, clusterInit = cl, covarList = covarList_directed)
  mySBM_missSBM$doVEM(trace = TRUE)

  ## correctness
  ## distance with blockmodels/sbm estiamtor
  expect_lt(rmse(mySBM_missSBM$connectParam$mean, mySBM_sbm$connectParam$mean), 0.1)
  expect_gt(ARI(mySBM_missSBM$memberships, mySBM_sbm$memberships), 0.8)
  expect_lt(rmse(mySBM_missSBM$loglik, mySBM_sbm$loglik), 0.05)

  ## distance to true values
  expect_lt(rmse(mySBM_missSBM$connectParam$mean, sampler_directed_cov$connectParam$mean), 0.1)
  expect_gt(ARI(mySBM_missSBM$memberships, sampler_directed_cov$memberships), 0.85)

})
