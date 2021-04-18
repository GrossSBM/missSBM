set.seed(1234)
library(sbm)
library(aricode)

N_cov   <- 100
Q <- 2
M <- 1
source("utils_test.R", local = TRUE)

test_that("SimpleSBM_fit 'Bernoulli' model, undirected, one covariate, dyad sampling, MAR", {

  sampler_undirected_cov$rNetwork(store = TRUE)

  adjMatrix  <- missSBM::observeNetwork(sampler_undirected_cov$networkData, "dyad", 0.8)
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix, covariates = covarList_undirected)
  cls <- net$clustering(1:(2*Q))
  cl <- cls[[Q]]

  mySBM <- missSBM:::SimpleSBM_fit_withCov$new(net, cl, covarList_undirected)
  mySBM$doVEM()

  ## correctness
  expect_lt(rmse(mySBM$connectParam$mean, sampler_undirected_cov$connectParam$mean), 0.075)
  expect_gt(ARI(mySBM$memberships       , sampler_undirected_cov$memberships), 0.95)

})

test_that("SimpleSBM_fit 'Bernoulli' model, directed, one covariate, dyad sampling, MAR", {

  sampler_directed_cov$rNetwork(store = TRUE)

  adjMatrix  <- missSBM::observeNetwork(sampler_directed_cov$networkData, "dyad", 0.8, covariates = covarList_directed)
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix, covariates = covarList_directed)
  cls <- net$clustering(1:(2*Q))
  cl <- cls[[Q]]

  mySBM <- missSBM:::SimpleSBM_fit_withCov$new(adjMatrix, cl, covarList_directed)
  mySBM$doVEM()

  ## correctness
  expect_lt(rmse(mySBM$connectParam$mean, sampler_directed_cov$connectParam$mean), 0.075)
  expect_gt(ARI(mySBM$memberships       , sampler_directed_cov$memberships), 0.95)

})
