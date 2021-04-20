set.seed(1234)
library(sbm)
library(aricode)

N_nocov   <- 200
Q <- 3
source("utils_test.R", local = TRUE)

test_that("SimpleSBM_fit 'Bernoulli' model, undirected, no covariate, dyad sampling, MAR", {

  sampler_undirected_nocov$rNetwork(store = TRUE)

  adjMatrix  <- missSBM::observeNetwork(sampler_undirected_nocov$networkData, "dyad", 0.5)
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix)
  cls <- net$clustering(1:(2*Q))
  cl <- cls[[Q]]

  mySBM <- missSBM:::SimpleSBM_fit_noCov$new(net, cl)
  mySBM$doVEM()

  ## correctness
  expect_lt(error(mySBM$connectParam$mean, sampler_undirected_nocov$connectParam$mean), 0.075)
  expect_gt(ARI(mySBM$memberships       , sampler_undirected_nocov$memberships), 0.95)

})

test_that("SimpleSBM_fit 'Bernoulli' model, undirected, no covariate, node sampling, MAR", {

  sampler_undirected_nocov$rNetwork(store = TRUE)

  adjMatrix  <- missSBM::observeNetwork(sampler_undirected_nocov$networkData, "node", 0.5)
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix)
  cls <- net$clustering(1:(2*Q))
  cl <- cls[[Q]]

  mySBM <- missSBM:::SimpleSBM_fit_noCov$new(net, cl)
  mySBM$doVEM()

  ## correctness
  expect_lt(error(mySBM$connectParam$mean, sampler_undirected_nocov$connectParam$mean), 0.05)
  expect_gt(ARI(mySBM$memberships       , sampler_undirected_nocov$memberships), 0.95)

})

test_that("SimpleSBM_fit 'Bernoulli' model, directed, no covariate, dyad sampling, MAR", {

  sampler_directed_nocov$rNetwork(store = TRUE)

  adjMatrix  <- missSBM::observeNetwork(sampler_directed_nocov$networkData, "dyad", 0.5)
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix)
  cls <- net$clustering(1:(2*Q))
  cl <- cls[[Q]]

  mySBM <- missSBM:::SimpleSBM_fit_noCov$new(net, cl)
  mySBM$doVEM()

  ## correctness
  expect_lt(error(mySBM$connectParam$mean, sampler_directed_nocov$connectParam$mean), 0.05)
  expect_gt(ARI(mySBM$memberships       , sampler_directed_nocov$memberships), 0.95)

})

test_that("SimpleSBM_fit 'Bernoulli' model, directed, no covariate, node samping, MAR", {

  sampler_directed_nocov$rNetwork(store = TRUE)

  adjMatrix  <- missSBM::observeNetwork(sampler_directed_nocov$networkData, "node", 0.5)
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix)
  cls <- net$clustering(1:(2*Q))
  cl <- cls[[Q]]

  mySBM <- missSBM:::SimpleSBM_fit_noCov$new(net, cl)
  mySBM$doVEM()

  ## correctness
  expect_lt(error(mySBM$connectParam$mean, sampler_directed_nocov$connectParam$mean), 0.05)
  expect_gt(ARI(mySBM$memberships       , sampler_directed_nocov$memberships), 0.95)

})
