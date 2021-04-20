context("test-misssbm_collection")

source("utils_test.R", local =TRUE)
sampler_undirected_nocov$rNetwork(store = TRUE)

test_that("missSBMcollection works", {

  adjMatrix  <- missSBM::observeNetwork(sampler_undirected_nocov$networkData, "dyad", .5, clusters = sampler_undirected_nocov$memberships)
  partlyObservedNet <- missSBM:::partlyObservedNetwork$new(adjMatrix)
  cl <- partlyObservedNet$clustering(1:3)

  ## Instantiate the collection of missSBM_fit
  collection <- missSBM_collection$new(
    partlyObservedNet  = partlyObservedNet,
    sampling    = "dyad",
    clusterInit = cl, 1, TRUE, TRUE)

  ## control parameter for the VEM
  control <- list(threshold = 1e-4, maxIter = 200, fixPointIter = 5, cores = 1, trace = 0)

  ## VEM Estimation on each element of the collection
  collection$estimate(control)
  expect_is(collection, "missSBM_collection")

  smooth(collection, "forward")
  expect_is(collection, "missSBM_collection")

  smooth(collection, "backward")
  expect_is(collection, "missSBM_collection")

  smooth(collection, "both")
  expect_is(collection, "missSBM_collection")
})

test_that("More smoothing tests", {

  adjMatrix  <- missSBM::observeNetwork(sampler_undirected_nocov$networkData, "dyad", .5)
  partlyObservedNet <- missSBM:::partlyObservedNetwork$new(adjMatrix)
  cl <- partlyObservedNet$clustering(1:4)

  ## Instantiate the collection of missSBM_fit
  collection <- missSBM_collection$new(
    partlyObservedNet  = partlyObservedNet,
    sampling    = "dyad",
    clusterInit = cl, 1, TRUE, TRUE)

  ## control parameter for the VEM
  control <- list(threshold = 1e-4, maxIter = 200, fixPointIter = 5, cores = 1, trace = 0)

  ## VEM Estimation on each element of the collection
  collection$estimate(control)
  expect_is(collection, "missSBM_collection")

  smooth(collection, "forward", control = list(iterates = 2))
  expect_is(collection, "missSBM_collection")

  smooth(collection, "backward", control = list(iterates = 2))
  expect_is(collection, "missSBM_collection")

  smooth(collection, "both", control = list(iterates = 2))
  expect_is(collection, "missSBM_collection")
})
