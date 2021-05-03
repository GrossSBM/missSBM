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
  control <- list(threshold = 1e-2, maxIter = 50, fixPointIter = 3, cores = 1, trace = 0, iterates = 0, smoothing = "both")

  ## VEM Estimation on each element of the collection
  collection$estimate(control)
  expect_is(collection, "missSBM_collection")

  control$iterates  <- 1

  control$smoothing <- "forward"
  collection$smooth(control)
  expect_is(collection, "missSBM_collection")

  control$smoothing <- "backward"
  collection$smooth(control)
  expect_is(collection, "missSBM_collection")

  control$smoothing <- "both"
  collection$smooth(control)
  expect_is(collection, "missSBM_collection")
})

