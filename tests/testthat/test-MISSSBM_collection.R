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
    clusterInit = cl,
    control = list(useCov = FALSE, trace = TRUE))

  ## control parameter for the VEM
  control <- list(threshold = 1e-2, maxIter = 50, fixPointIter = 3, trace = 0, iterates = 0)

  ## VEM Estimation on each element of the collection
  collection$estimate(control)
  expect_is(collection, "missSBM_collection")

  control$iterates  <- 1

  collection$explore(control, direction = "forward")
  expect_is(collection, "missSBM_collection")

  collection$explore(control, direction = "backward")
  expect_is(collection, "missSBM_collection")

  collection$explore(control, direction = "both")
  expect_is(collection, "missSBM_collection")
})

