context("test-misssbm_collection")

set.seed(1890718)
### A SBM model : ###
N <- 100
Q <- 3
pi <- rep(1, Q)/Q           # block proportion
theta <- diag(.45, Q, Q) + .05 # connectivity matrix
directed <- FALSE              # if the network is directed or not

### Draw a SBM model
mySBM <- missSBM::simulate(N, pi, theta, directed) # simulation of ad Bernoulli non-directed SBM
A <- mySBM$netMatrix             # the adjacency matrix

test_that("missSBMcollection works", {

  adjMatrix  <- missSBM::sample(A, "dyad", .5, clusters = mySBM$memberships)
  sampledNet <- missSBM:::sampledNetwork$new(adjMatrix)

  ## Instantiate the collection of missSBM_fit
  collection <- missSBM_collection$new(
    sampledNet  = sampledNet,
    vBlocks     = 1:4,
    sampling    = "dyad",
    clusterInit = 'hierarchical', 1, TRUE, TRUE)

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

  adjMatrix  <- missSBM::sample(A, "dyad", .5, clusters = mySBM$memberships)
  sampledNet <- missSBM:::sampledNetwork$new(adjMatrix)


  ## Instantiate the collection of missSBM_fit
  collection <- missSBM_collection$new(
    sampledNet  = sampledNet,
    vBlocks     = 1:4,
    sampling    = "dyad",
    clusterInit = 'hierarchical', 1, TRUE, TRUE)

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
