context("test-test-top-level-function-misssbm")

library(aricode)

set.seed(1890719)
### A SBM model : ###
N <- 80
Q <- 3
alpha <- rep(1, Q)/Q           # block proportion
theta <- diag(.45, Q, Q) + .05 # connectivity matrix
directed <- FALSE         # if the network is directed or not

### Draw a SBM model
mySBM <- missSBM::simulate(N, alpha, theta, directed) # simulation of ad Bernoulli non-directed SBM
A <- mySBM$adjacencyMatrix # the adjacency matrix

test_that("missSBM and class missSBM-fit are coherent", {

  l_psi <- list(
    "dyad" = c(.6),
    "node" = c(.6),
    "double-standard" = c(0.4, 0.8),
    "block-node" = c(.3, .8, .5),
    "block-dyad" = mySBM$connectParam,
    "degree" = c(.01, .01)
  )

  for (k in seq_along(l_psi)) {

    sampling <- names(l_psi)[k]
    adjMatrix <- missSBM::sample(A, sampling, l_psi[[k]], clusters = mySBM$memberships)
    sampledNet <- missSBM:::sampledNetwork$new(adjMatrix)
    ## control parameter for the VEM
    control <- list(threshold = 1e-4, maxIter = 200, fixPointIter = 5, trace = 0, clusterInit = "spectral")

    ## Perform inference with internal classes
    missSBM <- missSBM:::missSBM_fit$new(sampledNet, Q, sampling, "spectral", TRUE)
    out_missSBM <- missSBM$doVEM(control)

    ## Perform inference with the top level function
    collection <- missSBM::estimate(
      adjacencyMatrix = adjMatrix,
      vBlocks     = Q,
      sampling    = sampling,
      control     = control
    )

    expect_is(collection, "missSBM_collection")
    expect_gte(aricode::ARI(collection$models[[1]]$fittedSBM$memberships, missSBM$fittedSBM$memberships), 1)

  }

})

test_that("missSBM with a collection of models", {

  l_psi <- list(
    "dyad" = c(.75),
    "node" = c(.75),
    "double-standard" = c(0.4, 0.8),
    "block-node" = c(.3, .8, .5),
    "block-dyad" = mySBM$connectParam#,
#    "degree" = c(.01, .01)
  )

  for (k in seq_along(l_psi)) {

    sampling <- names(l_psi)[k]
    adjMatrix <- missSBM::sample(A, sampling, l_psi[[k]], clusters = mySBM$memberships)

    control <- list(threshold = 1e-4, maxIter = 200, fixPointIter = 5, trace = 0)

    cat("\nSampling:", sampling)

    ## Perform inference with the top level function
    collection <- missSBM::estimate(
      adjacencyMatrix = adjMatrix,
      vBlocks     = 1:5,
      sampling    = sampling,
      control     = control
    )

    expect_is(collection, "missSBM_collection")
    expect_true(is.data.frame(collection$optimizationStatus))

    expect_true(is.data.frame(collection$optimizationStatus))

    #expect_equal(collection$bestModel$fittedSBM$nbBlocks, Q)
    #expect_true(which.min(collection$ICL) == Q)

  }
})
