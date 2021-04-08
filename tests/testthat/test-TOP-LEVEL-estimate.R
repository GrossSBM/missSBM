context("test-test-top-level-function-misssbm")

library(aricode)

set.seed(1890719)
### A SBM model : ###
N <- 80
Q <- 3
pi <- rep(1, Q)/Q           # block proportion
theta <- list(mean = diag(.45,Q) + .05 ) # connectivity matrix

### Draw a SBM model
mySBM <- sbm::sampleSimpleSBM(N, pi, theta) # simulation of ad Bernoulli non-directed SBM
A <- mySBM$networkData # the adjacency matrix
### UGLY FIX - FIXME
diag(A) <- 0

test_that("missSBM and class missSBM-fit are coherent", {

  l_psi <- list(
    "dyad" = c(.6),
    "node" = c(.6),
    "double-standard" = c(0.4, 0.8),
    "block-node" = c(.3, .8, .5),
    "block-dyad" = mySBM$connectParam$mean,
    "degree" = c(.01, .01)
  )

  for (k in seq_along(l_psi)) {

    sampling <- names(l_psi)[k]
    adjMatrix <- missSBM::observeNetwork(A, sampling, l_psi[[k]], clusters = mySBM$memberships)
    partlyObservedNet <- missSBM:::partlyObservedNetwork$new(adjMatrix)
    cl0 <- partlyObservedNet$clustering(Q)[[Q]]

    ## control parameter for the VEM
    control <- list(threshold = 1e-4, maxIter = 200, fixPointIter = 5, trace = 0)

    ## Perform inference with internal classes
    missSBM <- missSBM:::missSBM_fit$new(partlyObservedNet, sampling, cl0, TRUE)
    out_missSBM <- missSBM$doVEM(control)

    ## Perform inference with the top level function
    collection <- estimateMissSBM(
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
    "block-dyad" = mySBM$connectParam$mean#,
#    "degree" = c(.01, .01)
  )

  for (k in seq_along(l_psi)) {

    sampling <- names(l_psi)[k]
    adjMatrix <- missSBM::observeNetwork(A, sampling, l_psi[[k]], clusters = mySBM$memberships)

    control <- list(threshold = 1e-4, maxIter = 200, fixPointIter = 5, trace = 0)

    cat("\nSampling:", sampling)

    ## Perform inference with the top level function
    collection <- estimateMissSBM(
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
