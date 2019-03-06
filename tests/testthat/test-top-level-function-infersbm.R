context("test-test-top-level-function-infersbm")

library(aricode)

set.seed(1890718)
### A SBM model : ###
N <- 300
Q <- 5
alpha <- rep(1,Q)/Q       # mixture parameter
pi <- diag(.45,Q) + .05   # connectivity matrix
directed <- FALSE         # if the network is directed or not

### Draw a SBM model
mySBM <- simulateSBM(N, alpha, pi, directed) # simulation of ad Bernoulli non-directed SBM
A <- mySBM$adjacencyMatrix             # the adjacency matrix

test_that("inferSBM and class missSBM-fit are coherent", {

  psi <- 0.3
  sampledNet <- samplingSBM(A, "dyad", psi)

  ## control parameter for the VEM
  control <- list(threshold = 1e-4, maxIter = 200, fixPointIter = 5, trace = FALSE)

  ## Perform inference with internal classes
  missSBM <- missSBM:::missingSBM_fit$new(sampledNet, Q, "dyad")
  out_missSBM <- missSBM$doVEM(control)

  ## Perform inference with the top level function
  collection <- inferSBM(
    adjacencyMatrix = sampledNet$adjacencyMatrix,
    vBlocks         = Q,
    sampling        = "dyad",
    smoothing       = "none",
    control_VEM     = control
  )

  expect_true(is.missSBMcollection(collection))
  expect_equivalent(collection[[1]], missSBM)

})

test_that("inferSBM with a collection of models", {

  psi <- 0.75
  sampledNet <- samplingSBM(A, "dyad", psi)

  ## control parameter for the VEM
  control <- list(threshold = 1e-4, maxIter = 200, fixPointIter = 5, trace = FALSE)

  ## Perform inference with the top level function
  collection <- inferSBM(
    adjacencyMatrix = sampledNet$adjacencyMatrix,
    vBlocks         = 1:7,
    sampling        = "dyad",
    smoothing       = "none",
    control_VEM     = control
  )

  expect_true(is.missSBMcollection(collection))
  expect_equal(getBestModel(collection)$fittedSBM$nBlocks, Q)
  expect_true(which.min(ICL(collection)) == Q)

  expect_true(is.data.frame(optimizationStatus(collection)))
})

# test_that("infer SBM with node sampling works", {
#
# })
#
# test_that("infer SBM with double standard sampling works", {
#
# })
#
# test_that("infer SBM with block sampling works", {
#
# })
#
# test_that("infer SBM with degree sampling works", {
#
# })
