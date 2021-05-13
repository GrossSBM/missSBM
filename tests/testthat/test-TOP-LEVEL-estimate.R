context("test-test-top-level-function-misssbm")

N_nocov <- 80
source("utils_test.R", local = TRUE)
sampler_undirected_nocov$rNetwork(store = TRUE)
A  <- sampler_undirected_nocov$networkData
cl <- sampler_undirected_nocov$memberships

l_psi <- list(
  "dyad" = c(.7),
  "node" = c(.7),
  "double-standard" = c(0.4, 0.8),
  "block-node" = c(.3, .8, .5),
  "block-dyad" = sampler_undirected_nocov$connectParam$mean#,
  #    "degree" = c(.01, .01)
)

test_that("missSBM and class missSBM-fit are coherent", {

  for (k in seq_along(l_psi)) {

    sampling <- names(l_psi)[k]
    adjMatrix <- missSBM::observeNetwork(A, sampling, l_psi[[k]], clusters = cl)
    partlyObservedNet <- missSBM:::partlyObservedNetwork$new(adjMatrix)
    cl0 <- partlyObservedNet$clustering(1:Q)[[Q]]

    ## control parameter for the VEM
    control <- list(threshold = 1e-3, maxIter = 100, fixPointIter = 5, trace = 0, exploration = "none")

    ## Perform inference with internal classes
    missSBM <- missSBM:::missSBM_fit$new(partlyObservedNet, sampling, cl0, TRUE)
    out_missSBM <- missSBM$doVEM(control)

    ## Perform inference with the top level function
    collection <- estimateMissSBM(
      adjacencyMatrix = adjMatrix,
      vBlocks     = 1:3,
      sampling    = sampling,
      control     = control
    )

    expect_is(collection, "missSBM_collection")
    expect_gte(aricode::ARI(collection$models[[Q]]$fittedSBM$memberships, missSBM$fittedSBM$memberships), 1)

  }

})

