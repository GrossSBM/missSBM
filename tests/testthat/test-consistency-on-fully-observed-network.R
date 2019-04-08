context("test-consistency-on-fully-observed-network")

library(igraph)

test_that("SBM_fit and missSBMfit are coherent", {
  data("war_graphs")

  ## adjacency matrix without missing values
  A <- war_graphs$beligerent %>%  as_adj(sparse = FALSE)

  ## coherence of sampledNetwork object
  sampledNet <- sampledNetwork$new(A)
  expect_equal(A, sampledNet$adjacencyMatrix)
  expect_equal(ncol(A), sampledNet$nNodes)
  expect_equal(ncol(A) * (ncol(A) - 1)/2, sampledNet$nDyads)
  expect_equal(length(sampledNet$missingDyads), 0)
  expect_equal(sampledNet$dyads, sampledNet$observedDyads)
  expect_equal(rep(TRUE, ncol(A)), sampledNet$observedNodes)

  ## initial clustering
  Q <- 3
  cl0   <- missSBM:::init_hierarchical(A, Q)

  control <- list(threshold = 1e-4, maxIter = 200, fixPointIter = 5, trace = TRUE)

  ## using SBM_fit class
  my_SBM <- missSBM:::SBM_fit_nocovariate$new(adjacencyMatrix = A, clusterInit = cl0)
  my_SBM$doVEM(control$threshold, control$maxIter, control$fixPointIter, control$trace)
  my_SBM$vICL

  ## using missSBM_fit class
  my_missSBM <- missSBM:::missSBM_fit$new(sampledNet = sampledNet, Q, netSampling = "node", clusterInit = cl0)
  my_missSBM$doVEM(control)
  my_missSBM$fittedSBM$vICL


  ## using missSBM_collection class
  my_collection <- missSBM_collection$new(
      sampledNet  = sampledNet,
      vBlocks     = Q,
      sampling    = "node",
      clusterInit = cl0,
      cores       = 1,
      trace       = TRUE
  )
  my_collection$estimate(control)

  expect_equivalent(my_SBM, my_missSBM$fittedSBM)
  expect_equivalent(my_SBM, my_collection$bestModel$fittedSBM)
  expect_lt(my_SBM$vICL, my_collection$ICL) ## different due an addition df for sampling
})
