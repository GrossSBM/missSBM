context("test-consistency-on-fully-observed-network")

library(igraph)

test_that("SimpleSBM_fit_missSBM and missSBMfit are coherent", {
  data("war")

  ## adjacency matrix without missing values
  A <- war$belligerent %>%  as_adj(sparse = FALSE)

  ## coherence of partlyObservedNetwork object
  partlyObservedNet <- missSBM:::partlyObservedNetwork$new(A)
  expect_equal(A, partlyObservedNet$networkData)
  expect_equal(ncol(A), partlyObservedNet$nbNodes)
  expect_equal(ncol(A) * (ncol(A) - 1)/2, partlyObservedNet$nbDyads)
  expect_equal(length(partlyObservedNet$missingDyads), 0)
  expect_equal(rep(TRUE, ncol(A)), partlyObservedNet$observedNodes)

  ## initial clustering
  Q <- 3
  cl0   <- missSBM:::init_hierarchical(A, Q)

  control <- list(threshold = 1e-4, maxIter = 200, fixPointIter = 5, cores = 1, trace = 1)

  ## using SBM_fit class
  my_SBM <- missSBM:::SimpleSBM_fit_missSBM$new(adjacencyMatrix = A, clusterInit = cl0)
  my_SBM$doVEM(control$threshold, control$maxIter, control$fixPointIter, control$trace)
  my_SBM$reorder()
  my_SBM$ICL

  ## using missSBM_fit class
  my_missSBM <- missSBM:::missSBM_fit$new(partlyObservedNet = partlyObservedNet, Q, netSampling = "node", clusterInit = cl0, useCov = TRUE)
  my_missSBM$doVEM(control)
  my_missSBM$fittedSBM$ICL


  ## using missSBM_collection class
  my_collection <- missSBM_collection$new(
      partlyObservedNet  = partlyObservedNet,
      vBlocks     = Q,
      sampling    = "node",
      clusterInit = cl0,
      cores       = 1,
      trace       = TRUE,
      useCov      = TRUE
  )
  my_collection$estimate(control)

  # expect_equivalent(my_SBM, my_missSBM$fittedSBM)
  # expect_equivalent(my_SBM, my_collection$bestModel$fittedSBM)
  expect_lt(my_SBM$ICL, my_collection$ICL) ## different due an addition df for sampling
})
