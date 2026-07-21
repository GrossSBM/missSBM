context("test missSBM_collection's handling of degenerate (collapsed-class) models")

source("utils_test.R", local = TRUE)
sampler_undirected_nocov$rNetwork(store = TRUE)

## builds a small, cheaply-fitted collection at nbBlocks = 1:5 and force-degenerates the fits at
## `degenerate_at` (by zeroing out their highest class in tau, as in test-missSBM_fit-polish.R),
## bypassing collection$estimate()'s automatic repair() so the forced state is under our control
build_forced_collection <- function(degenerate_at = integer(0)) {
  adjMatrix <- missSBM::observeNetwork(sampler_undirected_nocov$networkData, "dyad", .9,
                                        clusters = sampler_undirected_nocov$memberships)
  partlyObservedNet <- missSBM:::partlyObservedNetwork$new(adjMatrix)
  cl <- partlyObservedNet$clustering(1:5)
  collection <- missSBM_collection$new(
    partlyObservedNet = partlyObservedNet, sampling = "dyad", clusterInit = cl,
    control = list(useCov = FALSE, trace = FALSE))

  control <- list(threshold = 1e-2, maxIter = 50, fixPointIter = 3, trace = FALSE)
  for (Q in degenerate_at) {
    fit <- collection$models[[Q]]
    fit$doVEM(control)
    state <- fit$fittedSBM$get_state()
    state$Z[, Q] <- 0
    state$Z <- state$Z / rowSums(state$Z)
    fit$fittedSBM$set_state(state)
  }
  for (Q in setdiff(1:5, degenerate_at)) collection$models[[Q]]$doVEM(control)

  collection
}

test_that("occupiedBlocks/degenerate reflect a forced degenerate state", {
  collection <- build_forced_collection(degenerate_at = c(4, 5))
  expect_equal(collection$vBlocks, 1:5)
  expect_equal(collection$occupiedBlocks, c(1, 2, 3, 3, 4))
  expect_equal(collection$degenerate, c(FALSE, FALSE, FALSE, TRUE, TRUE))
})

test_that("bestModel never picks a degenerate model when a non-degenerate one exists", {
  collection <- build_forced_collection(degenerate_at = 5)
  ## the collection's own selection, computed only over non-degenerate models, must agree with
  ## a from-scratch argmin restricted the same way
  icl <- collection$ICL
  expected <- collection$models[[which(!collection$degenerate)[which.min(icl[!collection$degenerate])]]]
  expect_identical(collection$bestModel, expected)
  expect_false(collection$bestModel$degenerate)
})

test_that("bestModel falls back to a degenerate model when every model is degenerate", {
  collection <- build_forced_collection(degenerate_at = 1:5)
  ## nothing better available: must not error, must still return some model
  expect_s3_class(collection$bestModel, "missSBM_fit")
})

test_that("flag_degenerate_tail() caps forwardLimit only for a persistent run, and lifts it once cured", {
  collection <- build_forced_collection(degenerate_at = c(4, 5))
  priv <- environment(collection$explore)$private

  expect_warning(priv$flag_degenerate_tail(2), "does not appear to support")
  expect_equal(priv$forwardLimit, 3) # run of 2 consecutive degenerate models first complete at index 5 => limit = 5 - 2

  ## re-running with the same persistent pattern must not warn again (already flagged)
  expect_no_warning(priv$flag_degenerate_tail(2))

  ## simulate a cure (e.g. a later pass fixed it) by forcing full occupancy directly -- avoids
  ## depending on whether repair()'s VEM refit actually converges to a valid Q-block clustering
  ## for this (possibly over-specified) synthetic data, which is not what this test is about
  for (Q in c(4, 5)) {
    fit <- collection$models[[Q]]
    labels <- missSBM:::repair_empty_classes(fit$fittedSBM$memberships, Q)
    state <- fit$fittedSBM$get_state()
    state$Z <- missSBM:::clustering_indicator(labels)
    fit$fittedSBM$set_state(state)
  }
  priv$flag_degenerate_tail(2)
  expect_null(priv$forwardLimit)
})

test_that("estimateMissSBM()'s exploration stops growing forward into a persistently degenerate range", {
  skip_on_cran() # network-scale reproduction, too slow for CRAN's routine checks

  frenchblog <- igraph::delete_vertices(missSBM::frenchblog2007, which(igraph::degree(missSBM::frenchblog2007) == 0))
  frenchblog <- igraph::delete_vertices(frenchblog, 61:igraph::vcount(frenchblog))
  blog <- igraph::as_adjacency_matrix(frenchblog, sparse = FALSE)

  set.seed(3052008)
  sbm_full <- estimateMissSBM(blog, 1:6, "node", control = missSBM_param(trace = FALSE, iterates = 1, polish = FALSE))
  samplingParameters <- ifelse(sbm_full$bestModel$fittedSBM$blockProp < 0.1, 0.2, 0.8)
  blog_obs <- observeNetwork(blog, sampling = "block-node", parameters = samplingParameters,
                              clusters = sbm_full$bestModel$fittedSBM$memberships)

  expect_warning(
    res <- estimateMissSBM(blog_obs, 1:16, "block-node", control = missSBM_param(trace = FALSE, iterates = 1)),
    "does not appear to support"
  )

  ## a genuinely over-specified tail: at least the top model is flagged degenerate
  expect_true(res$degenerate[length(res$degenerate)])
  ## bestModel must not be one of the degenerate ones
  expect_false(res$bestModel$degenerate)
})
