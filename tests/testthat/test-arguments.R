context("test the checks of the arguments of estimateMissSBM() and observeNetwork()")

set.seed(1)
full <- sbm::sampleSimpleSBM(30, c(.5, .5), list(mean = diag(.4, 2) + .05))$networkData
A    <- observeNetwork(full, "dyad", .8)
ctrl <- missSBM_param(trace = FALSE, iterates = 0, polish = FALSE)

test_that("estimateMissSBM() rejects invalid arguments with informative messages", {
  expect_error(estimateMissSBM(A[, -1], 1:2, "dyad", control = ctrl), "square")
  expect_error(estimateMissSBM(as.data.frame(A), 1:2, "dyad", control = ctrl), "must be a matrix")
  expect_error(estimateMissSBM(A, 1:2, "foo", control = ctrl), "sampling must be one of")
  expect_error(estimateMissSBM(A, 1:2, c("dyad", "node"), control = ctrl), "sampling must be one of")
  expect_error(estimateMissSBM(A, 1:2, "covar-node", control = ctrl), "needs covariates")
  for (vBlocks in list(0:2, 2.5, 30, "2"))
    expect_error(estimateMissSBM(A, vBlocks, "dyad", control = ctrl), "vBlocks must be")
  expect_error(estimateMissSBM(A, 1:2, "dyad", covariates = rnorm(30), control = ctrl), "must be a list")
  expect_error(estimateMissSBM(A, 1:2, "dyad", covariates = list(rnorm(10)), control = ctrl), "size N = 30")
  expect_error(estimateMissSBM(A, 1:2, "block-node", covariates = list(rnorm(30)), control = ctrl), "only available")
  ctrl_cl <- ctrl; ctrl_cl$clusterInit <- list(rep(1:2, 15))
  expect_error(estimateMissSBM(A, 1:2, "dyad", control = ctrl_cl), "one clustering per number of blocks")
  ctrl_cl$clusterInit <- list(rep(1, 30), rep(1:2, 10))
  expect_error(estimateMissSBM(A, 1:2, "dyad", control = ctrl_cl), "N = 30 block labels")
})

test_that("estimateMissSBM() sorts clusterInit with vBlocks", {
  ctrl_cl <- ctrl; ctrl_cl$clusterInit <- list(rep(1:2, 15), rep(1, 30))
  collection <- expect_warning(estimateMissSBM(A, 2:1, "dyad", control = ctrl_cl), "sorted")
  expect_equal(sapply(collection$models, function(m) m$fittedSBM$nbBlocks), 1:2)
})

test_that("observeNetwork() rejects invalid arguments with informative messages", {
  expect_error(observeNetwork(full[, -1], "dyad", .5), "square")
  expect_error(observeNetwork(full, "foo", .5), "sampling must be one of")
  expect_error(observeNetwork(full, "covar-node", 1), "needs covariates")
  expect_error(observeNetwork(full, "covar-node", 1, covariates = list(rnorm(10))), "size N = 30")
})
