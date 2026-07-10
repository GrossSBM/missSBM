context("test missSBM_fit's split/merge methods")

set.seed(1)
N <- 150; Q <- 4
sampler <- sbm::SimpleSBM$new("bernoulli", N, FALSE, rep(1/Q, Q), list(mean = diag(.4, Q) + .05))
sampler$rMemberships(store = TRUE)
sampler$rNetwork(store = TRUE)
net <- missSBM:::partlyObservedNetwork$new(sampler$networkData)
fit <- missSBM:::missSBM_fit$new(net, "dyad", sampler$memberships, TRUE)
fit$doVEM(list(threshold = 1e-3, maxIter = 50, fixPointIter = 3, trace = FALSE))

test_that("split() returns a new missSBM_fit with one more block, self is untouched", {
  candidate <- fit$split(1)
  expect_s3_class(candidate, "missSBM_fit")
  expect_equal(candidate$fittedSBM$nbBlocks, fit$fittedSBM$nbBlocks + 1)
  expect_equal(fit$fittedSBM$nbBlocks, Q) # self untouched
})

test_that("split(in_place = TRUE) mutates self instead of returning a new object", {
  fit_clone <- missSBM:::missSBM_fit$new(net, "dyad", sampler$memberships, TRUE)
  fit_clone$doVEM(list(threshold = 1e-3, maxIter = 50, fixPointIter = 3, trace = FALSE))
  q_before <- fit_clone$fittedSBM$nbBlocks
  out <- fit_clone$split(1, in_place = TRUE)
  expect_identical(out, fit_clone)
  expect_equal(fit_clone$fittedSBM$nbBlocks, q_before + 1)
})

test_that("merge() returns a new missSBM_fit with one fewer block, self is untouched", {
  candidate <- fit$merge(c(1, 2))
  expect_s3_class(candidate, "missSBM_fit")
  expect_equal(candidate$fittedSBM$nbBlocks, fit$fittedSBM$nbBlocks - 1)
  expect_equal(fit$fittedSBM$nbBlocks, Q) # self untouched
})

test_that("candidates_split() returns one trial-fitted candidate per splittable cluster", {
  control <- list(threshold = 1e-3, maxIter = 50, fixPointIter = 3, trace = FALSE)
  candidates <- fit$candidates_split(control)
  expect_true(length(candidates) > 0)
  expect_true(all(sapply(candidates, function(m) inherits(m, "missSBM_fit"))))
  expect_true(all(sapply(candidates, function(m) m$fittedSBM$nbBlocks) == Q + 1))
  expect_true(all(sapply(candidates, function(m) is.finite(m$ICL))))
})

test_that("candidates_merge() returns one trial-fitted candidate per (capped) cluster pair", {
  control <- list(threshold = 1e-3, maxIter = 50, fixPointIter = 3, trace = FALSE)
  candidates <- fit$candidates_merge(control)
  expect_equal(length(candidates), choose(Q, 2))
  expect_true(all(sapply(candidates, function(m) m$fittedSBM$nbBlocks) == Q - 1))

  ## max_candidates caps the (combinatorial) number of pairs tried
  capped <- fit$candidates_merge(control, max_candidates = 2)
  expect_equal(length(capped), 2)
})

test_that("estimateMissSBM()'s default exploration still finds a sensible number of blocks", {
  adj <- missSBM::observeNetwork(sampler$networkData, "dyad", 0.85)
  collection <- estimateMissSBM(adj, vBlocks = 2:6, sampling = "dyad",
                                 control = list(trace = FALSE, iterates = 1, exploration = "both"))
  expect_s3_class(collection, "missSBM_collection")
  expect_true(all(is.finite(collection$ICL)))
})
