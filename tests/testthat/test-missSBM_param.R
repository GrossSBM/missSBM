context("test missSBM_param() and the deprecated list()-control gate")

test_that("missSBM_param() returns a classed list with the documented defaults", {
  ctrl <- missSBM_param()
  expect_s3_class(ctrl, "missSBM_param")
  expect_equal(ctrl$threshold, 1e-2)
  expect_equal(ctrl$maxIter, 50)
  expect_equal(ctrl$fixPointIter, 3)
  expect_equal(ctrl$imputation, "median")
  expect_identical(ctrl$similarity, l1_similarity)
  expect_true(ctrl$useCov)
  expect_null(ctrl$clusterInit)
  expect_true(ctrl$polish)
  expect_equal(ctrl$iterates, 1)
  expect_equal(ctrl$maxMergeCandidates, 30)
  expect_true(ctrl$stopOnDegenerate)
  expect_equal(ctrl$maxConsecutiveDegenerate, 2)
  expect_false(ctrl$warmChain)
  expect_true(ctrl$trace)
})

test_that("missSBM_param() overrides only the arguments passed", {
  ctrl <- missSBM_param(trace = FALSE, iterates = 0)
  expect_false(ctrl$trace)
  expect_equal(ctrl$iterates, 0)
  expect_equal(ctrl$threshold, 1e-2) # untouched default
})

test_that("missSBM_param()'s imputation is validated via match.arg()", {
  expect_equal(missSBM_param(imputation = "average")$imputation, "average")
  expect_error(missSBM_param(imputation = "not-a-strategy"))
})

test_that("missSBM_param() rejects unknown argument names", {
  expect_error(missSBM_param(notAnOption = TRUE), "unused argument")
})

test_that("estimateMissSBM() errors informatively when control is a plain list()", {
  adjMatrix <- missSBM::observeNetwork(
    sbm::sampleSimpleSBM(30, c(1/2, 1/2), list(mean = diag(.4, 2) + .05))$networkData,
    "dyad", 0.9)
  expect_error(
    estimateMissSBM(adjMatrix, 1:2, "dyad", control = list(trace = FALSE)),
    "missSBM_param"
  )
})

test_that("estimateMissSBM() accepts the default (no control argument) and missSBM_param() explicitly", {
  adjMatrix <- missSBM::observeNetwork(
    sbm::sampleSimpleSBM(30, c(1/2, 1/2), list(mean = diag(.4, 2) + .05))$networkData,
    "dyad", 0.9)
  expect_no_error(estimateMissSBM(adjMatrix, 1:2, "dyad"))
  expect_no_error(estimateMissSBM(adjMatrix, 1:2, "dyad", control = missSBM_param(trace = FALSE)))
})
