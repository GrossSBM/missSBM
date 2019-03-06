context("test-consistency-of-misssbm-with-dyad-sampling")

library(aricode)

referenceResults <- readRDS(system.file("extdata", "referenceResults.rds", package = "missSBM"))
refAlgo <- referenceResults$dyad
truth   <- referenceResults$true_sbm

tol_ref   <- 1e-3
tol_truth <- 1e-2
tol_ARI   <- .9
error <- function(beta1, beta2) {
  sum((sort(beta1) - sort(beta2))^2)/truth$Q
}


test_that("check consistency against Tim's code for dyad sampling", {

  # for (sampling in c("dyad", "node", "block", "") )

  missSBM_out <- missSBM::inferSBM(
    adjacencyMatrix = referenceResults$dyad$Y,
    vBlocks = referenceResults$true_sbm$Q,
    sampling = "dyad",
    control_VEM = list(threshold = 1e-4, maxIter = 200, fixPointIter = 5, trace = FALSE)
  )
  newAlgo <- missSBM_out[[1]]

  ## alpha
  expect_lt(error(newAlgo$fittedSBM$mixtureParam, refAlgo$alpha), tol_ref)
  expect_lt(error(newAlgo$fittedSBM$mixtureParam, truth$alpha), tol_truth)
  expect_lt(error(refAlgo$alpha, truth$alpha), tol_truth)

  ## pi
  expect_lt(error(newAlgo$fittedSBM$connectParam, refAlgo$pi), tol_ref)
  expect_lt(error(newAlgo$fittedSBM$connectParam, truth$pi), tol_truth)
  expect_lt(error(refAlgo$pi, truth$pi), tol_truth)

  ## clustering
  expect_gt(ARI(newAlgo$fittedSBM$memberships, refAlgo$cl), tol_ARI)
  expect_gt(ARI(newAlgo$fittedSBM$memberships, truth$cl), tol_truth)
  expect_gt(ARI(refAlgo$cl, truth$cl), tol_truth)

})
