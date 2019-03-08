context("test-consistency-of-misssbm-with-dyad-sampling")

library(aricode)

referenceResults <- readRDS(system.file("extdata", "referenceResults.rds", package = "missSBM"))
truth   <- referenceResults$true_sbm

tol_ref   <- 1e-3
tol_truth <- 1e-2
tol_ARI   <- .9
error <- function(beta1, beta2) {
  sum((sort(beta1) - sort(beta2))^2)/truth$Q
}

test_that("check consistency against Tim's code for dyad, node, double standard and block sampling", {


  for (sampling in c("dyad", "node", "twoStd", "block")) {

    refAlgo <- referenceResults[[sampling]]

    sampling <- ifelse(sampling == "twoStd", "double-standard", sampling)
    sampling <- ifelse(sampling == "block", "block-node", sampling)
    missSBM_out <- missSBM::missSBM(
      adjacencyMatrix = refAlgo$Y,
      vBlocks = truth$Q,
      sampling = sampling,
      trace = FALSE
    )
    newAlgo <- missSBM_out[[1]]

    ## alpha
    expect_lt(error(sort(newAlgo$fittedSBM$mixtureParam), sort(refAlgo$alpha)), tol_ref)
    expect_lt(error(sort(newAlgo$fittedSBM$mixtureParam), sort(truth$alpha)), tol_truth)
    expect_lt(error(sort(refAlgo$alpha), sort(truth$alpha)), tol_truth)

    ## pi
    expect_lt(error(newAlgo$fittedSBM$connectParam, refAlgo$pi), tol_ref)
    expect_lt(error(newAlgo$fittedSBM$connectParam, truth$pi), tol_truth)
    expect_lt(error(refAlgo$pi, truth$pi), tol_truth)

    ## clustering
    expect_gt(ARI(newAlgo$fittedSBM$memberships, refAlgo$cl), tol_ARI)
    expect_gt(ARI(newAlgo$fittedSBM$memberships, truth$cl), tol_truth)
    expect_gt(ARI(refAlgo$cl, truth$cl), tol_truth)

    if (!(sampling %in% c("dyad", "node"))) {
      ## psi
      expect_lt(error(sort(newAlgo$fittedSampling$parameters), sort(refAlgo$psi_chap)), tol_ref)
      expect_lte(
        error(sort(newAlgo$fittedSampling$parameters), sort(refAlgo$psi_true)),
        error(sort(refAlgo$psi_chap), sort(refAlgo$psi_true)),
        tol_ref
      )
    }
  }
})

# test_that("check consistency against Tim's code for dyad and node sampling with covariates", {
#
#
#   for (sampling in c("covDyad", "covNode")) {
#
#     refAlgo <- referenceResults[[sampling]]
#
#     missSBM_out <- missSBM::missSBM(
#       adjacencyMatrix = refAlgo$Y,
#       vBlocks = truth$Q,
#       sampling = ifelse(sampling == "covDyad", "dyad", "node"),
#       trace = FALSE,
#       covarMatrix =
#     )
#     newAlgo <- missSBM_out[[1]]
#
#     ## alpha
#     expect_lt(error(sort(newAlgo$fittedSBM$mixtureParam), sort(refAlgo$alpha)), tol_ref)
#     expect_lt(error(sort(newAlgo$fittedSBM$mixtureParam), sort(truth$alpha)), tol_truth)
#     expect_lt(error(sort(refAlgo$alpha), sort(truth$alpha)), tol_truth)
#
#     ## pi
#     expect_lt(error(newAlgo$fittedSBM$connectParam, missSBM:::logistic(refAlgo$gamma_covDyad)), tol_ref)
#     expect_lt(error(newAlgo$fittedSBM$connectParam, truth$pi), tol_truth)
#     expect_lt(error(refAlgo$pi, truth$pi), tol_truth)
#
#     ## clustering
#     expect_gt(ARI(newAlgo$fittedSBM$memberships, refAlgo$cl), tol_ARI)
#     expect_gt(ARI(newAlgo$fittedSBM$memberships, truth$cl), tol_truth)
#     expect_gt(ARI(refAlgo$cl, truth$cl), tol_truth)
#
#     if (!(sampling %in% c("dyad", "node"))) {
#       ## psi
#       expect_lt(error(sort(newAlgo$fittedSampling$parameters), sort(refAlgo$psi_chap)), tol_ref)
#       expect_lte(
#         error(sort(newAlgo$fittedSampling$parameters), sort(refAlgo$psi_true)),
#         error(sort(refAlgo$psi_chap), sort(refAlgo$psi_true)),
#         tol_ref
#       )
#     }
#   }
# })
