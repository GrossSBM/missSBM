context("test consistency missSBM top-level function")

library(aricode)
source("utils_test.R")

referenceResults <- readRDS(system.file("extdata", "referenceResults.rds", package = "missSBM"))

test_that("check consistency against Tim's code for dyad, node, double standard and block sampling", {

  tol_ref   <- 1e-3
  tol_truth <- 5e-2
  tol_ARI   <- .9
  truth   <- referenceResults$true_sbm

  for (sampling in c("dyad", "node", "double-standard", "block-node")) {

    refAlgo <- referenceResults[[sampling]]

    missSBM_out <- missSBM::missSBM(
      adjacencyMatrix = refAlgo$Y,
      vBlocks = truth$Q,
      sampling = sampling,
      trace = FALSE
    )
    newAlgo <- missSBM_out[[1]]

    ## alpha
    expect_lt(error(newAlgo$fittedSBM$mixtureParam, refAlgo$alpha, sort = TRUE), tol_ref)
    expect_lt(error(newAlgo$fittedSBM$mixtureParam, truth$alpha, sort = TRUE), tol_truth)
    expect_lt(error(refAlgo$alpha, truth$alpha, sort = TRUE), tol_truth)

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
      expect_lt(error(newAlgo$fittedSampling$parameters, refAlgo$psi_chap, sort = TRUE), tol_ref)
      expect_lt(error(newAlgo$fittedSampling$parameters, refAlgo$psi_vrai, sort = TRUE), tol_truth)
      expect_lt(error(refAlgo$psi_chap, refAlgo$psi_vrai, sort = TRUE), tol_truth)
    }
  }
})

test_that("check consistency against Tim's code for dyad and node sampling with covariates", {

  truth   <- referenceResults$true_sbm_wCov
  tol_ref   <- 1e-2
  tol_truth <- 1e-2
  tol_ARI   <- .7

  for (sampling in c("covDyad", "covNode")) {

    refAlgo <- referenceResults[[sampling]]

    missSBM_out <- missSBM::missSBM(
      adjacencyMatrix = refAlgo$Y,
      vBlocks = truth$Q,
      sampling = ifelse(sampling == "covDyad", "dyad", "node"),
      trace = FALSE,
      covarMatrix = truth$covMatrix,
    )
    newAlgo <- missSBM_out[[1]]

    ## alpha
    expect_lt(error(newAlgo$fittedSBM$mixtureParam, refAlgo$alpha, sort = TRUE), tol_ref)
    expect_lt(error(newAlgo$fittedSBM$mixtureParam, truth$alpha, sort = TRUE), tol_truth)
    ## refAlgo is poor for alpha and clsutering, but good for estimating gamma !!
    ## expect_lt(error(refAlgo$alpha, truth$alpha, sort = TRUE), tol_truth)

    ## PI
    # expect_lt(error(newAlgo$fittedSBM$connectProb, refAlgo$pi), tol_ref)
    # expect_lt(error(missSBM:::logistic(newAlgo$fittedSBM$connectParam), missSBM:::logistic(truth$gamma)), tol_truth)
    # expect_lt(error(missSBM:::logistic(refAlgo$gamma_covDyad), missSBM:::logistic(truth$gamma)), tol_truth)

    ## clustering
##    expect_gt(ARI(newAlgo$fittedSBM$memberships, refAlgo$cl), tol_ARI)
    expect_gt(ARI(newAlgo$fittedSBM$memberships, truth$cl), tol_ARI)
   ##  expect_gt(ARI(refAlgo$cl, truth$cl), tol_ARI)

  }
})
