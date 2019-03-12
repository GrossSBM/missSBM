rm(list = ls())
library(missSBM)
source("tests/testthat/utils_test.R")

referenceResults <- readRDS(system.file("extdata", "referenceResults.rds", package = "missSBM"))

truth   <- referenceResults$true_sbm_cov
tol_ref   <- 1e-2
tol_truth <- 1e-2
tol_ARI   <- .7

refAlgo <- referenceResults[["dyad-covariates"]]

control <- list(threshold = 1e-4, maxIter = 100, fixPointIter = 3, trace = TRUE)

sampledNet <- sampleNetwork(refAlgo$sampledNet, "dyad", truth$covarParam, covarMatrix = refAlgo$covarMatrix)
newAlgo <- missSBM:::missingSBM_fit$new(sampledNet, truth$nBlocks, "dyad", covarMatrix = refAlgo$covarMatrix, clusterInit = "spectral")
missSBM_out <- newAlgo$doVEM(control)

## connectivity parameters (pi)
err_new <- error(newAlgo$fittedSBM$connectParam, truth$connectParam, sort = TRUE)
err_old <- error(refAlgo$connectParam          , truth$connectParam, sort = TRUE)
err_gap <- error(newAlgo$fittedSBM$connectParam, refAlgo$connectParam, sort = TRUE)
# if (err_new < err_old) {
#   expect_lt(err_new, tol_truth)
#   cat(" new better on connectivity")
# } else {
#   expect_lt(err_old, tol_ref)
#   expect_lt(err_old, tol_ref)
# }

## mixture parameters (alpha)
err_new <- error(newAlgo$fittedSBM$mixtureParam, truth$mixtureParam, sort = TRUE)
err_old <- error(refAlgo$mixtureParam          , truth$mixtureParam, sort = TRUE)
gap_old <- error(newAlgo$fittedSBM$mixtureParam, refAlgo$mixtureParam, sort = TRUE)
# if (err_new < err_old) {
#   expect_lt(err_new, tol_truth)
#   cat(" new better on mixture")
# } else {
#   expect_lt(gap_old, tol_ref)
#   expect_lt(gap_old, tol_ref)
# }


## clustering
ARI_new <- ARI(newAlgo$fittedSBM$memberships, truth$memberships)
ARI_old <- ARI(refAlgo$memberships, truth$memberships)
ARI_gap <- ARI(refAlgo$memberships, newAlgo$fittedSBM$memberships)
# expect_gt(ARI_new, tol_ARI)
# expect_gt(ARI_old, tol_ARI)
# expect_gt(ARI_gap, tol_ARI)
