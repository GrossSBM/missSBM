context("test consistency missSBM top-level function")

library(aricode)
error <- function(beta1, beta2, sort = FALSE) {
  if (sort)
    err <- sum((sort(beta1) - sort(beta2))^2)/length(beta2)
  else
    err <- sum((beta1 - beta2)^2)/length(beta2)
  err
}

referenceResults <- readRDS(system.file("extdata", "referenceResults.rds", package = "missSBM"))

test_that("check consistency against Tim's code for dyad, node, double standard and block sampling", {

  tol_ref   <- 1e-3
  tol_truth <- 1e-3
  tol_ARI   <- .8
  truth   <- referenceResults$true_sbm

  cat("\nsampling: ")
  for (sampling in c("dyad", "node", "double-standard", "block-node")) {

    refAlgo <- referenceResults[[sampling]]
    cat(sampling)

    missSBM_out <- estimateMissSBM(
      adjacencyMatrix = refAlgo$sampledNet,
      vBlocks = truth$nBlocks,
      sampling = sampling,
      control = list(trace = 0)
    )
    newAlgo <- missSBM_out$bestModel

    ## mixture parameters (pi)
    err_new <- error(newAlgo$fittedSBM$blockProp, truth$mixtureParam, sort = TRUE)
    err_old <- error(refAlgo$mixtureParam       , truth$mixtureParam, sort = TRUE)
    gap_old <- error(newAlgo$fittedSBM$blockProp, refAlgo$mixtureParam, sort = TRUE)
    if (err_new < err_old) {
      expect_lt(err_new, tol_truth)
      cat(" new better on mixture")
    } else {
      expect_lt(gap_old, tol_ref)
      expect_lt(gap_old, tol_ref)
    }

    ## connectivity parameters (theta)
    err_new <- error(newAlgo$fittedSBM$connectParam$mean, truth$connectParam, sort = TRUE)
    err_old <- error(refAlgo$connectParam               , truth$connectParam, sort = TRUE)
    err_gap <- error(newAlgo$fittedSBM$connectParam$mean, refAlgo$connectParam, sort = TRUE)
    if (err_new < err_old) {
      expect_lt(err_new, tol_truth)
      cat(" new better on connectivity")
    } else {
      expect_lt(err_old, tol_ref)
      expect_lt(err_gap, tol_ref)
    }

    ## clustering
    ARI_new <- ARI(newAlgo$fittedSBM$memberships, truth$memberships)
    ARI_old <- ARI(refAlgo$memberships, truth$memberships)
    ARI_gap <- ARI(refAlgo$memberships, newAlgo$fittedSBM$memberships)
    expect_gt(ARI_new, tol_ARI)
    expect_gt(ARI_old, tol_ARI)
    expect_gt(ARI_gap, tol_ARI)

    if (!(sampling %in% c("dyad", "node"))) {
      ## psi
      err_new  <- error(newAlgo$fittedSampling$parameters, refAlgo$true_samplingParam, sort = TRUE)
      err_old  <- error(refAlgo$samplingParam            , refAlgo$true_samplingParam, sort = TRUE)
      err_gap  <- error(newAlgo$fittedSampling$parameters, refAlgo$samplingParam, sort = TRUE)
      if (err_new < err_old) {
        expect_lt(err_new, tol_truth*10)
        cat(" new better on sampling parameters")
      } else {
        expect_lt(err_old, tol_ref)
        expect_lt(err_gap, tol_ref)
      }
    }
    cat("\n")
  }
})

test_that("check consistency against Tim's code for dyad and node sampling with covariates", {

  truth   <- referenceResults$true_sbm_cov
  tol_ref   <- 1e-3
  tol_truth <- 1e-3
  tol_ARI   <- .7

  covarMatrix <- referenceResults$`dyad-covariates`$covarMatrix
  covarArray  <- missSBM:::getCovarArray(covarMatrix, missSBM:::l1_similarity)

  referenceResults$`dyad-covariates`$covariates <- covarArray
  referenceResults$`node-covariates`$covariates <- covarMatrix
  for (sampling in c("dyad-covariates", "node-covariates")) {

    refAlgo <- referenceResults[[sampling]]

    missSBM_out <- estimateMissSBM(
      adjacencyMatrix = refAlgo$sampledNet,
      vBlocks     = truth$nBlocks,
      sampling    = ifelse(sampling == "dyad-covariates", "covar-dyad", "covar-node"),
      covariates  = refAlgo$covariates,
      control     = list(trace = 0)
    )
    newAlgo <- missSBM_out$bestModel

    ## mixture parameters (pi)
    err_new <- error(newAlgo$fittedSBM$blockProp, truth$mixtureParam, sort = TRUE)
    err_old <- error(refAlgo$mixtureParam       , truth$mixtureParam, sort = TRUE)
    err_gap <- error(newAlgo$fittedSBM$blockProp, refAlgo$mixtureParam, sort = TRUE)
    if (err_new < err_old) {
      expect_lt(err_new, tol_truth)
      cat(" new better on mixture")
    } else {
      expect_lt(err_old, tol_ref)
      expect_lt(err_gap, tol_ref)
    }

    ## connectivity parameters (theta)
    err_new <- error(newAlgo$fittedSBM$connectParam$mean, missSBM:::.logistic(truth$connectParam), sort = TRUE)
    err_old <- error(missSBM:::.logistic(refAlgo$connectParam)    , missSBM:::.logistic(truth$connectParam), sort = TRUE)
    err_gap <- error(newAlgo$fittedSBM$connectParam$mean, missSBM:::.logistic(refAlgo$connectParam), sort = TRUE)
    if (err_new < err_old) {
      expect_lt(err_new, tol_truth*3)
      cat(" new better on connectivity")
    } else {
      expect_lt(err_new, 10*tol_ref)
      expect_lt(err_old, 10*tol_ref)
    }

    ## clustering
    ARI_new <- ARI(newAlgo$fittedSBM$memberships, truth$memberships)
    ARI_old <- ARI(refAlgo$memberships, truth$memberships)
    ARI_gap <- ARI(refAlgo$memberships, newAlgo$fittedSBM$memberships)
    expect_gt(ARI_new, tol_ARI)
    expect_gt(ARI_old, tol_ARI)
    expect_gt(ARI_gap, tol_ARI)

  }
})
