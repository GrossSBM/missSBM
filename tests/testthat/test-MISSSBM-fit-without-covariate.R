context("test-sbm-with-missing-values")

library(aricode)

set.seed(1890718)


### A SBM model : ###
N <- 400
Q <- 5
alpha <- rep(1,Q)/Q       # mixture parameter
pi <- diag(.45,Q) + .05   # connectivity matrix
directed <- FALSE         # if the network is directed or not

### Draw a SBM model
sbm <- simulateSBM(N, alpha, pi, directed) # simulation of ad Bernoulli non-directed SBM

truth <- list(
  "sbm"  = sbm,
  "dyad"            = list(psi = 0.3, class = "dyadSampling_fit"),
  "node"            = list(psi = 0.2, class = "nodeSampling_fit"),
  "double-standard" = list(psi =  c(.3, .6), class = "doubleStandardSampling_fit"),
  "block-node"      = list(psi = c(.1, .3, .2, .5, .7), class = "blockSampling_fit")
)

error <- function(beta1, beta2, sort = FALSE) {
  if (sort)
    err <- sum((sort(beta1) - sort(beta2))^2)/Q
  else
    err <- sum((beta1 - beta2)^2)/Q
  err
}

test_that("missSBM-fit works and is consistent for all samplings", {

  ## Consistency
  tol_truth <- 1e-2
  tol_ARI   <- .8

  cat("Tested sampling:")
  for (sampling in c("dyad", "node", "double-standard", "block-node")) {
    cat("\n -", sampling)
    refSampling <- truth[[sampling]]

    ## sampled the network
    sampledNet <- sampleNetwork(A, sampling, refSampling$psi, truth$sbm$memberships)

    ## Perform inference
    missSBM <- missSBM:::missingSBM_fit$new(sampledNet, Q, sampling)
    out <- missSBM$doVEM(control)

    ## Sanity check
    expect_is(missSBM, "missingSBM_fit")
    expect_is(missSBM$fittedSBM, "SBM_fit_nocovariate")
    expect_is(missSBM$fittedSampling, refSampling$class)
    expect_is(missSBM$sampledNetwork, "sampledNetwork")
    expect_equal(out, missSBM$monitoring)

    ## Optimization success
    expect_gt(diff(range(out$objective, na.rm = TRUE)), 0)

    ## SBM: parameters estimation
    expect_lt(error(missSBM$fittedSBM$connectParam, truth$sbm$connectParam), tol_truth)
    expect_lt(error(missSBM$fittedSBM$mixtureParam, truth$sbm$mixtureParam, sort = TRUE), tol_truth)

    ## sampling design: parameters estimation
    expect_lt(error(missSBM$fittedSampling$parameters, refSampling$psi, sort = TRUE), tol_truth)

    ## clustering
    expect_gt(ARI(missSBM$fittedSBM$memberships, mySBM$memberships), tol_ARI)
  }

})

# test_that("miss SBM with degree sampling works", {
#
#   psi <- c(-5, .1)
#   sampledNet <- sampleNetwork(A, "degree", psi)
#   ## Perform inference
#   missSBM <- missSBM:::missingSBM_fit$new(sampledNet, Q, "degree")
#   out <- missSBM$doVEM(control)
#
#   ## Sanity check
#   expect_is(missSBM, "missingSBM_fit")
#   expect_is(missSBM$fittedSBM, "SBM_fit_nocovariate")
#   expect_is(missSBM$fittedSampling, "degreeSampling_fit")
#   expect_is(missSBM$sampledNetwork, "sampledNetwork")
#
#   ## Consistency
#   tol <- 1e-2
#   ## Optimization success
#   expect_gte(diff(range(out$objective, na.rm = TRUE)), 0)
#   ## SBM: parameters estimation
#   ## FIXME: expect_lt(sum((missSBM$fittedSBM$connectParam - mySBM$connectParam)^2)/(Q*(Q + 1)/2), tol)
#   ## sampling design: parameters estimation
#   ## FIXME: this does work!!! expect_lt(sum((sort(missSBM$fittedSampling$parameters) - sort(psi))^2/Q), tol)
#   ## clustering
#   tol <- .9
#   expect_gt(ARI(missSBM$fittedSBM$memberships, mySBM$memberships), tol)
#
# })
#
#
#
