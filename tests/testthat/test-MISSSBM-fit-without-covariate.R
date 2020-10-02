context("test missSBM-fit without covariate")

library(aricode)
source("utils_test.R")

set.seed(1890718)
### A SBM model : ###
N <- 200
Q <- 3
pi <- rep(1, Q)/Q           # block proportion
theta <- diag(.45, Q, Q) + .05 # connectivity matrix
directed <- FALSE         # if the network is directed or not

### Draw a SBM model
sbm <- missSBM::simulate(N, pi, theta, directed) # simulation of ad Bernoulli non-directed SBM

samplings <- list(
  list(name = "dyad", psi = 0.5, class = "dyadSampling_fit"),
  list(name = "node", psi = 0.5, class = "nodeSampling_fit"),
  list(name = "double-standard", psi =  c(.3, .6), class = "doubleStandardSampling_fit"),
  list(name = "block-node", psi = c(.3, .5, .7), class = "blockSampling_fit")
)

## control parameter for the VEM
control <- list(threshold = 1e-4, maxIter = 200, fixPointIter = 5, trace = FALSE)

test_that("missSBM-fit works and is consistent for all samplings", {

  ## Consistency
  tol_truth <- 1e-2
  tol_ARI   <- .8

##  cat("Tested sampling:")
  for (sampling in samplings) {
##    cat("\n -", sampling$name)

    ## sampled the network
    adjMatrix  <- missSBM::sample(sbm$netMatrix, sampling$name, sampling$psi, sbm$memberships)
    sampledNet <- missSBM:::sampledNetwork$new(adjMatrix)

    ## Perform inference
    missSBM <- missSBM:::missSBM_fit$new(sampledNet, Q, sampling$name, "hierarchical", FALSE)
    out <- missSBM$doVEM(control)

    ## Sanity check
    expect_is(missSBM, "missSBM_fit")
    expect_is(missSBM$fittedSBM, "SBM_fit_nocovariate")
    expect_is(missSBM$fittedSampling, sampling$class)
    expect_is(missSBM$sampledNetwork, "sampledNetwork")
    expect_equal(out, missSBM$monitoring)

    ## Optimization success
    expect_gt(diff(range(out$objective, na.rm = TRUE)), 0)

    ## SBM: parameters estimation
    expect_lt(error(missSBM$fittedSBM$connectParam, sbm$connectParam), tol_truth)
    expect_lt(error(missSBM$fittedSBM$blockProp, sbm$blockProp, sort = TRUE), tol_truth)

    ## sampling design: parameters estimation
    expect_lt(error(missSBM$fittedSampling$parameters, sampling$psi, sort = TRUE), tol_truth)

    ## clustering
    expect_gt(ARI(missSBM$fittedSBM$memberships, sbm$memberships), tol_ARI)
  }

})

# test_that("miss SBM with degree sampling works", {
#
#   psi <- c(-5, .1)
#   sampledNet <- sampleNetwork(A, "degree", psi)
#   ## Perform inference
#   missSBM <- missSBM:::missSBM_fit$new(sampledNet, Q, "degree", "hierarchical)
#   out <- missSBM$doVEM(control)
#
#   ## Sanity check
#   expect_is(missSBM, "missSBM_fit")
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
