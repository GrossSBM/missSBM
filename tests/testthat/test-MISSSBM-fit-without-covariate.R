context("test missSBM-fit without covariate")

library(aricode)
error <- function(beta1, beta2, sort = FALSE) {
  if (sort)
    err <- sum((sort(beta1) - sort(beta2))^2)/length(beta2)
  else
    err <- sum((beta1 - beta2)^2)/length(beta2)
  err
}

set.seed(1890718)
### A SBM model : ###
N <- 200
Q <- 3
pi <- rep(1, Q)/Q           # block proportion
theta <- list(mean = diag(.45, Q, Q) + .05) # connectivity matrix
directed <- FALSE         # if the network is directed or not

### Draw a SBM model
sbm <- sbm::sampleSimpleSBM(N, pi, theta) # simulation of ad Bernoulli non-directed SBM

samplings <- list(
  list(name = "dyad", psi = 0.5, class = "dyadSampling_fit"),
  list(name = "node", psi = 0.5, class = "nodeSampling_fit"),
  list(name = "double-standard", psi =  c(.3, .6), class = "doubleStandardSampling_fit")#,
#  list(name = "block-node", psi = c(.3, .5, .7), class = "blockSampling_fit")
)

## control parameter for the VEM
control <- list(threshold = 1e-3, maxIter = 200, fixPointIter = 5, trace = FALSE)

test_that("missSBM-fit works and is consistent for all samplings", {

  ## Consistency
  tol_truth <- 5e-3
  tol_ARI   <- .9

  cat("Tested sampling:")
  for (sampling in samplings) {
    cat("\n -", sampling$name)

    ## sampled the network
    adjMatrix  <- missSBM::observeNetwork(sbm$networkData, sampling$name, sampling$psi, sbm$memberships)
    myNet <- missSBM:::partlyObservedNetwork$new(adjMatrix)
    cl <- myNet$clustering(1:Q)

    ## Perform inference
    missSBM <- missSBM:::missSBM_fit$new(myNet, sampling$name, cl[[Q]], FALSE)
    out <- missSBM$doVEM(control)

    ## Sanity check
    expect_true(inherits(missSBM, "missSBM_fit"))
    expect_true(inherits(missSBM$fittedSBM, "SimpleSBM_fit"))
    expect_is(missSBM$fittedSampling, sampling$class)
    expect_equal(out, missSBM$monitoring)

    ## Optimization success
    expect_gte(diff(range(out$objective, na.rm = TRUE)), 0)

    ## SBM: parameters estimation
    expect_lt(error(missSBM$fittedSBM$connectParam$mean, sbm$connectParam$mean), tol_truth)
    expect_lt(error(missSBM$fittedSBM$blockProp, sbm$blockProp, sort = TRUE), tol_truth)

    ## sampling design: parameters estimation
    expect_lt(error(missSBM$fittedSampling$parameters, sampling$psi, sort = TRUE), tol_truth)

    ## clustering
    expect_gt(ARI(missSBM$fittedSBM$memberships, sbm$memberships), tol_ARI)
  }

})

