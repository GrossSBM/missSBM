context("Reordering of the blocks: the parameters of the sampling follow")

## run_VEM() permutes the blocks of the SBM once converged (SimpleSBM_fit$reorder()). The
## parameters of 'block-node' and 'block-dyad' are indexed by the blocks: they must follow,
## otherwise fittedSampling$parameters and log_lambda (used by polish()) refer to the former
## labels.

set.seed(3)
n <- 300; Q <- 3
sbm_sim <- sbm::sampleSimpleSBM(n, c(.5, .3, .2), list(mean = matrix(.05, Q, Q) + diag(c(.1, .3, .5))))
ctrl <- list(threshold = 1e-4, maxIter = 50, fixPointIter = 3, trace = FALSE)

test_that("block-node: rho is aligned with the blocks of the SBM after reordering", {
  A   <- observeNetwork(sbm_sim$networkData, "block-node", c(.2, .5, .9), clusters = sbm_sim$memberships)
  pon <- missSBM:::partlyObservedNetwork$new(A)
  ## the three labelings of the same partition: reorder() relabels at least one of them
  for (perm in list(c(1, 2, 3), c(2, 1, 3), c(3, 2, 1))) {
    fit <- missSBM_fit$new(pon, "block-node", perm[sbm_sim$memberships])
    fit$doVEM(ctrl)
    cl <- fit$fittedSBM$memberships
    ## near-hard memberships: rho_q is the proportion of observed nodes in block q
    expect_equal(fit$fittedSampling$parameters, as.numeric(tapply(pon$observedNodes, cl, mean)),
                 tolerance = 1e-3)
    ## log_lambda of an observed node is log(rho), in the order of the blocks of the SBM
    i <- which(pon$observedNodes)[1]
    expect_equal(as.numeric(fit$fittedSampling$log_lambda[i, ]), log(fit$fittedSampling$parameters))
  }
})

test_that("block-dyad: psi is aligned with the blocks of the SBM after reordering", {
  psi <- matrix(.3, Q, Q); diag(psi) <- c(.5, .7, .9)
  A   <- observeNetwork(sbm_sim$networkData, "block-dyad", psi, clusters = sbm_sim$memberships)
  pon <- missSBM:::partlyObservedNetwork$new(A)
  R   <- 1 * !is.na(A); diag(R) <- 0
  for (perm in list(c(1, 2, 3), c(2, 1, 3), c(3, 2, 1))) {
    fit <- missSBM_fit$new(pon, "block-dyad", perm[sbm_sim$memberships])
    fit$doVEM(ctrl)
    Z <- fit$fittedSBM$indMemberships
    ## near-hard memberships: psi_ql is the proportion of observed dyads between blocks q and l
    empirical <- (t(Z) %*% R %*% Z) / (colSums(Z) %o% colSums(Z) - diag(colSums(Z)))
    expect_equal(unname(fit$fittedSampling$parameters), unname(empirical), tolerance = 2e-2)
  }
})
