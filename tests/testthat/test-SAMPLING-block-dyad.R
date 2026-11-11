context("block-dyad sampling: the SBM part is fitted on the observed dyads")

## With few observed dyads, fitting the SBM on the imputations of the missing ones lets them feed
## back on themselves: the connectivity of a block could drift to 0 and the clustering was lost
## (see NEWS). The missing dyads are ignorable given the blocks.

for (seed in c(4, 12)) {
  test_that(paste("block-dyad sampling: the clustering and the connectivity are recovered, seed", seed), {
    set.seed(seed)
    n <- 200; Q <- 3
    pi <- matrix(.05, Q, Q); diag(pi) <- c(.6, .3, .15); pi[1, 2] <- pi[2, 1] <- .15
    psi <- matrix(.3, Q, Q); diag(psi) <- c(.9, .5, .1)
    sbm <- sbm::sampleSimpleSBM(n, c(.2, .3, .5), list(mean = pi))
    A_obs <- missSBM::observeNetwork(sbm$networkData, "block-dyad", psi, clusters = sbm$memberships)

    fit <- estimateMissSBM(A_obs, vBlocks = Q, sampling = "block-dyad",
                           control = missSBM_param(trace = FALSE, iterates = 0, threshold = 1e-5, maxIter = 300))$models[[1]]

    expect_gt(aricode::ARI(fit$fittedSBM$memberships, sbm$memberships), 0.85)
    expect_gt(min(diag(fit$fittedSBM$connectParam$mean)), 0.1)
    ## the objective increases along the optimization
    expect_true(all(diff(fit$monitoring$elbo) >= -1e-8))
  })
}

## the expectation of the sampling log-likelihood is E_q[log p(R | Z)], the expectation of the log
## (consistent with the update of psi), not log(Z psi Z^T), the log of the expectation (Jensen)
for (directed in c(FALSE, TRUE)) {
  test_that(paste("block-dyad sampling: vExpec is the expectation of the log,", if (directed) "directed" else "undirected"), {
    set.seed(7)
    n <- 80; Q <- 3
    pi <- matrix(.05, Q, Q); diag(pi) <- c(.6, .3, .15)
    psi <- matrix(.3, Q, Q); diag(psi) <- c(.8, .5, .2)
    sbm <- sbm::sampleSimpleSBM(n, rep(1/3, 3), list(mean = pi), directed = directed)
    A_obs <- missSBM::observeNetwork(sbm$networkData, "block-dyad", psi, clusters = sbm$memberships)
    pon <- missSBM:::partlyObservedNetwork$new(A_obs)
    Z <- matrix(runif(n * Q), n, Q); Z <- Z / rowSums(Z) # soft memberships, where log E and E log differ
    fit <- missSBM:::blockDyadSampling_fit$new(pon, Z)
    R <- as.matrix(pon$samplingMatrix) != 0
    S <- as.matrix(pon$samplingMatrixBar) != 0
    psi_hat <- fit$parameters
    E_log <- sum((Z %*% log(psi_hat) %*% t(Z))[R]) + sum((Z %*% log(1 - psi_hat) %*% t(Z))[S])
    log_E <- sum(log(Z %*% psi_hat %*% t(Z))[R]) + sum(log(1 - Z %*% psi_hat %*% t(Z))[S])
    expect_equal(fit$vExpec, E_log, tolerance = 1e-10)
    expect_lt(fit$vExpec, log_E)
  })
}
