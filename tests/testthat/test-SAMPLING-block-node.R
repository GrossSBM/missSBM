context("block-node sampling: the SBM part is fitted on the observed dyads")

## With a node-centered sampling, the dyads between two unobserved nodes are never observed.
## Fitting the SBM on their imputations lets them feed back on themselves: the connectivity of a
## block of unobserved nodes could drift to 0 and the clustering was lost (see NEWS).

for (seed in c(3, 12)) {
  test_that(paste("block-node sampling: the clustering and the connectivity are recovered, seed", seed), {
    set.seed(seed)
    n <- 200; Q <- 3
    pi <- matrix(.05, Q, Q); diag(pi) <- c(.6, .3, .15); pi[1, 2] <- pi[2, 1] <- .15
    sbm <- sbm::sampleSimpleSBM(n, c(.2, .3, .5), list(mean = pi))
    A_obs <- missSBM::observeNetwork(sbm$networkData, "block-node", c(.8, .5, .2), clusters = sbm$memberships)

    fit <- estimateMissSBM(A_obs, vBlocks = Q, sampling = "block-node",
                           control = missSBM_param(trace = FALSE, iterates = 0, threshold = 1e-5, maxIter = 300))$models[[1]]

    expect_gt(aricode::ARI(fit$fittedSBM$memberships, sbm$memberships), 0.85)
    expect_gt(min(diag(fit$fittedSBM$connectParam$mean)), 0.1)
    ## the objective increases along the optimization
    expect_true(all(diff(fit$monitoring$elbo) >= -1e-8))
  })
}
