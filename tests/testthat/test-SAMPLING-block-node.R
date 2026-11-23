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

## The SBM part of block-node and block-dyad is fitted on the observed dyads, but their ICL must
## include the imputed dyads as for the other designs, so that the designs can be compared: with a
## single block, block-node is node sampling and block-dyad is dyad sampling, and the ICLs coincide.
test_that("block-node and block-dyad: the ICL is comparable with node and dyad sampling", {
  set.seed(5)
  n <- 120
  pi <- matrix(.05, 2, 2); diag(pi) <- c(.4, .2)
  sbm <- sbm::sampleSimpleSBM(n, c(.4, .6), list(mean = pi))
  ctrl <- missSBM_param(trace = FALSE, iterates = 0, threshold = 1e-8, maxIter = 500)

  A_node <- missSBM::observeNetwork(sbm$networkData, "node", .6)
  ICL_node       <- estimateMissSBM(A_node, 1, "node",       control = ctrl)$ICL
  ICL_block_node <- estimateMissSBM(A_node, 1, "block-node", control = ctrl)$ICL
  expect_equal(ICL_block_node, ICL_node, tolerance = 1e-4)

  A_dyad <- missSBM::observeNetwork(sbm$networkData, "dyad", .6)
  ICL_dyad       <- estimateMissSBM(A_dyad, 1, "dyad",       control = ctrl)$ICL
  ICL_block_dyad <- estimateMissSBM(A_dyad, 1, "block-dyad", control = ctrl)$ICL
  expect_equal(ICL_block_dyad, ICL_dyad, tolerance = 1e-4)
})
