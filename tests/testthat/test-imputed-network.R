context("imputed network: the imputations at the final parameters")

## The designs that integrate the missing dyads out given the blocks compute no imputation along
## the VEM: imputedNetwork computes them on demand, at the final parameters.

set.seed(1)
Q <- 3
pi <- matrix(.05, Q, Q); diag(pi) <- c(.6, .3, .15)
sim <- sbm::sampleSimpleSBM(100, c(.25, .35, .4), list(mean = pi))
par <- list(dyad = .5, node = .5, "block-node" = c(.8, .5, .2),
            "block-dyad" = {p <- matrix(.3, Q, Q); diag(p) <- c(.7, .5, .3); p})

test_that("imputedNetwork: logistic(tau logit(pi) tau') at the missing dyads, the data elsewhere", {
  for (design in names(par)) {
    A <- observeNetwork(sim$networkData, design, par[[design]], clusters = sim$memberships)
    fit <- estimateMissSBM(A, Q, design, control = missSBM_param(trace = FALSE, iterates = 0))$models[[1]]
    tau <- fit$fittedSBM$probMemberships
    expected <- plogis(tau %*% qlogis(fit$fittedSBM$connectParam$mean) %*% t(tau))
    Y <- as.matrix(fit$imputedNetwork)
    miss <- is.na(A); diag(miss) <- FALSE
    expect_equal(Y[miss], expected[miss])
    obs <- !is.na(A); diag(obs) <- FALSE
    expect_equal(Y[obs], A[obs])
    expect_true(isSymmetric(Y))
  }
})
