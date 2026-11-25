context("double-standard sampling: the missing dyads are integrated out given the blocks")

## A dyad is observed with probability p0 if it is 0, p1 if it is 1. Given the blocks, a dyad
## between blocks q and l is missing with probability m_ql = pi_ql (1 - p1) + (1 - pi_ql) (1 - p0):
## the missing dyads are integrated out, and the posterior probability that a missing dyad is an
## edge is nu_ql = pi_ql (1 - p1) / m_ql.

sim_ds <- function(seed, n = 200, psi = c(.2, .8), directed = FALSE) {
  set.seed(seed)
  Q  <- 3
  pi <- matrix(.05, Q, Q); diag(pi) <- c(.6, .3, .15); pi[1, 2] <- pi[2, 1] <- .15
  sbm <- sbm::sampleSimpleSBM(n, c(.2, .3, .5), list(mean = pi), directed = directed)
  list(sbm = sbm, A = observeNetwork(sbm$networkData, "double-standard", psi))
}
ctrl_vem <- list(threshold = 1e-4, maxIter = 100, fixPointIter = 3, trace = FALSE)

## the variational bound, computed dyad by dyad: sum_i sum_q tau_iq log(alpha_q) plus, for each
## dyad, sum_{q,l} tau_iq tau_jl log p(Y_ij^o, R_ij | Z_i = q, Z_j = l)
brute_force_bound <- function(fit, A, directed) {
  tau   <- fit$fittedSBM$probMemberships
  theta <- fit$fittedSBM$connectParam$mean
  alpha <- fit$fittedSBM$blockProp
  p0 <- fit$fittedSampling$parameters[1]; p1 <- fit$fittedSampling$parameters[2]
  res <- sum(tau %*% log(alpha))
  n <- nrow(A)
  for (i in 1:n) for (j in 1:n) {
    if (i == j || (!directed && j < i)) next
    logp <- if (is.na(A[i, j])) log(theta * (1 - p1) + (1 - theta) * (1 - p0))
      else if (A[i, j] == 1) log(theta) + log(p1) else log(1 - theta) + log(p0)
    res <- res + sum(outer(tau[i, ], tau[j, ]) * logp)
  }
  res
}

test_that("the bound integrates the missing dyads out exactly (undirected and directed)", {
  for (directed in c(FALSE, TRUE)) {
    sim <- sim_ds(1, n = 30, directed = directed)
    fit <- missSBM_fit$new(missSBM:::partlyObservedNetwork$new(sim$A), "double-standard", sim$sbm$memberships)
    fit$doVEM(list(threshold = 1e-4, maxIter = 3, fixPointIter = 1, trace = FALSE))
    ## soft memberships: the check covers the bilinear forms in tau
    expect_true(any(fit$fittedSBM$probMemberships > 1e-3 & fit$fittedSBM$probMemberships < 1 - 1e-3))
    expect_equal(fit$vExpec, brute_force_bound(fit, sim$A, directed), tolerance = 1e-10)
    expect_equal(fit$loglik, fit$vExpec + fit$fittedSBM$entropy)
  }
})

test_that("the clustering and the sampling parameters are recovered, the bound increases", {
  for (psi in list(c(.2, .8), c(.6, .3))) {
    sim <- sim_ds(2, psi = psi)
    fit <- estimateMissSBM(sim$A, 3, "double-standard",
                           control = missSBM_param(trace = FALSE, iterates = 0))$models[[1]]
    expect_gt(aricode::ARI(fit$fittedSBM$memberships, sim$sbm$memberships), 0.9)
    expect_lt(max(abs(fit$fittedSampling$parameters - psi)), 0.1)
    expect_true(all(diff(fit$monitoring$elbo) >= -1e-8))
  }
})

test_that("the joint M-step on (theta, p0, p1) is iterated to convergence", {
  sim <- sim_ds(3)
  fit <- missSBM_fit$new(missSBM:::partlyObservedNetwork$new(sim$A), "double-standard", sim$sbm$memberships)
  fit$doVEM(ctrl_vem)
  ## a fixed point: one more M-step from the current parameters leaves them unchanged
  theta <- fit$fittedSBM$connectParam$mean
  psi   <- fit$fittedSampling$parameters
  psi_new <- fit$fittedSBM$update_parameters(psi)
  expect_equal(psi_new, psi, tolerance = 1e-6)
  expect_equal(fit$fittedSBM$connectParam$mean, theta, tolerance = 1e-6)
})

test_that("the imputed network averages the posterior given the blocks at the missing dyads", {
  sim <- sim_ds(4, n = 60)
  fit <- missSBM_fit$new(missSBM:::partlyObservedNetwork$new(sim$A), "double-standard", sim$sbm$memberships)
  fit$doVEM(ctrl_vem)
  tau <- fit$fittedSBM$probMemberships
  expected <- tau %*% fit$fittedSampling$blockImputation %*% t(tau)
  Y <- as.matrix(fit$imputedNetwork)
  miss <- is.na(sim$A); diag(miss) <- FALSE
  expect_equal(Y[miss], expected[miss])
  obs <- !is.na(sim$A); diag(obs) <- FALSE
  expect_equal(Y[obs], sim$A[obs])
  expect_true(isSymmetric(Y))
})

test_that("missing at random: p0 and p1 are estimated close to each other", {
  sim <- sim_ds(5, psi = c(.5, .5))
  fit <- missSBM_fit$new(missSBM:::partlyObservedNetwork$new(sim$A), "double-standard", sim$sbm$memberships)
  fit$doVEM(ctrl_vem)
  expect_lt(abs(diff(fit$fittedSampling$parameters)), 0.05)
})

test_that("with a single block: not identifiable, fitted as dyad sampling with a warning", {
  sim <- sim_ds(6, n = 100)
  ctrl <- missSBM_param(trace = FALSE, iterates = 0)
  expect_warning(
    col <- estimateMissSBM(sim$A, 1:2, "double-standard", control = ctrl),
    "not identifiable"
  )
  expect_silent(estimateMissSBM(sim$A, 2, "double-standard", control = ctrl))
  fit <- col$models[[1]]
  rate <- mean(!is.na(sim$A[upper.tri(sim$A)]))
  expect_equal(fit$fittedSampling$parameters, c(rate, rate))
  expect_equal(fit$fittedSampling$df, 1)
  expect_equal(col$models[[2]]$fittedSampling$df, 2)
  ## the same ICL as dyad sampling
  ICL_ds   <- suppressWarnings(estimateMissSBM(sim$A, 1, "double-standard", control = ctrl))$ICL
  ICL_dyad <- estimateMissSBM(sim$A, 1, "dyad", control = ctrl)$ICL
  expect_equal(ICL_ds, ICL_dyad, tolerance = 1e-8)
})
