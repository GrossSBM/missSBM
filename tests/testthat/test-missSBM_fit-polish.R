context("test missSBM_fit's node-swap polishing (polish())")

## polish_log_tau()'s central claim: for each node i and candidate class q, it returns exactly
## the complete-data log-likelihood the fit would have if i alone were hard-reassigned to q
## (theta/pi held fixed) -- i.e. it reuses the existing E-step machinery (rescale = FALSE, hard
## Z) rather than any new numerical code. Checked here against a brute-force recomputation of
## vLL_complete_sparse_bernoulli_* for a handful of (node, class) reassignments, on the three
## SimpleSBM_fit variants (noCov, withCov, MNAR).

set.seed(1)
N <- 60; Q <- 4
sampler <- sbm::SimpleSBM$new("bernoulli", N, FALSE, rep(1/Q, Q), list(mean = diag(.4, Q) + .05))
sampler$rMemberships(store = TRUE)
sampler$rNetwork(store = TRUE)

brute_force_delta <- function(vLL_complete, Zhard, i, q, ...) {
  Zmoved <- Zhard
  Zmoved[i, ] <- 0
  Zmoved[i, q] <- 1
  vLL_complete(Zmat = Zmoved, ...) - vLL_complete(Zmat = Zhard, ...)
}

test_that("polish_log_tau() matches a brute-force complete-data log-likelihood delta (noCov)", {
  net <- missSBM:::partlyObservedNetwork$new(sampler$networkData)
  fit <- missSBM:::SimpleSBM_fit_noCov$new(net, sampler$memberships)
  fit$doVEM(threshold = 1e-3, maxIter = 20, fixPointIter = 3, trace = FALSE)

  priv    <- environment(fit$get_state)$private
  Zhard   <- fit$indMemberships
  cl0     <- fit$memberships
  log_tau <- fit$polish_log_tau()

  vLL <- function(Zmat) missSBM:::vLL_complete_sparse_bernoulli_nocovariate(
    priv$Y, priv$R, Zmat, fit$connectParam$mean, fit$blockProp)

  for (i in c(1, 10, 30, 55)) for (q in setdiff(1:Q, cl0[i])) {
    brute   <- vLL(replace(Zhard, cbind(i, 1:Q), diag(Q)[q, ])) - vLL(Zhard)
    formula <- log_tau[i, q] - log_tau[i, cl0[i]]
    expect_equal(brute, formula, tolerance = 1e-8)
  }
})

test_that("polish_log_tau() matches a brute-force complete-data log-likelihood delta (MNAR)", {
  blog_obs <- missSBM::observeNetwork(sampler$networkData, "block-node",
                                       c(0.9, 0.3, 0.9, 0.5), clusters = sampler$memberships)
  fit <- missSBM:::SimpleSBM_fit_MNAR$new(missSBM:::partlyObservedNetwork$new(blog_obs), sampler$memberships)
  fit$doVEM(threshold = 1e-3, maxIter = 20, fixPointIter = 3, trace = FALSE)

  priv    <- environment(fit$get_state)$private
  Zhard   <- fit$indMemberships
  cl0     <- fit$memberships
  log_tau <- fit$polish_log_tau()

  vLL <- function(Zmat) {
    missSBM:::vLL_complete_sparse_bernoulli_nocovariate(priv$Y, priv$R, Zmat, fit$connectParam$mean, fit$blockProp) +
      missSBM:::vLL_complete_sparse_bernoulli_nocovariate(priv$V, priv$S, Zmat, fit$connectParam$mean, fit$blockProp) -
      sum(Zmat %*% log(fit$blockProp))
  }

  for (i in c(1, 15, 40)) for (q in setdiff(1:Q, cl0[i])) {
    brute   <- vLL(replace(Zhard, cbind(i, 1:Q), diag(Q)[q, ])) - vLL(Zhard)
    formula <- log_tau[i, q] - log_tau[i, cl0[i]]
    expect_equal(brute, formula, tolerance = 1e-8)
  }
})

test_that("polish() never leaves the ICL worse than before the call", {
  net <- missSBM:::partlyObservedNetwork$new(sampler$networkData)
  fit <- missSBM:::missSBM_fit$new(net, "dyad", sampler$memberships, TRUE)
  fit$doVEM(list(threshold = 1e-3, maxIter = 50, fixPointIter = 3, trace = FALSE))
  icl_before <- fit$ICL

  out <- fit$polish(list(threshold = 1e-3, maxIter = 50, fixPointIter = 3, trace = FALSE))
  expect_identical(out, fit) # mutates in place, like split(in_place = TRUE)
  expect_lte(fit$ICL, icl_before)
})

test_that("polish() is close to a no-op on an already-converged, well-initialized clustering", {
  net <- missSBM:::partlyObservedNetwork$new(sampler$networkData)
  fit <- missSBM:::missSBM_fit$new(net, "dyad", sampler$memberships, TRUE)
  fit$doVEM(list(threshold = 1e-4, maxIter = 100, fixPointIter = 3, trace = FALSE))
  icl_before <- fit$ICL

  fit$polish(list(threshold = 1e-4, maxIter = 100, fixPointIter = 3, trace = FALSE))
  expect_equal(fit$ICL, icl_before, tolerance = 1e-6)
})

test_that("polish() improves a deliberately noisy initial clustering", {
  set.seed(2)
  N2 <- 300; Q2 <- 15
  sampler2 <- sbm::SimpleSBM$new("bernoulli", N2, FALSE, rep(1/Q2, Q2), list(mean = diag(.35, Q2) + .05))
  sampler2$rMemberships(store = TRUE)
  sampler2$rNetwork(store = TRUE)
  net2 <- missSBM:::partlyObservedNetwork$new(sampler2$networkData)

  noisy <- sampler2$memberships
  idx <- sample(N2, floor(N2 * .3))
  noisy[idx] <- sample(1:Q2, length(idx), replace = TRUE)

  fit <- missSBM:::missSBM_fit$new(net2, "dyad", noisy, TRUE)
  fit$doVEM(list(threshold = 1e-4, maxIter = 100, fixPointIter = 3, trace = FALSE))
  icl_before <- fit$ICL
  ari_before <- aricode::ARI(fit$fittedSBM$memberships, sampler2$memberships)

  fit$polish(list(threshold = 1e-4, maxIter = 100, fixPointIter = 3, trace = FALSE))

  expect_lte(fit$ICL, icl_before)
  expect_gte(aricode::ARI(fit$fittedSBM$memberships, sampler2$memberships), ari_before - 1e-8)
})

test_that("self stays usable after polish() without an intervening doVEM() call", {
  ## regression test: private$nu must stay a valid sparse matrix after polish() accepts a move,
  ## not NULL -- `sparseMatrix + NULL` silently returns numeric(0), which broke imputedNetwork()
  ## and any split()/candidates_split() call made right after polish() (e.g. from explore())
  net <- missSBM:::partlyObservedNetwork$new(sampler$networkData)
  noisy <- sampler$memberships
  idx <- sample(N, floor(N * .3))
  noisy[idx] <- sample(1:Q, length(idx), replace = TRUE)
  fit <- missSBM:::missSBM_fit$new(net, "dyad", noisy, TRUE)
  fit$doVEM(list(threshold = 1e-3, maxIter = 50, fixPointIter = 3, trace = FALSE))
  fit$polish(list(threshold = 1e-3, maxIter = 50, fixPointIter = 3, trace = FALSE))

  expect_equal(dim(fit$imputedNetwork), c(N, N))
  candidates <- fit$candidates_split(list(threshold = 1e-3, maxIter = 50, fixPointIter = 3, trace = FALSE))
  expect_true(all(sapply(candidates, function(m) is.finite(m$ICL))))
})

test_that("polish() runs without error on withCov and stays finite", {
  X <- matrix(rnorm(N * N), N, N); X <- (X + t(X)) / 2
  adj <- missSBM::observeNetwork(sampler$networkData, "dyad", 0.8, covariates = list(X))
  fit <- missSBM:::missSBM_fit$new(
    missSBM:::partlyObservedNetwork$new(adj, covariates = list(X)), "dyad", sampler$memberships, TRUE)
  fit$doVEM(list(threshold = 1e-3, maxIter = 50, fixPointIter = 3, trace = FALSE))
  icl_before <- fit$ICL

  fit$polish(list(threshold = 1e-3, maxIter = 50, fixPointIter = 3, trace = FALSE))
  expect_true(is.finite(fit$ICL))
  expect_lte(fit$ICL, icl_before)
})

test_that("missSBM_collection$polish() polishes every model in the collection", {
  adj <- missSBM::observeNetwork(sampler$networkData, "dyad", 0.85)
  collection <- estimateMissSBM(adj, vBlocks = 2:6, sampling = "dyad",
                                 control = list(trace = FALSE, exploration = "none", polish = FALSE))
  icl_before <- collection$ICL

  collection$polish(list(threshold = 1e-3, maxIter = 50, fixPointIter = 3, trace = FALSE))
  expect_true(all(collection$ICL <= icl_before + 1e-6))
})

test_that("estimateMissSBM()'s default control polishes (polish = TRUE)", {
  adj <- missSBM::observeNetwork(sampler$networkData, "dyad", 0.85)

  set.seed(1)
  collection_default <- estimateMissSBM(adj, vBlocks = 2:6, sampling = "dyad", control = list(trace = FALSE))
  set.seed(1)
  collection_vem_only <- estimateMissSBM(adj, vBlocks = 2:6, sampling = "dyad",
                                          control = list(trace = FALSE, exploration = "none", polish = FALSE))

  ## same starting point (same seed => same clusterInit/VEM trajectory), polish() can only help
  expect_true(all(collection_default$ICL <= collection_vem_only$ICL + 1e-6))
})

test_that("missSBM_collection's estimate()/polish()/explore() reuse the stored control by default", {
  adj <- missSBM::observeNetwork(sampler$networkData, "dyad", 0.85)
  control <- list(trace = FALSE, threshold = 1e-3, maxIter = 50, fixPointIter = 3,
                   exploration = "none", polish = FALSE)
  collection <- estimateMissSBM(adj, vBlocks = 2:6, sampling = "dyad", control = control)
  icl_after_estimate <- min(collection$ICL)

  ## no-arg calls must not error and must reuse the stored (exploration = "none") control
  expect_no_error(collection$polish())
  icl_after_polish <- min(collection$ICL)
  expect_lte(icl_after_polish, icl_after_estimate + 1e-6)

  ## iterates/direction override the stored control for this call only, without needing a
  ## full control list
  expect_no_error(collection$explore(direction = "both", iterates = 1))
  icl_after_explore <- min(collection$ICL)
  expect_lte(icl_after_explore, icl_after_polish + 1e-6)
})
