context("test SQUAREM acceleration of the VEM (SimpleSBM_fit_noCov)")

## SimpleSBM_fit_noCov is the only class validated so far to support SQUAREM acceleration
## (run_VEM() dispatches to run_VEM_accelerated() when get_flat_state()/set_flat_state() are
## supplied, see R/utils_missSBM.R and R/R6Class-simpleSBM_fit.R). Correctness (same fixed
## point) is not *guaranteed* for every start -- SQUAREM changes the optimization trajectory,
## and the SBM likelihood is genuinely multimodal for larger Q -- but it should coincide with
## plain VEM in typical, non-adversarial cases while needing markedly fewer iterations.

test_that("supports_acceleration() is opted in only where validated", {
  net <- missSBM:::partlyObservedNetwork$new(matrix(c(NA,1,0,1,NA,1,0,1,NA), 3, 3))
  expect_true(missSBM:::SimpleSBM_fit_noCov$new(net, c(1,1,2))$supports_acceleration())
  expect_false(missSBM:::SimpleSBM_fit_MNAR$new(net, c(1,1,2))$supports_acceleration())
  ## SimpleSBM_fit_withCov does not override supports_acceleration, so it inherits the
  ## FALSE default straight from the (virtual) SimpleSBM_fit base class -- checked on the
  ## generator directly rather than instantiating (which requires actual covariates)
  expect_null(missSBM:::SimpleSBM_fit_withCov$public_methods$supports_acceleration)
})

test_that("get_flat_state()/set_flat_state() round-trip theta/pi exactly", {
  net <- missSBM:::partlyObservedNetwork$new(matrix(c(NA,1,0,1,1,NA,1,0,0,1,NA,1,1,0,1,NA), 4, 4))
  fit <- missSBM:::SimpleSBM_fit_noCov$new(net, c(1,1,2,2))
  p <- fit$get_flat_state()
  theta0 <- fit$connectParam$mean
  pi0    <- fit$blockProp
  fit$set_flat_state(p)
  expect_equal(fit$connectParam$mean, theta0, tolerance = 1e-10)
  expect_equal(fit$blockProp, pi0, tolerance = 1e-10)
})

test_that("accelerated VEM reaches the same fixed point as plain VEM, in markedly fewer iterations", {
  set.seed(1)
  N <- 400; Q <- 20
  sampler <- sbm::SimpleSBM$new("bernoulli", N, FALSE, rep(1/Q, Q), list(mean = diag(.35, Q) + .05))
  sampler$rMemberships(store = TRUE)
  sampler$rNetwork(store = TRUE)
  net <- missSBM:::partlyObservedNetwork$new(sampler$networkData)

  ## perturb the true membership to force a non-trivial number of VEM iterations
  noisy <- sampler$memberships
  idx <- sample(N, floor(N * .3))
  noisy[idx] <- sample(1:Q, length(idx), replace = TRUE)

  fit_accelerated <- missSBM:::SimpleSBM_fit_noCov$new(net, noisy)
  fit_classic      <- missSBM:::SimpleSBM_fit_noCov$new(net, noisy)

  out_accelerated <- fit_accelerated$doVEM(threshold = 1e-6, maxIter = 300, fixPointIter = 3, trace = FALSE)
  out_classic <- missSBM:::run_VEM_classic(
    control    = list(threshold = 1e-6, maxIter = 300, fixPointIter = 3, trace = FALSE),
    init_stop  = fit_classic$nbBlocks <= 1,
    e_step     = function(fixPointIter) for (i in seq.int(fixPointIter)) fit_classic$update_blocks(),
    m_step     = function() fit_classic$update_parameters(),
    get_loglik = function() fit_classic$loglik,
    get_theta  = function() fit_classic$connectParam$mean,
    snapshot   = function() fit_classic$get_state(),
    restore    = function(state) fit_classic$set_state(state),
    reorder    = function() fit_classic$reorder()
  )

  expect_lt(nrow(out_accelerated), nrow(out_classic))
  expect_equal(tail(out_accelerated$elbo, 1), tail(out_classic$elbo, 1), tolerance = 1e-3)
  expect_gt(aricode::ARI(fit_accelerated$memberships, fit_classic$memberships), 0.99)
})
