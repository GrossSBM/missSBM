context("test the Newton-Raphson M-step for the covariate Bernoulli model (M_step_sparse_bernoulli_covariates)")

## The M-step objective for the covariate Bernoulli model is a weighted logistic regression,
## globally concave in (Gamma, beta) jointly: Newton-Raphson should converge in a handful of
## iterations to the unique global maximum, for both the directed and undirected (symmetric
## parameterization of Gamma) cases, and remain robust when the curvature nearly vanishes
## (e.g. once block memberships are nearly hard, as commonly happens late in a VEM run).

simulate_case <- function(seed, N, Q, K, sym, obs_rate = 1) {
  set.seed(seed)
  Z <- matrix(runif(N * Q), N, Q); Z <- Z / rowSums(Z)
  X <- array(rnorm(N * N * K), dim = c(N, N, K))
  Yfull <- matrix(rbinom(N * N, 1, 0.3), N, N); diag(Yfull) <- 0
  Rfull <- matrix(rbinom(N * N, 1, obs_rate), N, N); diag(Rfull) <- 0
  if (sym) {
    Yfull[lower.tri(Yfull)] <- t(Yfull)[lower.tri(Yfull)]
    Rfull[lower.tri(Rfull)] <- t(Rfull)[lower.tri(Rfull)]
  }
  Y <- as(Yfull, "dgCMatrix")
  R <- as(Rfull, "dgCMatrix")
  if (sym) { Y[lower.tri(Y)] <- 0; R[lower.tri(R)] <- 0 }
  list(Y = Y, R = R, X = X, Z = Z, Q = Q, K = K)
}

test_that("Newton solver converges quickly and matches an independent BFGS optimum", {
  for (sym in c(FALSE, TRUE)) {
    case <- simulate_case(1, N = 25, Q = 3, K = 2, sym = sym, obs_rate = 0.7)
    init <- list(Gamma = matrix(0, case$Q, case$Q), beta = rep(0, case$K))

    res <- missSBM:::M_step_sparse_bernoulli_covariates(
      init, case$Y, case$R, case$X, case$Z, sym, 50L, 1e-12
    )
    expect_equal(res$status, 1L)
    expect_lt(res$iterations, 10L)

    Gamma_hat <- log(res$theta$mean / (1 - res$theta$mean))
    if (sym) expect_equal(Gamma_hat, t(Gamma_hat))

    ## independent reference: maximize the same objective with BFGS
    Robs <- as.matrix(case$R) != 0
    objective <- function(Gamma, beta) {
      loglik <- 0
      idx <- which(Robs, arr.ind = TRUE)
      for (r in seq_len(nrow(idx))) {
        i <- idx[r, 1]; j <- idx[r, 2]
        mu  <- sum(beta * case$X[i, j, ])
        eta <- Gamma + mu
        w   <- outer(case$Z[i, ], case$Z[j, ])
        loglik <- loglik + sum(w * (case$Y[i, j] * eta - log1p(exp(eta))))
      }
      loglik
    }
    if (sym) {
      pairs <- which(row(matrix(0, case$Q, case$Q)) <= col(matrix(0, case$Q, case$Q)))
      rc <- arrayInd(pairs, c(case$Q, case$Q))
      pack <- function(par) {
        Gm <- matrix(0, case$Q, case$Q)
        for (p in seq_along(pairs)) { Gm[rc[p,1],rc[p,2]] <- par[p]; Gm[rc[p,2],rc[p,1]] <- par[p] }
        list(Gamma = Gm, beta = par[(length(pairs)+1):length(par)])
      }
      npar <- length(pairs) + case$K
    } else {
      pack <- function(par) list(Gamma = matrix(par[1:(case$Q^2)], case$Q, case$Q), beta = par[(case$Q^2+1):length(par)])
      npar <- case$Q^2 + case$K
    }
    negobj <- function(par) { pp <- pack(par); -objective(pp$Gamma, pp$beta) }
    opt <- optim(rep(0, npar), negobj, method = "BFGS", control = list(reltol = 1e-14, maxit = 2000))
    ref <- pack(opt$par)

    expect_equal(Gamma_hat, ref$Gamma, tolerance = 1e-5)
    expect_equal(as.numeric(res$beta), ref$beta, tolerance = 1e-5)
  }
})

test_that("Newton solver is robust when block memberships are nearly hard (vanishing curvature)", {
  set.seed(99)
  N <- 30; Q <- 3; K <- 2
  Z <- matrix(1e-6, N, Q)
  cl <- sample(1:Q, N, replace = TRUE)
  for (i in 1:N) Z[i, cl[i]] <- 1 - (Q - 1) * 1e-6
  case <- simulate_case(99, N = N, Q = Q, K = K, sym = TRUE, obs_rate = 1)
  case$Z <- Z

  init <- list(Gamma = matrix(0, Q, Q), beta = rep(0, K))
  expect_no_error(
    res <- missSBM:::M_step_sparse_bernoulli_covariates(init, case$Y, case$R, case$X, case$Z, TRUE, 50L, 1e-10)
  )
  expect_equal(res$status, 1L)
})
