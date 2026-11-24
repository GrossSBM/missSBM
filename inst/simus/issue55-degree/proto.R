## Issue #55: prototypes of three estimators for the degree sampling (undirected networks),
## written on top of the internals of missSBM (C++ E/M steps of the SBM, degreeSampling_fit),
## without changing the package.
##
## - stage 1: the SBM part on the observed dyads (VEM of a MAR design), what missSBM does now;
## - stage 2: (nu, a, b, zeta) given the SBM part;
## - stage 3: the coordinate ascent of the bound J_zeta in all the variables (the thesis):
##   tau and pi are also updated on the network completed by nu.
##
## A = stages 1 and 2; B = stage 3 from the usual initialization; C = A, then stage 3.
##
## Every stage monitors J_zeta, the bound of the Formularius (eq-bound-degree): factorized
## q(Z) q(Y^m) and the Jaakkola-Jordan bound on the sampling. (The package currently uses, for
## the sampling term, log logistic(a + b E[D_i]), a plug-in which is not a lower bound.)

E_step <- missSBM:::E_step_sparse_bernoulli_nocovariate
M_step <- missSBM:::M_step_sparse_bernoulli_nocovariate
vLL    <- missSBM:::vLL_complete_sparse_bernoulli_nocovariate
h      <- missSBM:::h
cb     <- missSBM:::check_boundaries
xlogx  <- missSBM:::xlogx
softmax_rows <- function(log_tau) {
  log_tau <- log_tau - apply(log_tau, 1, max)
  tau <- exp(log_tau); tau / rowSums(tau)
}

## the state of a fit: Z (tau), theta, alpha, the sampling object and nu
new_state <- function(pon, cl) {
  sbm <- missSBM:::SimpleSBM_fit_noCov$new(pon, cl) # pi from the observed dyads
  Z   <- sbm$probMemberships
  samp <- missSBM:::degreeSampling_fit$new(pon, Z, sbm$connectParam$mean)
  list(pon = pon, Y = pon$networkData, R = pon$samplingMatrix, S = pon$samplingMatrixBar,
       Z = Z, theta = sbm$connectParam$mean, alpha = sbm$blockProp,
       samp = samp, nu = NULL, sbm = sbm)
}

## the imputations of the SBM at the missing dyads: logistic(sum_ql tau_iq tau_jl logit pi_ql)
PI_at_missing <- function(st) {
  missSBM:::.bilinear_at_pattern(st$Z, missSBM:::.logit(st$theta), st$S, missSBM:::.logistic)
}

## the sampling part of J_zeta (undirected): -sum_{N^m} (a + b D_i) + sum_i [log s(zeta_i) +
## (a + b D_i - zeta_i)/2 + h(zeta_i)(a^2 + 2ab D_i + b^2 D2_i - zeta_i^2)], D_i = E[D_i],
## D2_i = E[D_i^2], computed from nu
J_sampling <- function(st) {
  p <- environment(st$samp$update_parameters)$private
  nu <- st$nu
  D  <- p$D_obs + Matrix::rowSums(nu) + Matrix::colSums(nu)
  V  <- nu; V@x <- nu@x * (1 - nu@x)
  D2 <- Matrix::rowSums(V) + Matrix::colSums(V) + D^2
  a <- p$psi[1]; b <- p$psi[2]; z <- p$ksi
  - sum((a + b * D)[!p$N_obs]) +
    sum(log(missSBM:::.logistic(z)) + (a + b * D - z) / 2 + h(z) * (a^2 + 2 * a * b * D + b^2 * D2 - z^2))
}
entropy_nu  <- function(nu) - sum(xlogx(nu@x) + xlogx(1 - nu@x))
entropy_tau <- function(Z) - sum(xlogx(Z))

## J_zeta: E[log p(Y^o, Y^m | Z)] + E[log p(Z)] + sampling part + H(tau) + H(nu)
J_zeta <- function(st) {
  vLL(st$Y + st$nu, st$R + st$S, st$Z, st$theta, st$alpha) + J_sampling(st) +
    entropy_tau(st$Z) + entropy_nu(st$nu)
}
## the bound of stage 1 (observed dyads, no sampling)
J_obs <- function(st) vLL(st$Y, st$R, st$Z, st$theta, st$alpha) + entropy_tau(st$Z)

## generic loop: one_step(st) returns the updated state; stops on the relative change of the
## objective; records it, and counts its decreases
iterate <- function(st, one_step, objective, tol = 1e-7, maxit = 500) {
  obj <- objective(st)
  for (k in seq_len(maxit)) {
    st <- one_step(st)
    obj <- c(obj, objective(st))
    if (abs(diff(tail(obj, 2))) < tol * abs(tail(obj, 1))) break
  }
  st$trace <- obj
  st$decreases <- sum(diff(obj) < -1e-8 * abs(obj[-1]))
  st
}

## stage 1: VEM of the SBM on the observed dyads (tau fixed point, then pi and alpha)
stage1 <- function(st, fixPointIter = 3, ...) {
  iterate(st, function(st) {
    for (k in seq_len(fixPointIter))
      st$Z <- softmax_rows(E_step(st$Y, st$R, st$Z, st$theta, st$alpha, rescale = FALSE))
    m <- M_step(st$Y, st$R, st$Z, TRUE)
    st$theta <- cb(m$theta$mean); st$alpha <- cb(as.numeric(m$pi))
    st
  }, J_obs, ...)
}

## one update of the sampling part: nu given the SBM, then zeta, (a, b)
update_sampling <- function(st) {
  st$nu <- st$samp$update_imputation(PI_at_missing(st))
  st$samp$update_parameters(st$nu)
  st
}

## stage 2: (nu, a, b, zeta) given the SBM part
stage2 <- function(st, ...) {
  if (is.null(st$nu)) st <- update_sampling(st)
  iterate(st, update_sampling, J_zeta, ...)
}

## stage 3: coordinate ascent of J_zeta: nu, (zeta, a, b), tau on the completed network, then pi
## and alpha on the completed network
stage3 <- function(st, fixPointIter = 3, ...) {
  if (is.null(st$nu)) st <- update_sampling(st)
  Rc <- st$R + st$S
  iterate(st, function(st) {
    for (k in seq_len(fixPointIter)) {
      st <- update_sampling(st)
      st$Z <- softmax_rows(E_step(st$Y + st$nu, Rc, st$Z, st$theta, st$alpha, rescale = FALSE))
    }
    m <- M_step(st$Y + st$nu, Rc, st$Z, TRUE)
    st$theta <- cb(m$theta$mean); st$alpha <- cb(as.numeric(m$pi))
    st$samp$update_parameters(st$nu)
    st
  }, J_zeta, ...)
}

## the ICL of degree: -2 (J_zeta - H(tau)) + penalty (as in the package since #54)
ICL_degree <- function(st) {
  Q <- ncol(st$Z); n <- nrow(st$Z)
  pen <- Q * (Q + 1) / 2 * log(sum(st$R)) + (Q - 1) * log(n) + 2 * log(n)
  -2 * (J_zeta(st) - entropy_tau(st$Z)) + pen
}

## the three estimators, from the same initial clustering
fit_A <- function(pon, cl) stage2(stage1(new_state(pon, cl)))
fit_B <- function(pon, cl) stage3(new_state(pon, cl))
fit_C <- function(stA) {stA$samp <- stA$samp$clone(); stage3(stA)} # R6: copy the sampling object
