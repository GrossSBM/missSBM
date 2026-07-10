zero <- .Machine$double.eps


available_samplings <- c("dyad", "covar-dyad", "node", "covar-node", "block-node", "block-dyad", "double-standard", "degree","snowball")

available_samplings_covariates <- c("dyad", "covar-dyad", "node", "covar-node")

#' L1-similarity
#'
#' Compute l1-similarity between two vectors
#' @param x a vector
#' @param y a vector
#'
#' @return a vector equal to -abs(x-y)
#' @export
l1_similarity <- function(x, y) {-abs(x - y)}

clustering_indicator <- function(clustering) {
  K <- max(unique(clustering))
  N  <- length(clustering)
  Z <- matrix(0, N, K)
  Z[cbind(seq.int(N), clustering)] <- 1
  Z
}

array2list <-function(X) {
  if (is.null(X)) {
    L <- list()
  } else {
    D <- dim(X)[3]
    L <- vector("list", length = D)
    for (i in 1:D) L[[i]] <- X[,,i]
  }
  L
}

format_covariates <- function(covariates, similarity) {
  if (length(covariates) > 0) {
    # Conversion of covariates to an array
    covariates <- simplify2array(covariates)
    # if a list of vector (covariates node-centered), will be a matrix
    # and thus must be node centered
    if (is.matrix(covariates)) {
      covarMatrix <- covariates
      covarArray  <- getCovarArray(covarMatrix, similarity)
    }
    # if a list of matrix (covariates dyad-centered), will be a 3-dimensional array
    # and thus must be dyad centered
    if (length(dim(covariates)) == 3) {
      covarMatrix <- NULL
      covarArray  <- covariates
    }
  } else {
    covarMatrix <- NULL
    covarArray  <- NULL
  }
  res <- list(Matrix = covarMatrix, Array = covarArray)
  res
}

getCovarArray <- function(X, s) {
  if (is.null(X))
    return(NULL)
  N <- nrow(X)
  M <- ncol(X)
  if (identical(s, l1_similarity)) {
    ## fast, fully vectorized path for the default similarity -abs(x - y):
    ## avoids an O(N^2) R-level double loop
    phi <- vapply(seq_len(M), function(m) -abs(outer(X[, m], X[, m], "-")), matrix(0, N, N))
  } else {
    phi <- array(dim = c(N, N, M))
    for (i in 1:N)
      for (j in 1:N)
        phi[i,j,] <- s(X[i, ], X[j, ])
  }
  phi
}

#'
#' @importFrom Matrix drop0
dropNA <- function(x) {
    if(!is(x, "matrix")) stop("x needs to be a matrix!")

    zeros <- which(x==0, arr.ind=TRUE)
    ## keep zeros
    x[is.na(x)] <- 0
    x[zeros] <- NA
    x <- drop0(x)
    x[zeros] <- 0
    x
}

.logistic <- function(x) {1/(1 + exp(-x))}
.logit    <- function(x) {log(x/(1 - x))}
xlogx     <- function(x) ifelse(x < .Machine$double.eps, 0, x*log(x))

.softmax <- function(x) {
  b <- max(x)
  exp(x - b) / sum(exp(x - b))
}

h <- function(x) {-.5 * (.logistic(x) - 0.5) / x}

#' Shared variational-EM driver
#'
#' Runs the E-step/M-step loop common to [`SimpleSBM_fit`] and [`missSBM_fit`], monitoring
#' convergence and stepping back to the previous state if the objective decreases. Behavior
#' (fixed-point E-step, convergence criterion, step-back) is factorized here so both classes
#' stay in sync; what differs between them (how a state is snapshotted/restored, how many models
#' are involved) is passed in as closures.
#'
#' @param control a list with components \code{threshold}, \code{maxIter}, \code{fixPointIter}, \code{trace}
#' @param init_stop logical, should the loop be skipped altogether (e.g. a single block, nothing to estimate)
#' @param e_step a function(fixPointIter) performing the (repeated) variational E-step, called for its side effects
#' @param m_step a function() performing the M-step, called for its side effects
#' @param get_loglik a function() returning the current (variational) log-likelihood
#' @param get_theta a function() returning the current connectivity parameter, used to assess convergence
#' @param snapshot a function() capturing the current mutable state, to be passed to \code{restore()}
#' @param restore a function(state) resetting the mutable state to a previous \code{snapshot()}
#' @param reorder a function() called once convergence is reached, to permute group labels
#' @param get_flat_state optional function() returning a flat, unconstrained-space parameter
#'   vector to accelerate with SQUAREM (see \code{run_VEM_accelerated()}); if \code{NULL}
#'   (together with \code{set_flat_state}), the plain (unaccelerated) driver is used -- the
#'   default, and the only path exercised by callers that don't opt in
#' @param set_flat_state optional function(p) rebuilding the model from a flat vector as
#'   returned by \code{get_flat_state()}
#' @return a data.frame with the convergence monitoring (iteration, delta_parameters, delta_objective, elbo)
#' @noRd
run_VEM <- function(control, init_stop, e_step, m_step, get_loglik, get_theta, snapshot, restore, reorder,
                     get_flat_state = NULL, set_flat_state = NULL) {
  if (is.null(get_flat_state) || is.null(set_flat_state)) {
    return(run_VEM_classic(control, init_stop, e_step, m_step, get_loglik, get_theta, snapshot, restore, reorder))
  }
  run_VEM_accelerated(control, init_stop, e_step, m_step, get_loglik, get_theta, snapshot, restore, reorder,
                       get_flat_state, set_flat_state)
}

run_VEM_classic <- function(control, init_stop, e_step, m_step, get_loglik, get_theta, snapshot, restore, reorder) {

  delta_par <- vector("numeric", control$maxIter); delta_par[1] <- NA
  delta_obj <- vector("numeric", control$maxIter); delta_obj[1] <- NA
  objective <- vector("numeric", control$maxIter); objective[1] <- get_loglik()
  iterate <- 1; stop <- init_stop

  while (!stop) {
    iterate <- iterate + 1
    if (control$trace) cat(" iteration #:", iterate, "\r")

    theta_old <- get_theta()
    state_old <- snapshot()

    ## Variational E-step
    e_step(control$fixPointIter)

    ## M-step
    m_step()

    ## Assess convergence
    objective[iterate] <- get_loglik()
    delta_par[iterate] <- sqrt(sum((get_theta() - theta_old)^2)) / sqrt(sum((theta_old)^2))
    delta_obj[iterate] <- abs(objective[iterate] - objective[iterate-1]) / abs(objective[iterate])
    stop <- (iterate > control$maxIter) | ((delta_par[iterate] < control$threshold) & (delta_obj[iterate] < control$threshold))

    ## Step back if the objective decreases
    if (objective[iterate] < objective[iterate-1] & iterate > 2) {
      restore(state_old)
      iterate <- iterate - 1
      stop <- TRUE
    }
  }
  reorder()
  if (control$trace) cat("\n")

  data.frame(iteration = 1:iterate, delta_parameters = delta_par[1:iterate], delta_objective = delta_obj[1:iterate], elbo = objective[1:iterate])
}

#' SQUAREM-accelerated variant of run_VEM_classic()
#'
#' Every cycle runs up to two plain VEM steps (p0 -> p1 -> p2 in the model's flat, unconstrained
#' parameterization), then attempts a SQUAREM extrapolation (Varadhan & Roland, 2008) beyond p2,
#' stabilized by one more E-step/M-step pass. Step-back is per *cycle*: rolled back only if the
#' whole cycle ends up worse than before it started.
#' @noRd
run_VEM_accelerated <- function(control, init_stop, e_step, m_step, get_loglik, get_theta, snapshot, restore, reorder,
                                 get_flat_state, set_flat_state) {

  delta_par <- vector("numeric", control$maxIter); delta_par[1] <- NA
  delta_obj <- vector("numeric", control$maxIter); delta_obj[1] <- NA
  objective <- vector("numeric", control$maxIter); objective[1] <- get_loglik()
  iterate <- 1; stop <- init_stop

  while (!stop && iterate < control$maxIter) {
    iterate0  <- iterate
    theta_ref <- get_theta()
    state_ref <- snapshot()
    obj_ref   <- objective[iterate]

    p0 <- get_flat_state()
    e_step(control$fixPointIter); m_step()
    iterate <- iterate + 1
    objective[iterate] <- get_loglik()
    if (control$trace) cat(" iteration #:", iterate, "\r")

    if (abs(objective[iterate] - obj_ref) >= control$threshold && iterate < control$maxIter) {
      p1 <- get_flat_state()
      e_step(control$fixPointIter); m_step()
      iterate <- iterate + 1
      objective[iterate] <- get_loglik()
      if (control$trace) cat(" iteration #:", iterate, "\r")
      p2   <- get_flat_state()
      obj2 <- objective[iterate]

      if (try_squarem_step(p0, p1, p2, obj2, get_flat_state, set_flat_state,
                            e_step, m_step, get_loglik, control$fixPointIter)) {
        iterate <- iterate + 1
        objective[iterate] <- get_loglik()
        if (control$trace) cat(" iteration #:", iterate, "\r")
      }
    }

    delta_par[iterate] <- sqrt(sum((get_theta() - theta_ref)^2)) / sqrt(sum((theta_ref)^2))
    delta_obj[iterate] <- abs(objective[iterate] - obj_ref) / abs(objective[iterate])
    stop <- (delta_par[iterate] < control$threshold) & (delta_obj[iterate] < control$threshold)

    ## Step back the whole cycle if it made things worse overall
    if (objective[iterate] < obj_ref) {
      restore(state_ref)
      iterate <- iterate0
      stop <- TRUE
    }
  }
  reorder()
  if (control$trace) cat("\n")

  data.frame(iteration = 1:iterate, delta_parameters = delta_par[1:iterate], delta_objective = delta_obj[1:iterate], elbo = objective[1:iterate])
}

#' One SQUAREM extrapolation attempt from three consecutive plain-VEM flat states p0 -> p1 -> p2
#' (model currently at p2, get_loglik() == obj2). Uses Varadhan & Roland's steplength "S3" with
#' backtracking (halving the distance to alpha = -1, i.e. no extrapolation) until the stabilized
#' candidate is at least as good as obj2. Returns TRUE if the model was left at the
#' accelerated+stabilized state, FALSE if restored to p2.
#' @noRd
try_squarem_step <- function(p0, p1, p2, obj2, get_flat_state, set_flat_state, e_step, m_step, get_loglik, fixPointIter) {
  r <- p1 - p0
  v <- (p2 - p1) - r
  vv <- sum(v^2)
  if (vv < 1e-14) return(FALSE) # already essentially in the linear/converged regime

  alpha <- -sqrt(sum(r^2) / vv) # steplength "S3" of Varadhan & Roland (2008)
  accept_slack <- 1e-8 * (1 + abs(obj2)) # tolerate floating-point-scale noise only

  success <- FALSE
  for (attempt in seq_len(10)) {
    if (alpha >= -1) break
    p_new <- p0 - 2 * alpha * r + alpha^2 * v
    set_flat_state(p_new)
    e_step(fixPointIter)
    m_step()
    if (get_loglik() >= obj2 - accept_slack) { success <- TRUE; break }
    alpha <- 0.5 * (alpha - 1) # halve the distance to alpha = -1 (no extrapolation)
  }
  if (!success) {
    set_flat_state(p2)
    e_step(fixPointIter)
  }
  success
}

check_boundaries <- function(x, zero = .Machine$double.eps) {
  x[is.nan(x)] <- zero
  x[x > 1 - zero] <- 1 - zero
  x[x <     zero] <-     zero
  x
}

#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL
