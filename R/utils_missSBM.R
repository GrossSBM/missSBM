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
  phi <- array(dim = c(N, N, M))
  for (i in 1:N)
    for (j in 1:N)
      phi[i,j,] <- s(X[i, ], X[j, ])
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
#' @return a data.frame with the convergence monitoring (iteration, delta_parameters, delta_objective, elbo)
#' @noRd
run_VEM <- function(control, init_stop, e_step, m_step, get_loglik, get_theta, snapshot, restore, reorder) {

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
