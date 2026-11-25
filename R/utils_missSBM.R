available_samplings <- c("dyad", "covar-dyad", "node", "covar-node", "block-node", "block-dyad", "double-standard", "degree","snowball")

available_samplings_covariates <- c("dyad", "covar-dyad", "node", "covar-node")

## checks of the arguments of the user-facing functions, with informative messages
.check_adjacency <- function(adjacencyMatrix) {
  if (!(is.matrix(adjacencyMatrix) || inherits(adjacencyMatrix, "dgCMatrix")))
    stop("adjacencyMatrix must be a matrix, possibly with NA (or a dgCMatrix).", call. = FALSE)
  if (nrow(adjacencyMatrix) != ncol(adjacencyMatrix))
    stop("adjacencyMatrix must be a square matrix.", call. = FALSE)
}

.check_sampling <- function(sampling, covariates) {
  if (!is.character(sampling) || length(sampling) != 1 || !(sampling %in% available_samplings))
    stop("sampling must be one of ", paste0("'", available_samplings, "'", collapse = ", "), ".", call. = FALSE)
  if (sampling %in% c("covar-dyad", "covar-node") && length(covariates) == 0)
    stop("sampling '", sampling, "' needs covariates.", call. = FALSE)
}

.check_covariates <- function(covariates, nbNodes) {
  if (!is.list(covariates))
    stop("covariates must be a list of size-N vectors (node-centered) or of N x N matrices (dyad-centered).", call. = FALSE)
  ok <- vapply(covariates, function(x)
    if (is.null(dim(x))) length(x) == nbNodes else length(dim(x)) == 2 && all(dim(x) == nbNodes), logical(1))
  if (!all(ok))
    stop("each covariate must be a vector of size N = ", nbNodes, " or an N x N matrix.", call. = FALSE)
}

.check_vBlocks <- function(vBlocks, nbNodes) {
  if (!is.numeric(vBlocks) || length(vBlocks) == 0 || anyNA(vBlocks) ||
      any(vBlocks != round(vBlocks)) || any(vBlocks < 1) || any(vBlocks >= nbNodes))
    stop("vBlocks must be a vector of numbers of blocks: integers between 1 and N - 1 = ", nbNodes - 1, ".", call. = FALSE)
}

.check_clusterInit <- function(clusterInit, nbBlocks, nbNodes) {
  if (!is.list(clusterInit) || length(clusterInit) != nbBlocks)
    stop("clusterInit must be a list with one clustering per number of blocks in vBlocks.", call. = FALSE)
  if (!all(vapply(clusterInit, function(cl) is.numeric(cl) && length(cl) == nbNodes, logical(1))))
    stop("each clustering of clusterInit must be a vector of N = ", nbNodes, " block labels.", call. = FALSE)
}

#' L1-similarity
#'
#' Compute l1-similarity between two vectors
#' @param x a vector
#' @param y a vector
#'
#' @return a vector equal to -abs(x-y)
#' @examples
#' l1_similarity(1:5, 5:1)
#' @export
l1_similarity <- function(x, y) {-abs(x - y)}

clustering_indicator <- function(clustering) {
  K <- max(unique(clustering))
  N  <- length(clustering)
  Z <- matrix(0, N, K)
  Z[cbind(seq.int(N), clustering)] <- 1
  Z
}

## fills empty classes in `labels` (values in 1:K) by moving a node from the currently-largest
## class into each one, so no already-valid class is emptied in the process
repair_empty_classes <- function(labels, K) {
  stopifnot(length(labels) >= K)
  counts <- tabulate(labels, nbins = K)
  for (lab in which(counts == 0)) {
    donor  <- which.max(counts)
    victim <- sample(which(labels == donor), 1)
    labels[victim] <- lab
    counts[donor]  <- counts[donor] - 1
    counts[lab]    <- 1
  }
  labels
}

## TRUE if a missSBM_fit's clustering has fewer occupied classes than its structural nbBlocks
is_degenerate <- function(fit) {
  fit$occupiedBlocks < fit$fittedSBM$nbBlocks
}

## the missSBM_fit with the smallest ICL among a list of candidates (the criterion used
## throughout to rank models)
best_by_icl <- function(candidates) {
  candidates[[which.min(sapply(candidates, function(m) m$ICL))]]
}

## future_lapply() with this package's policy for parallel candidate generation: draws are
## reproducible (future.seed) but workers are assigned candidates in random order, so that a
## worker that draws a slow candidate doesn't systematically hold up the same downstream slot
future_lapply_shuffled <- function(X, FUN, ...) {
  future_lapply(X, FUN, ..., future.seed = TRUE, future.scheduling = structure(TRUE, ordering = "random"))
}

## logical N x N mask selecting each dyad exactly once: all off-diagonal pairs if directed, only
## the upper triangle (i < j) if undirected
## the patterns of the observed dyads (R), of the missing ones (S) and of the edges (Y) of a
## network given as a dense matrix with NA or as a dgCMatrix, as sparse matrices on the valid
## dyads (upper triangle if undirected), and whether it is directed, without any n x n temporary
## (no mask, copy or transpose of the input). A dense matrix is read in place in C++
## (network_patterns_dense()); a sparse one (possibly with NA), column by column. Also checks that the
## observed values (diagonal included) are 0 or 1, both present.
.network_patterns <- function(A) {
  n <- ncol(A)
  if (inherits(A, "sparseMatrix")) {
    directed <- !Matrix::isSymmetric(A)
    obs <- miss <- edges <- vector("list", n)
    ## the diagonal is not a dyad, but its values are checked too (as by the former check)
    d <- Matrix::diag(A); d <- d[!is.na(d)]
    if (any(d != 0 & d != 1)) stop("the network must be binary (0, 1, or NA for missing)")
    values <- c(zero = any(d == 0), one = any(d == 1))
    for (j in seq_len(n)) {
      i <- if (directed) seq_len(n)[-j] else seq_len(j - 1)
      v <- as.numeric(A[i, j])
      na <- is.na(v); v <- v[!na]
      if (any(v != 0 & v != 1)) stop("the network must be binary (0, 1, or NA for missing)")
      values <- values | c(any(v == 0), any(v == 1))
      obs[[j]]   <- i[!na]
      miss[[j]]  <- i[na]
      edges[[j]] <- i[!na][v != 0]
    }
    pointers <- function(rows) c(0L, cumsum(lengths(rows)))
    p <- list(R_i = unlist(obs) - 1L, R_p = pointers(obs), S_i = unlist(miss) - 1L, S_p = pointers(miss),
              Y_i = unlist(edges) - 1L, Y_p = pointers(edges),
              directed = directed, binary = TRUE, has_zero = values[["zero"]], has_one = values[["one"]])
  } else {
    p <- network_patterns_dense(A)
  }
  if (!p$binary) stop("the network must be binary (0, 1, or NA for missing)")
  if (!(p$has_zero && p$has_one)) stop("the network must have both 0 and 1 among its observed dyads")
  ## the indices come column by column, rows increasing: the compressed format directly, without
  ## the sort of the triplets of sparseMatrix() (about 3 times the result)
  pattern <- function(i, p) new("dgCMatrix", i = i, p = as.integer(p), x = rep(1, length(i)), Dim = c(n, n))
  list(R = pattern(p$R_i, p$R_p), S = pattern(p$S_i, p$S_p), Y = pattern(p$Y_i, p$Y_p), directed = p$directed)
}

valid_dyads <- function(dim, directed) {
  if (directed) .row(dim) != .col(dim) else .row(dim) < .col(dim)
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

.logistic <- function(x) {1/(1 + exp(-x))}
.logit    <- function(x) {log(x/(1 - x))}
xlogx <- function(x) {
  out <- x * log(x)
  out[x < .Machine$double.eps] <- 0
  out
}

## values of a dense matrix at the structural nonzeros of a dgCMatrix pattern, packaged as a
## dgCMatrix with that same pattern. Avoids `as(dense, "dgCMatrix") * pattern`: that Ops-based
## route goes through Matrix's generic sparse-arithmetic dispatch (intersection of nonzero
## indices, class promotion) even though the dense operand has no actual sparsity to exploit.
## `fun` is applied *after* extraction, not to the whole dense matrix: when pattern is much
## sparser than dense (e.g. a minority of dyads missing), this keeps the elementwise transform
## (typically the costly part, e.g. exp()-based) down to O(nnz(pattern)) instead of O(N^2).
.mask_dense_at_pattern <- function(dense, pattern, fun = identity) {
  cidx <- rep.int(seq_len(ncol(pattern)), diff(pattern@p))
  pattern@x <- fun(dense[cbind(pattern@i + 1L, cidx)])
  pattern
}

## values of Z %*% M %*% t(Z) at the structural nonzeros of a dgCMatrix pattern, packaged as a
## dgCMatrix with that same pattern, without forming the N x N matrix nor any (nnz x Q)
## temporary: O(nnz(pattern) * Q) time, the size of the result in memory (C++).
.bilinear_at_pattern <- function(Z, M, pattern, fun = identity) {
  Z <- as.matrix(Z); storage.mode(Z) <- "double"
  pattern@x <- fun(bilinear_at_pattern_cpp(t(Z %*% M), t(Z), pattern@i, pattern@p))
  pattern
}

.softmax <- function(x) {
  b <- max(x)
  exp(x - b) / sum(exp(x - b))
}

## the softmax of each row of a matrix of log-probabilities (up to a constant per row)
.softmax_rows <- function(log_tau) {t(apply(log_tau, 1, .softmax))}

## the MAP memberships of the probabilities tau, as an indicator matrix (empty classes kept)
.map_indicator <- function(tau) {
  Z <- matrix(0, nrow(tau), ncol(tau))
  Z[cbind(seq_len(nrow(tau)), max.col(tau, ties.method = "first"))] <- 1
  Z
}

## the counts of a (sparse) n x n matrix X by pair of blocks, weighted by the memberships Z:
## (Z' X Z)_ql = sum_ij Z_iq X_ij Z_jl, a Q x Q matrix
.block_counts <- function(Z, X) {as.matrix(t(Z) %*% X %*% Z)}

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
