zero <- .Machine$double.eps


available_samplings <- c("dyad", "covar-dyad", "node", "covar-node", "block-node", "block-dyad", "double-standard", "degree","snowball")

available_samplings_covariates <- c("dyad", "covar-dyad", "node", "covar-node")

l1_similarity <- function(x, y) {-abs(x - y)}

bar <- function(X) {
  X.bar <- 1 - X ; diag(X.bar) <- 0
  X.bar
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
  if (!is.null(covariates)) {
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

quad_form <- function(A,x) {t(x) %*% A %*% x}

.logistic <- function(x) {1/(1 + exp(-x))}
.logit    <- function(x) {log(x/(1 - x))}

.softmax <- function(x) {
  b <- max(x)
  exp(x - b) / sum(exp(x - b))
}

h <- function(x) {-.5 * (.logistic(x) - 0.5) / x}

xlogx <- function(x) ifelse(x < .Machine$double.eps, 0, x*log(x))

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
