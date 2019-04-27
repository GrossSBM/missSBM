zero <- .Machine$double.eps

available_samplings <- c("dyad", "covar-dyad", "node", "covar-node", "block-node", "block-dyad", "double-standard", "degree")

available_samplings_covariates <- c("dyad", "covar-dyad", "node", "covar-node")

l1_similarity <- function(x, y) {-abs(x - y)}

bar <- function(X) {
  X.bar <- 1 - X ; diag(X.bar) <- 0
  X.bar
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

logistic <- function(x) {1/(1 + exp(-x))}
logit    <- function(x) {log(x/(1 - x))}

.softmax <- function(x) {
  b <- max(x)
  exp(x - b) / sum(exp(x - b))
}

h <- function(x) {-.5 * (logistic(x) - 0.5) / x}

xlogx <- function(x) ifelse(x < .Machine$double.eps, 0, x*log(x))

check_boundaries <- function(x, zero = .Machine$double.eps) {
  x[is.nan(x)] <- zero
  x[x > 1 - zero] <- 1 - zero
  x[x <     zero] <-     zero
  x
}
#'
#' #' @importFrom graphics box image par
#' image_NA <- function(z,  zlim = c(0,1), col = c("white", "midnightblue"), na.color = 'gray', outside.below.color = 'black', outside.above.color = 'white', ...)
#' {
#'   zstep <- (zlim[2] - zlim[1]) / length(col); # step in the color palette
#'   newz.below.outside <- zlim[1] - 2 * zstep # new z for values below zlim
#'   newz.above.outside <- zlim[2] + zstep # new z for values above zlim
#'   newz.na <- zlim[2] + 2 * zstep # new z for NA
#'
#'   z[which(z < zlim[1])] <- newz.below.outside # we affect newz.below.outside
#'   z[which(z > zlim[2])] <- newz.above.outside # we affect newz.above.outside
#'   z[which(is.na(z > zlim[2]))] <- newz.na # same for newz.na
#'
#'   zlim[1] <- zlim[1] - 2 * zstep # extend lower limit to include below value
#'   zlim[2] <- zlim[2] + 2 * zstep # extend top limit to include the two new values above and na
#'
#'   col <- c(outside.below.color, col[1], col, outside.above.color, na.color) #correct by including col[1] at bottom of range
#'
#'   par(mar = c(2.1,8.1,3.1,3.1))
#'   image(z[nrow(z):1,],  zlim = zlim, col = col, xaxt = "n", yaxt = "n", ...) # we finally call image(...)
#'   box()
#'   # par(mar=c(5.1,4.1,4.1,2.1))
#' }
#'
