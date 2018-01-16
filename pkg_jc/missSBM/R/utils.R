zero <- .Machine$double.eps

available_samplings <- c("dyad", "node", "snowball", "degree", "block", "double_standard")

bar <- function(X) {
  X.bar <- 1 - X ; diag(X.bar) <- 0
  X.bar
}

logistic <- function(x) 1/(1 + exp(-x))

logx <- function(x) ifelse(x < .Machine$double.eps, 0, log(x))

log1mx <- function(x) ifelse(x > 1 - .Machine$double.eps, 0, log(1 - x))

init_spectral <- function(X, K) {

  ## basic handling of missing values
  if (anyNA(X)) X[is.na(X)] <- 0

  ## handling lonely souls
  cl.final <- rep(NA, ncol(X))
  unconnected <- which(rowSums(X) == 0)
  connected <- setdiff(1:ncol(X), unconnected)

  X <- X[connected,connected]
  if (K > 1) {

    ## Normalized Laplacian
    D <- colSums(X)
    L <- diag(rep(1,ncol(X))) -
      diag(D^(-1/2)) %*% X %*% diag(D^(-1/2))

    ## Absolute eigenvalue in order
    E <- order(-abs(eigen(L)$values))

    ## Go into eigenspace
    U <- eigen(L)$vectors[,E]
    U <- U[,c((ncol(U) - K + 1):ncol(U))]
    U <- U / rowSums(U^2)^(1/2)

    ## Applying the K-means in the eigenspace
    cl <- kmeans(U, K, nstart = 10, iter.max = 30)$cluster
  } else {
    cl <- as.factor(rep(1,ncol(X)))
  }

  ## handing lonely souls
  cl.final[connected] <- cl
  cl.final[unconnected] <- which.min(rowsum(D,cl))

  as.factor(cl.final)
}

init_hierarchical <- function(X, K) {

  ## basic handling of missing values
  if (anyNA(X)) X[is.na(X)] <- 0

  D  <- as.matrix(dist(X, method = "manhattan"))
  D[X == 1] <- D[X == 1] - 2 ## does not make sense for Poisson models ..
  cl0 <- cutree(hclust(as.dist(D), method = "ward.D"), K)
  as.factor(cl0)
}

init_kmeans <- function(X, K) {
  cl0 <- kmeans(X, K)$cl
  as.factor(cl0)
}

#' @export
image_NA <- function(z,  zlim = c(0,1), col = c("white", "midnightblue"), na.color = 'gray', outside.below.color = 'black', outside.above.color = 'white', ...)
{
  zstep <- (zlim[2] - zlim[1]) / length(col); # step in the color palette
  newz.below.outside <- zlim[1] - 2 * zstep # new z for values below zlim
  newz.above.outside <- zlim[2] + zstep # new z for values above zlim
  newz.na <- zlim[2] + 2 * zstep # new z for NA

  z[which(z < zlim[1])] <- newz.below.outside # we affect newz.below.outside
  z[which(z > zlim[2])] <- newz.above.outside # we affect newz.above.outside
  z[which(is.na(z > zlim[2]))] <- newz.na # same for newz.na

  zlim[1] <- zlim[1] - 2 * zstep # extend lower limit to include below value
  zlim[2] <- zlim[2] + 2 * zstep # extend top limit to include the two new values above and na

  col <- c(outside.below.color, col[1], col, outside.above.color, na.color) #correct by including col[1] at bottom of range

  par(mar = c(2.1,8.1,3.1,3.1))
  image(z[nrow(z):1,],  zlim = zlim, col = col, xaxt = "n", yaxt = "n", ...) # we finally call image(...)
  box()
  # par(mar=c(5.1,4.1,4.1,2.1))
}
