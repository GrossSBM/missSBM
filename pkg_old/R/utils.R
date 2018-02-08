zero <- .Machine$double.eps

available_samplings <- c("MAREdge", "MARNode", "snowball", "starDegree", "class", "doubleStandard")

bar <- function(X) {
  X.bar <- 1 - X ; diag(X.bar) <- 0
  X.bar
}

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
