#' @importFrom stats sd
init_spectral <- function(X, K) {

  ## basic handling of missing values
  ## handling lonely souls
  n <- ncol(X)
  cl0 <- rep(NA_integer_, n)
  unconnected <- which(rowSums(X, na.rm = TRUE) == 0)
  connected <- setdiff(1:n, unconnected)

  if (K > 1) {
    if (length(connected) == 0) {
      cl0 <- as.integer(base::sample(1:K, n, replace = TRUE))
    } else {
      X <- X[connected,connected]

      D <- rowSums(X, na.rm = TRUE)
      A <- X; A[is.na(A)] <- mean(A, na.rm = TRUE)
      A <- 1/(1 + exp(-A/sd(A)))

      L <- diag(D^(-1/2)) %*% A %*% diag(D^(-1/2))
      U <- eigen(L, symmetric = TRUE)$vectors[,1:K]

      ## Applying the K-means in the eigenspace
      cl <- as.integer(kmeans(U, K, iter.max = 100, nstart = 10)$cl)

      ## handing lonely souls
      cl0[connected] <- cl
      cl0[unconnected] <- which.min(rowsum(D,cl))
    }
  } else {
    cl0 <- rep(1L,n)
  }
  cl0
}

#' @importFrom ape additive
#' @importFrom stats as.dist cutree dist hclust
init_hierarchical <- function(X, K) {
  if (K > 1) {
    D <- as.matrix(dist(X, method = "manhattan"))
    D <- as.dist(ape::additive(D))
    cl0 <- cutree(hclust(D, method = "ward.D2"), K)
  } else {
    cl0 <- rep(1L,nrow(X))
  }
  cl0
}

#' @importFrom stats kmeans
init_kmeans <- function(X, K) {
  if (K > 1) {
    D  <- as.matrix(dist(X, method = "euclidean"))
    # D[which(X == 1)] <- D[which(X == 1)] - 2
    cl0 <- as.integer(kmeans(ape::additive(D), K, nstart = 50, iter.max = 100)$cl)
  } else {
    cl0 <- rep(1L, nrow(X))
  }
  cl0
}

clustering_indicator <- function(clustering) {
  nbBlocks <- length(unique(clustering))
  nbNodes  <- length(clustering)
  if (nbBlocks > 1) {
      Z <- matrix(0,nbNodes, nbBlocks)
      Z[cbind(seq.int(nbNodes), clustering)] <- 1
  } else {
    Z <- matrix(1,nbNodes, nbBlocks)
  }
  Z
}

