#' @importFrom stats binomial glm.fit residuals
init_clustering <- function(adjacencyMatrix, nBlocks, covarArray = NULL, clusterInit = "spectral") {

  N <- nrow(adjacencyMatrix)

  if (nBlocks > 1) {
    if (!is.null(covarArray)) {
      y <- as.vector(adjacencyMatrix)
      X <- cbind(1, apply(covarArray, 3, as.vector))
      adjacencyMatrix <- matrix(NA, N, N)
      NAs <- is.na(y)
      adjacencyMatrix[!NAs] <- logistic(residuals(glm.fit(X[!NAs, ], y[!NAs], family = binomial())))
    }

    if (is.character(clusterInit)) {
      clusterInit <-
        switch(clusterInit,
               "spectral" = init_spectral(    adjacencyMatrix, nBlocks),
               "kmeans"   = init_kmeans(      adjacencyMatrix, nBlocks),
                            init_hierarchical(adjacencyMatrix, nBlocks)
        )
    } else if (is.numeric(clusterInit) | is.factor(clusterInit)) {
      clusterInit <- as.integer(clusterInit)
    } else {
      stop("unknown type for initial clustering")
    }
  } else {
    clusterInit <- rep(1L, N)
  }
  clusterInit
}

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
      cl <- as.integer(kmeans(U, K, iter.max = 100, nstart = 100)$cl)

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
  nBlocks <- length(unique(clustering))
  nNodes  <- length(clustering)
  Z <- matrix(0,nNodes, nBlocks)
  Z[cbind(seq.int(nNodes), clustering)] <- 1
  Z
}

