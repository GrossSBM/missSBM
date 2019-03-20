#' @importFrom stats binomial glm.fit residuals
init_clustering <- function(adjacencyMatrix, nBlocks, covarArray = NULL, clusterInit = "hierarchical") {

  N <- nrow(adjacencyMatrix)

  if (nBlocks > 1) {
    if (!is.null(covarArray)) {
      y <- as.vector(adjacencyMatrix)
      X <- apply(covarArray, 3, as.vector)
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
      ## Normalized Laplacian
      t <- sapply(1:length(connected), function(i) mean(1*(!is.na(X[i,]))))
      D <- colSums(X, na.rm = TRUE)/t
      A <- X; A[is.na(A)] <- mean(A, na.rm = TRUE)
      # A <- A + (1/4)*mean(rowSums(A))/nrow(A)*matrix(1, nrow(A), nrow(A))
      L <- diag(rep(1,length(t))) - diag(D^(-1/2)) %*% A %*% diag(D^(-1/2))
      ## Absolute eigenvalue in order
      E <- order(-abs(eigen(L)$values))

      ## Go into eigenspace
      U <- eigen(L)$vectors[,E]
      U <- U[,c((ncol(U) - K + 1):ncol(U))]
      U <- U / rowSums(U^2)^(1/2)
      U[is.na(U)] <- 0

      ## Applying the K-means in the eigenspace
      cl <- kmeans(U, K, nstart = 10, iter.max = 50)$cluster
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
    D  <- as.matrix(dist(X, method = "manhattan"))
    D[which(X == 1)] <- D[which(X == 1)] - 2
    D <- as.dist(ape::additive(D))
    cl0 <- cutree(hclust(as.dist(D), method = "ward.D"), K)
  } else {
    cl0 <- rep(1L,ncol(X))
  }
  cl0
}

#' @importFrom stats kmeans
init_kmeans <- function(X, K) {
  if (K > 1) {
    D  <- as.matrix(dist(X, method = "euclidean"))
    # D[which(X == 1)] <- D[which(X == 1)] - 2
    cl0 <- kmeans(ape::additive(D), K, nstart = 10, iter.max = 50)$cl
  } else {
    cl0 <- rep(1L, ncol(X))
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

