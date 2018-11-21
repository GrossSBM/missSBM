graphCAH_old <- function(X, Q) {
  D  <- as.matrix(dist(X, method="manhattan"))
  D[which(X == 1)] <- D[which(X == 1)] - 2
  is (anyNA(D))
  D <- as.dist(additive(D))
  return(cutree(hclust(as.dist(D), method="ward.D"), Q))
}

graphCAH_new <- function(X, Q) {
  D  <- as.matrix(dist(X, method="manhattan"))
  D[which(X == 1)] <- D[which(X == 1)] - 2
  is (anyNA(D))
  D <- as.dist(additive(D))
  return(cutree(hclust(as.dist(D), method="ward.D"), Q))
}

SpectralClustering_NAisZero <- function (A,Q) {
  if(Q > 1){
    ## basic handling of missing values
    if (anyNA(A)) A[is.na(A)] <- 0

    ## handling lonely souls
    cl.final <- rep(NA, ncol(A))
    unconnected <- which(rowSums(A) == 0)
    connected <- setdiff(1:ncol(A), unconnected)

    A <- A[connected,connected]
    if (Q > 1) {

      ## Normalized Laplacian
      D <- colSums(A)
      L <- diag(rep(1,ncol(A))) -
        diag(D^(-1/2)) %*% A %*% diag(D^(-1/2))

      ## Absolute eigenvalue in order
      E <- order(-abs(eigen(L)$values))

      ## Go into eigenspace
      U <- eigen(L)$vectors[,E]
      U <- U[,c((ncol(U)-Q+1):ncol(U))]
      U <- U / rowSums(U^2)^(1/2)

      ## Applying the K-means in the eigenspace
      cl <- kmeans(U, Q, nstart = 10, iter.max = 30)$cluster
    } else {
      cl <- as.factor(rep(1,ncol(A)))
    }

    ## handing lonely souls
    cl.final[connected] <- cl
    cl.final[unconnected] <- which.min(rowsum(D,cl))
    ##   cl.final[unconnected] <- sample(levels(factor(cl)), length(unconnected), rep=TRUE, prob=table(cl)/sum(table(cl)))

    return(as.factor(cl.final))
  } else {
    return(rep(1, ncol(A)))
  }
}

SpectralClustering_NAisMean <- function (A,Q) {
  if(Q > 1){
    ## basic handling of missing values
    if (anyNA(A)) A[is.na(A)] <- mean(A[!is.na(A)])

    ## handling lonely souls
    cl.final <- rep(NA, ncol(A))
    unconnected <- which(rowSums(A) == 0)
    connected <- setdiff(1:ncol(A), unconnected)

    A <- A[connected,connected]

    ## Normalized Laplacian
    D <- colSums(A)
    L <- diag(rep(1,ncol(A))) -
      diag(D^(-1/2)) %*% A %*% diag(D^(-1/2))

    ## Absolute eigenvalue in order
    E <- order(-abs(eigen(L)$values))

    ## Go into eigenspace
    U <- eigen(L)$vectors[,E]
    U <- U[,c((ncol(U)-Q+1):ncol(U))]
    U <- U / rowSums(U^2)^(1/2)

    ## Applying the K-means in the eigenspace
    cl <- kmeans(U, Q, nstart = 10, iter.max = 30)$cluster

    ## handing lonely souls
    cl.final[connected] <- cl
    cl.final[unconnected] <- which.min(rowsum(D,cl))
    ##   cl.final[unconnected] <- sample(levels(factor(cl)), length(unconnected), rep=TRUE, prob=table(cl)/sum(table(cl)))

    return(as.factor(cl.final))
  } else {
    return(rep(1, ncol(A)))
  }
}

SparseSpectralClustering_NAisZero <- function (A,Q) {

  ## basic handling of missing values
  if (anyNA(A)) A[is.na(A)] <- 0

  ## handling lonely souls
  cl.final <- rep(NA, ncol(A))
  unconnected <- which(rowSums(A) == 0)
  connected <- setdiff(1:ncol(A), unconnected)

  A <- A[connected,connected]
  if (Q > 1) {

    ## Normalized Laplacian
    D <- colSums(A)
    A <- A + (1/4)*mean(rowSums(A))/nrow(A)*matrix(1, nrow(A), nrow(A))
    L <- diag(rep(1,ncol(A))) -
      diag(D^(-1/2)) %*% A %*% diag(D^(-1/2))

    ## Absolute eigenvalue in order
    E <- order(-abs(eigen(L)$values))

    ## Go into eigenspace
    U <- eigen(L)$vectors[,E]
    U <- U[,c((ncol(U)-Q+1):ncol(U))]
    U <- U / rowSums(U^2)^(1/2)

    ## Applying the K-means in the eigenspace
    cl <- kmeans(U, Q, nstart = 10, iter.max = 30)$cluster
  } else {
    cl <- as.factor(rep(1,ncol(A)))
  }

  ## handing lonely souls
  cl.final[connected] <- cl
  cl.final[unconnected] <- which.min(rowsum(D,cl))
  ##   cl.final[unconnected] <- sample(levels(factor(cl)), length(unconnected), rep=TRUE, prob=table(cl)/sum(table(cl)))

  return(as.factor(cl.final))
}

SparseSpectralClustering_NAisMean <- function (A,Q) {

  ## basic handling of missing values
  if (anyNA(A)) A[is.na(A)] <- mean(A[!is.na(A)])

  ## handling lonely souls
  cl.final <- rep(NA, ncol(A))
  unconnected <- which(rowSums(A) == 0)
  connected <- setdiff(1:ncol(A), unconnected)

  A <- A[connected,connected]
  if (Q > 1) {

    ## Normalized Laplacian
    D <- colSums(A)
    A <- A + (1/4)*mean(rowSums(A))/nrow(A)*matrix(1, nrow(A), nrow(A))
    L <- diag(rep(1,ncol(A))) -
      diag(D^(-1/2)) %*% A %*% diag(D^(-1/2))

    ## Absolute eigenvalue in order
    E <- order(-abs(eigen(L)$values))

    ## Go into eigenspace
    U <- eigen(L)$vectors[,E]
    U <- U[,c((ncol(U)-Q+1):ncol(U))]
    U <- U / rowSums(U^2)^(1/2)

    ## Applying the K-means in the eigenspace
    cl <- kmeans(U, Q, nstart = 10, iter.max = 30)$cluster
  } else {
    cl <- as.factor(rep(1,ncol(A)))
  }

  ## handing lonely souls
  cl.final[connected] <- cl
  cl.final[unconnected] <- which.min(rowsum(D,cl))
  ##   cl.final[unconnected] <- sample(levels(factor(cl)), length(unconnected), rep=TRUE, prob=table(cl)/sum(table(cl)))

  return(as.factor(cl.final))
}


# library(missSBM)
# library(aricode)
# N <- 300
# Q <- 3
# alpha <- rep(1,Q)/Q
# pi <- diag(.25,Q) + .05
# directed <- FALSE
# mySBM <- simulateSBM(N, alpha, pi, directed)
#
# adjacencyMatrix <- mySBM$adjacencyMatrix
#
# samplingParameters <- .15
# sampling <- "dyad"
# sampledNetwork <- samplingSBM(adjacencyMatrix, sampling, samplingParameters)
#
#
# A <- sampledNetwork$adjacencyMatrix
#
# SpCl_NAiZero <- SpectralClustering_NAisZero(A,Q)
# SpCl_NAiMean <- SpectralClustering_NAisMean(A,Q)
# SparseSpCl_NAiZero <- SparseSpectralClustering_NAisZero(A,Q)
# SparseSpCl_NAiMean <- SparseSpectralClustering_NAisMean(A,Q)
#
# ARI(SpCl_NAiZero, mySBM$memberships)
# ARI(SpCl_NAiMean, mySBM$memberships)
# ARI(SparseSpCl_NAiZero, mySBM$memberships)
# ARI(SparseSpCl_NAiMean, mySBM$memberships)

