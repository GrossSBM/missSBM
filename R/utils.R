SpectralClustering <- function(X, K) {

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

  return(as.factor(cl.final))
}
