# spectral_clustering <- function(A, vBlocks) {
#   n <- ncol(A)
#   A <- A %*% t(A) ## get second order paths between  node
#   ## handling lonely souls
#   unconnected <- which(rowSums(abs(A)) == 0)
#   connected   <- setdiff(1:n, unconnected)
#   A <- A[connected,connected]
#   ## Absolute Spectral clustering with Normalized weighted Laplacian
#   d <- 1/sqrt(rowSums(abs(A)))
#   L <- sweep(sweep(A, 1, d, "*"), 2, d, "*")
#   U <- eigs_sym(L, max(vBlocks))$vectors[, 1:max(vBlocks), drop = FALSE]
#   res <- future_lapply(vBlocks, function(k) {
#     cl <- rep(1L, n)
#     if (k != 1) {
#       Un <- U[, 1:k, drop = FALSE]
#       Un <- sweep(Un, 1, sqrt(rowSums(Un^2)), "/")
#       Un[is.nan(Un)] <- 0
#       cl_ <- as.integer(
#         kmeans_missSBM(Un, k)
#       )
#       ## handing lonely souls
#       cl[connected] <- cl_
#       cl[unconnected] <- which.max(rowsum(d, cl_))
#     }
#     cl
#   }, future.seed = TRUE, future.scheduling = structure(TRUE, ordering = "random"))
#   res
# }

#' @importFrom stats dist
kmeans_missSBM <- function(coordinates, k)
{

    n   <- nrow(coordinates)
    dim <- ncol(coordinates)

    if (k == 1) {
      classif <- rep(1,n)
    } else {
      centroids <- matrix(0,k,dim)

      dists <- as.matrix( dist( coordinates ))

      choosen <- which(dists==max(dists),arr.ind = TRUE)[1,c('row','col')]

      while(length(choosen)<k)
      {
        choosen<-c(choosen,which.max(apply(dists[choosen,],2,min)))
      }

      centroids <- coordinates[choosen,]

      classif <- kmeans_cpp(coordinates, centroids)
      classif <- as.vector(1 + classif)
    }
    classif
}
