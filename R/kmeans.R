#' @importFrom stat dist
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
