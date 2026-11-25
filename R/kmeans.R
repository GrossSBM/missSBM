kmeans_missSBM <- function(coordinates, k)
{

    n   <- nrow(coordinates)
    dim <- ncol(coordinates)

    if (k == 1) {
      classif <- rep(1,n)
    } else {
      ## farthest-point seeding, without the n x n matrix of distances: the farthest pair,
      ## then incrementally, for each point, its distance to the nearest already-chosen centroid
      ## (the same choices as with as.matrix(dist(coordinates)), ties included)
      coordinates <- as.matrix(coordinates)
      storage.mode(coordinates) <- "double"
      points   <- t(coordinates) # one point per column, contiguous in memory
      choosen  <- farthest_pair_cpp(points)
      min_dist <- pmin(dist_to_row_cpp(points, choosen[1]), dist_to_row_cpp(points, choosen[2]))
      while(length(choosen)<k)
      {
        next_point <- which.max(min_dist)
        choosen <- c(choosen, next_point)
        min_dist <- pmin(min_dist, dist_to_row_cpp(points, next_point))
      }

      centroids <- coordinates[choosen,]

      classif <- kmeans_cpp(coordinates, centroids)
      classif <- as.vector(1 + classif)
    }
    classif
}
