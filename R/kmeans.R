#' @importFrom stats dist
kmeans_missSBM <- function(coordinates, k)
{

    n   <- nrow(coordinates)
    dim <- ncol(coordinates)

    if (k == 1) {
      classif <- rep(1,n)
    } else {
      dists <- as.matrix( dist( coordinates ))

      choosen <- as.integer(which(dists==max(dists),arr.ind = TRUE)[1,c('row','col')])

      ## farthest-point seeding: incrementally track, for each point, its distance to the
      ## nearest already-chosen centroid, instead of recomputing the min over all chosen
      ## rows at every iteration (turns an O(k^2 N) loop into O(k N))
      min_dist <- pmin(dists[choosen[1], ], dists[choosen[2], ])
      while(length(choosen)<k)
      {
        next_point <- which.max(min_dist)
        choosen <- c(choosen, next_point)
        min_dist <- pmin(min_dist, dists[next_point, ])
      }

      centroids <- coordinates[choosen,]

      classif <- kmeans_cpp(coordinates, centroids)
      classif <- as.vector(1 + classif)
    }
    classif
}
