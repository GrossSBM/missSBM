#' @import R6
#' @export
sampling_model <-
R6Class(classname = "sampling",
  ## fields
  private = list(
    name  = NULL, # type of sampling
    psi   = NULL  # vector of missing parameters (a.k.a. psi)
  ),
  public = list(
    ## methods
    initialize = function(type = NA, parameter = NA) {
      private$name <- type
      private$psi  <- parameters
    }
  ),
  active = list(
    type = function(value) {
      if (missing(value)) return(private$psi) else private$psi <- value
    },
    parameters = function(value) {
      if (missing(value)) return(private$psi) else private$psi <- value
    }
  )
)

sampling_model_doubleStandard <-
R6Class(classname = "sampling_model_doubleStandard",
  inherit = sampling_model,
  public = list(
    rSampling = function(adjMatrix) {
      samplingMatrix <- matrix(0, self$nNodes, self$nNodes)
      if (!self$directed) {
        areOne  <- (adjMatrix == 1) & upper.tri(adjMatrix)
        areZero <- (adjMatrix == 0) & upper.tri(adjMatrix)

        samplingMatrix[areOne]  <- runif(sum(areOne))  < self$missingParam[2]
        samplingMatrix[areZero] <- runif(sum(areZero)) < self$missingParam[1]
        samplingMatrix <- t(samplingMatrix) | samplingMatrix
        diag(samplingMatrix) <- 1
      } else {
        areOne  <- (adjMatrix == 1)
        areZero <- (adjMatrix == 0)

        samplingMatrix[areOne]  <- runif(sum(areOne))  < self$missingParam[2]
        samplingMatrix[areZero] <- runif(sum(areZero)) < self$missingParam[1]
        diag(samplingMatrix) <- 1
      }

      sampAdjMatrix  <- adjMatrix
      sampAdjMatrix[which(samplingMatrix == 0)] <- NA
      return(sampledNetwork$new(sampAdjMatrix, self$directed))
    }
  )
)

sampling_class <-
R6Class(classname = "sampling_class",
  inherit = sampling,
  public = list(
    rSampling = function(adjMatrix, blockVarParam) {
      samplingMatrix <- matrix(0, self$nNodes, self$nNodes)
      sampProb <- runif(self$nNodes) < self$missingParam[apply(blockVarParam, 1, which.max)]
      obsNodes <- which(runif(self$nNodes) < sampProb)

      samplingMatrix <- matrix(0,self$nNodes,self$nNodes) ; samplingMatrix[obsNodes,] <- 1
      diag(samplingMatrix) <- 1
      samplingMatrix <- (t(samplingMatrix) | samplingMatrix)*1

      sampAdjMatrix  <- adjMatrix
      sampAdjMatrix[which(samplingMatrix == 0)] <- NA
      return(sampledNetwork$new(sampAdjMatrix, self$directed))
    }
  )
)

sampling_starDegree <-
R6Class(classname = "sampling_starDegree",
  inherit = sampling,
  public = list(
    rSampling = function(adjMatrix) {
      samplingMatrix      <- matrix(0, self$nNodes, self$nNodes)
      sampProb            <- self$missingParam[1]+self$missingParam[2]*rowSums(adjMatrix)
      samprob             <- 1/(1+exp(-sampProb))
      obsNodes            <- which(runif(self$nNodes) < sampProb)

      samplingMatrix <- matrix(0,self$nNodes,self$nNodes) ; samplingMatrix[obsNodes,] <- 1
      diag(samplingMatrix) <- 1
      samplingMatrix <- (t(samplingMatrix) | samplingMatrix)*1

      sampAdjMatrix  <- adjMatrix ; sampAdjMatrix[which(samplingMatrix == 0)] <- NA
      return(sampledNetwork$new(sampAdjMatrix, self$directed))
    }
  ),
  private = list(
    g = function(x){
      return(-(1/(1+exp(-x)) - 0.5)/(0.5*x))
    }
  )
)

sampling_randomPairMAR <-
R6Class(classname = "sampling_randomPairMAR",
  inherit = sampling,
  public = list(
    rSampling = function(adjMatrix) {
      samplingMatrix <- matrix(0, self$nNodes, self$nNodes)
      if(!self$directed){
        edgeSamp <- sample(which(lower.tri(adjMatrix)), floor((self$nNodes*(self$nNodes-1)/2)*self$missingParam))
      } else {
        edgeSamp <- sample(which(lower.tri(adjMatrix) | upper.tri(adjMatrix)), floor((self$nNodes^2 - self$nNodes)*self$missingParam))
      }

      samplingMatrix <- matrix(0,self$nNodes,self$nNodes)
      samplingMatrix[edgeSamp] <- 1
      if(!self$directed){ samplingMatrix <- t(samplingMatrix) | samplingMatrix }
      diag(samplingMatrix) <- 1

      sampAdjMatrix  <- adjMatrix ; sampAdjMatrix[which(samplingMatrix == 0)] <- NA
      return(sampledNetwork$new(sampAdjMatrix, self$directed))
    }
  )
)

sampling_randomNodesMAR <-
R6Class(classname = "sampling_randomNodesMAR",
  inherit = sampling,
    public = list(
      rSampling = function(adjMatrix) {
        samplingMatrix <- matrix(0, self$nNodes, self$nNodes)
        obsNodes       <- sample(1:self$nNodes, floor((self$nNodes)*self$missingParam))

        samplingMatrix <- matrix(0,self$nNodes,self$nNodes) ; samplingMatrix[obsNodes,] <- 1
        diag(samplingMatrix) <- 1
        samplingMatrix <- (t(samplingMatrix) | samplingMatrix)*1

        sampAdjMatrix  <- adjMatrix ; sampAdjMatrix[which(samplingMatrix == 0)] <- NA
        return(sampledNetwork$new(sampAdjMatrix, self$directed))
      }
    )
)

sampling_snowball <-
R6Class(classname = "sampling_snowball",
  inherit = sampling,
  public = list(
    rSampling = function(adjMatrix) {
      samplingMatrix <- matrix(0, self$nNodes, self$nNodes)
      # browser()

      # First step :
      obsNodes       <- which(runif(self$nNodes) < self$missingParam)
      samplingMatrix[obsNodes,] <- 1
      samplingMatrix[,obsNodes] <- 1

      # Second step :
      obsNeighbours <- apply(adjMatrix[obsNodes, ], 1, function(x) which(x != 0))
      obsNeighbours <- unique(unlist(obsNeighbours))
      samplingMatrix[obsNeighbours, ] <- samplingMatrix[, obsNeighbours] <- 1

      diag(samplingMatrix) <- 1

      sampAdjMatrix  <- adjMatrix ; sampAdjMatrix[which(samplingMatrix == 0)] <- NA
      return(sampledNetwork$new(sampAdjMatrix, self$directed))
    }
  )
)

