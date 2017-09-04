#' an SBM model
#'
#' @field nNodes          number of nodes
#' @field nBlocks         number of blocks
#' @field blockProportion vector of block proportion (a.k.a. alpha)
#' @field modelParameters vector of model parameters (a.k.a. theta)
#'
#' @importFrom R6 R6Class
#' @export
sampling <-
R6Class(classname = "sampling",
  public = list(
    ## fields
    nNodes         = NULL, # number of nodes
    missingParam   = NULL, # vector of missing parameters (a.k.a. alpha)
    completeLogLik = NULL, #
    samplingMatrix = NULL, #
    directed       = FALSE, #
    ## methods
    initialize = function(nNodes=NA, missingParam=NA, directed = FALSE) {
      self$nNodes         <- nNodes
      self$missingParam   <- missingParam
      self$directed       <- directed
    }
  )
)


sampling_doubleStandard <-
R6Class(classname = "sampling_doublestandard",
  inherit = sampling,
  public = list(
   initialize = function(nNodes, missingParam, directed = FALSE) {
     super$initialize(nNodes, missingParam, directed)
   },
   rSampling = function(adjMatrix) {
     self$samplingMatrix <- matrix(0, self$nNodes, self$nNodes)
     if(!self$directed){
       areOne  <- (adjMatrix == 1) & upper.tri(adjMatrix)
       areZero <- (adjMatrix == 0) & upper.tri(adjMatrix)
       
       self$samplingMatrix[areOne]  <- runif(sum(areOne))  < self$missingParam[2]
       self$samplingMatrix[areZero] <- runif(sum(areZero)) < self$missingParam[1]
       self$samplingMatrix <- t(self$samplingMatrix) | self$samplingMatrix
       diag(self$samplingMatrix) <- 1
     } else {
       areOne  <- (adjMatrix == 1)
       areZero <- (adjMatrix == 0)
       
       self$samplingMatrix[areOne]  <- runif(sum(areOne))  < self$missingParam[2]
       self$samplingMatrix[areZero] <- runif(sum(areZero)) < self$missingParam[1]
       diag(self$samplingMatrix) <- 1
     }


     sampAdjMatrix  <- adjMatrix
     sampAdjMatrix[which(self$samplingMatrix == 0)] <- NA
     return(list(samplingMatrix = self$samplingMatrix, sampAdjMatrix = sampAdjMatrix))
   }
   # samplingLogLik = function(adjMatrix) {
   #   obsEdges  <- which(!is.na(adjMatrix) & (upper.tri(adjMatrix) | lower.tri(adjMatrix)), arr.ind = TRUE)
   #   missEdges <- which(is.na(adjMatrix) , arr.ind = TRUE)
   #   return((sum(log(self$missingParam[2]) * adjMatrix[obsEdges] + log(self$missingParam[1]) * (1-adjMatrix[obsEdges])) +
   #             sum(log(1-self$missingParam[2]) * adjMatrix[missEdges] + log(1-self$missingParam[1]) * (1-adjMatrix[missEdges])))/2)
   # },
  )
)

# Undirected : 
mySBM <- sampling_doubleStandard$new(10, c(1/4, 1/2))
samp <- mySBM$rSampling(matrix(1,10,10))
# Directed : 
mySBM <- sampling_doubleStandard$new(10, c(1/4, 1/2), directed = TRUE)
samp <- mySBM$rSampling(matrix(1,10,10))

