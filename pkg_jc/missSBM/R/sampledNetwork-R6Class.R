#' @import R6
#' @export
sampledNetwork <-
  R6::R6Class(classname = "sampledNetwork",
  public = list(
    adjacencyMatrix = NULL, # adjacency matrix
    directed        = NULL, # directed network of not
    dyads           = NULL, #
    missingDyads    = NULL, # array indices of missing dyads
    observedDyads   = NULL, # array indices of observed dyads
    samplingMatrix  = NULL, # matrix of observed and non-observed edges
    samplingVector  = NULL  # vector of observed and non-observed nodes
  ),
  active = list(
    ## percentage of observed dyads
    samplingRate = function(value) {
      if (!self$directed) {
        res <- (length(self$observedDyads) - self$nNodes)/(2*self$nDyads)
      } else {
        res <- (length(self$observedDyads) - self$nNodes)/self$nDyads
      }
      res
    },
    # number of nodes
    nNodes = function(value) {ncol(self$adjacencyMatrix)},
    # number of dyads
    nDyads = function(value) {
      ifelse(directed, self$nNodes*(self$nNodes - 1), self$nNodes*(self$nNodes - 1)/2)
    }
  )
)

sampledNetwork$set("public", "initialize",
function(adjacencyMatrix) {
  if (isSymmetric(adjacencyMatrix)) self$directed <- FALSE else self$directed <- TRUE

  self$adjacencyMatrix <- adjacencyMatrix
  self$nNodes          <- ncol(adjacencyMatrix)
  self$missingDyads  <- which( is.na(adjacencyMatrix))
  self$observedDyads <- which(!is.na(adjacencyMatrix))

  self$samplingMatrix  <- matrix(0, self$nNodes, self$nNodes)
  self$samplingMatrix[self$observedDyads] <- 1

  self$samplingVector <- rep(0, self$nNodes)
  self$samplingVector[which(!is.na(rowSums(adjacencyMatrix)))] <- 1

})

