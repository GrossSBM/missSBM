#' a sampled network
#' 
#' @field nNodes
#' @field nDyads
#' @field adjacencyMatrix
#' @field directed
#' @field missingDyads
#' @field observedDyads
#' @field samplingRate
#' 
#' @importFrom R6 R6Class
#' @export
sampledNetwork <- 
R6Class(classname = "sampledNetwork", 
  public = list(
    ## fields
    nNodes          = NULL, # number of nodes
    nDyads          = NULL, # number of dyads
    directed        = NULL, # directed network of not
    adjacencyMatrix = NULL, # adjacency matrix
    missingDyads    = NULL, # array indices of missing dyads
    observedDyads   = NULL, # array indices of observed dyads
    samplingRate    = NULL, # percentage of observed dyads
    samplingMatrix  = NULL, # matrix of observed and non-observed edges
  )
)  

sampledNetwork$set("public", "initialize",
function(adjacencyMatrix, directed) {
  ### TODO : check all arguments consistency
  if (!directed) stopifnot(isSymmetric(adjacencyMatrix))
  
  self$adjacencyMatrix <- adjacencyMatrix
  self$directed        <- directed
  self$nNodes          <- ncol(adjacencyMatrix)
  self$nDyads          <- ifelse(directed, 
                                 self$nNodes*(self$nNodes-1),
                                 self$nNodes*(self$nNodes-1)/2)
  # dyads are either defined on the whole matrix or on the lower triangular,
  # depending if the network is directed or not
  if (directed) {
    self$missingDyads  <- which( is.na(adjacencyMatrix))
    self$observedDyads <- which(!is.na(adjacencyMatrix))
  } else {
    self$missingDyads  <- which( is.na(adjacencyMatrix) & lower.tri(adjacencyMatrix) )
    self$observedDyads <- which(!is.na(adjacencyMatrix) & lower.tri(adjacencyMatrix) )
  }
  self$samplingRate    <- nrow(missingDyads)/nDyads
  
  
  self$samplingMatrix  <- matrix(0, self$nNodes, self$nNodes)
  self$samplingMatrix[observedDyads] <- 1
})

sampledNetwork$set("public", "plot",
plot = function(zlim=c(0,1), col=c("white", "midnightblue"), na.color='gray', outside.below.color='black', outside.above.color='white',...)
{
  z <- as.matrix(self$adjacencyMatrix)
  zstep <- (zlim[2] - zlim[1]) / length(col); # step in the color palette
  newz.below.outside <- zlim[1] - 2 * zstep # new z for values below zlim
  newz.above.outside <- zlim[2] + zstep # new z for values above zlim
  newz.na <- zlim[2] + 2 * zstep # new z for NA
    
  z[which(z<zlim[1])] <- newz.below.outside # we affect newz.below.outside
  z[which(z>zlim[2])] <- newz.above.outside # we affect newz.above.outside
  z[which(is.na(z>zlim[2]))] <- newz.na # same for newz.na
  
  zlim[1] <- zlim[1] - 2 * zstep # extend lower limit to include below value
  zlim[2] <- zlim[2] + 2 * zstep # extend top limit to include the two new values above and na
  
  col <- c(outside.below.color, col[1], col, outside.above.color, na.color) #correct by including col[1] at bottom of range
  
  par(mar=c(0.1,0.1,0.1,0.1))
  image(z[nrow(z):1,],  zlim=zlim, col=col, xaxt="n", yaxt="n") # we finally call image(...)
  box()
  par(mar=c(5.1,4.1,4.1,2.1))
})

sampledNetwork$set("public", "show",
function(){
  cat("\n=======================================")  
  cat("\n", ifelse(directed,"Directed","Undirected"), "network with", nNodes)  
  cat("\nPercentage of sampling rate:",samplingRate)
})
