#' @import R6
#' @export
SBM <- # this virtual class is the mother of all subtypes of SBM (either sample or fit)
R6Class(classname = "SBM",
  ## fields for internal use (refering to mathematical notations)
  private = list(
    N        = NULL, # number of nodes
    alpha    = NULL, # vector of block parameters
    pi       = NULL, # matrix of connectivity
    directed = NULL  # directed or undirected network
  ),
  public = list(
    ## constructor
    initialize = function(directed = FALSE, nNodes=NA, mixtureParam=NA, connectParam=NA) {
      private$N        <- nNodes
      private$alpha    <- mixtureParam
      private$pi       <- connectParam
      private$directed <- directed
    }
  ),
  active = list(
    ## active binding to access fields outside the class
    nNodes    = function(value) {private$N}        , # number of nodes
    nBlocks   = function(value) {length(private$alpha)}, # number of blocks
    nDyads    = function(value) {ifelse(private$directed, self$nNodes*(self$nNodes - 1), self$nNodes*(self$nNodes - 1)/2)},
    direction = function(value) {if (private$directed) "directed" else "undirected"} , # directed network or not
    ## the following fields may change if a SBM is fitted
    mixtureParam = function(value) {                    # vector of block parameters (a.k.a. alpha)
      if (missing(value)) return(private$alpha) else private$alpha <- value
    },
    connectParam = function(value) { # matrix of connectivity (a.k.a. pi)
      if (missing(value)) return(private$pi) else private$pi <- value
    },
    df_mixtureParams = function(value) {self$nBlocks - 1},
    df_connectParams = function(value) {ifelse(private$directed, self$nBlocks^2, self$nBlocks*(self$nBlocks + 1)/2)}
  )
)

