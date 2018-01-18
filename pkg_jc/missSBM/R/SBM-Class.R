#' @import R6
#' @export
SBM <- # this virtual class is the mother of all subtypes of SBM (either sample or fit)
R6Class(classname = "SBM",
  ## fields for internal use (refering to mathematical notations)
  private = list(
    N        = NULL, # number of nodes
    Q        = NULL, # number of blocks
    alpha    = NULL, # vector of block parameters (a.k.a. alpha)
    pi       = NULL, # matrix of connectivity (a.k.a. pi)
    directed = NULL, # directed or undirected network
    family   = NULL  # emission law
  ),
  public = list(
    ## constructor
    initialize = function(family = "Bernoulli", directed = FALSE, nNodes=NA, mixtureParam=NA, connectParam=NA) {
      private$N        <- nNodes
      private$alpha    <- mixtureParam
      private$pi       <- connectParam
      private$Q        <- length(mixtureParam)
      private$family   <- family
      private$directed <- directed
    }
  ),
  active = list(
    ## active binding to access fields outside the class
    nNodes       = function(value) {private$N}        , # number of nodes
    nBlocks      = function(value) {private$Q}        , # number of blocks
    direction    = function(value) {if (private$directed) "directed" else "undirected"} , # directed network or not
    emissionLaw  = function(value) {private$family}   , # emission law
    ## the following fields may change if a SBM is fitted
    mixtureParam = function(value) {                    # vector of block parameters (a.k.a. alpha)
      if (missing(value)) return(private$alpha) else private$alpha <- value
    },
    connectParam = function(value) {                    # matrix of connectivity (a.k.a. pi)
      if (missing(value)) return(private$pi) else private$pi <- value
    },
    df_mixtureParams = function(value) {private$Q - 1},
    df_connectParams = function(value) {ifelse(private$directed, private$Q^2, private$Q*(private$Q + 1)/2)}
  )
)

