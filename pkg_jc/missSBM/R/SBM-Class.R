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
    family   = NULL, # emission law
    d_law    = NULL, # the density of the emission law of the edges
    r_law    = NULL  # random generation for the emission law of the edges
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
      private$r_law <- switch(family,
            "Bernoulli" = function(n, prob) {rbinom(n, 1, prob)},
            "Poisson"   = function(n, prob) {rpois( n,    prob)})
      private$d_law <- switch(family,
            "Bernoulli" = function(x, prob) {dbinom(x, 1, prob)},
            "Poisson"   = function(x, prob) {dpois( x,    prob)})
    }
  ),
  active = list(
    ## active binding to access fields outside the class
    nNodes       = function(value) {private$N}        , # number of nodes
    nBlocks      = function(value) {private$Q}        , # number of blocks
    isDirected   = function(value) {private$directed} , # directed network or not
    ## the following fields may change if a SBM is fitted
    mixtureParam = function(value) {                    # vector of block parameters (a.k.a. alpha)
      if (missing(value)) return(private$alpha) else private$alpha <- value
    },
    connectParam = function(value) {                    # matrix of connectivity (a.k.a. pi)
      if (missing(value)) return(private$pi) else private$pi <- value
    }
  )
)

