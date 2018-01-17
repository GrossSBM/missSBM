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

#' @import R6
#' @export
SBM_sample <- # this class is used to sample from an SBM, thus has some additional fields
  ## and methods related to the blocks and the adjancecy matrix (Z and X)
R6Class(classname = "SBM_sample",
  inherit = SBM,
  ## fields for internal use (refering to mathematical notations)
  private = list(
    Z        = NULL, # the sampled indicator of blocks
    X        = NULL  # the sampled adjacency matrix
  ),
  public = list(
    ## constructor is the same as the above, so no need to specify initialize
    ## a method to generate a vector of clusters indicators
    rBlocks = function() {
      private$Z <- t(rmultinom(private$N, size = 1, prob = private$alpha))
    },
    ## a method to sample an adjacency matrix for the current SBM
    rAdjMatrix = function() {
      X <- matrix(private$r_law(private$N^2, private$Z %*% private$pi %*% t(private$Z)), private$N)
      if (!private$directed) X <- X * lower.tri(X) + t(X * lower.tri(X))
      diag(X) <- 0
      private$X <- X
    }
  ),
  active = list(
    blocks       = function(value) {                # indicator of blocks
      if (missing(value)) return(private$Z) else private$Z <- value
      }    ,
    clusters     = function(value) {                # vector of clusters
      if (!is.null(private$Z)) apply(private$Z, 1, which.max) else NA
    }    ,
    adjacencyMatrix = function(value) {
      if (missing(value)) return(private$X) else private$X <- value
    }
  )
)

