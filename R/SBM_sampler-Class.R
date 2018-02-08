#' R6 Class definition of an sampler for a Stochastics Bloc Model
#'
#' this class is used to sample from an SBM, thus has some additional fields and methods related to the blocks and the adjancecy matrix (Z and X)
#'
#' @include SBM-Class.R
#' @export
SBM_sampler <-
R6Class(classname = "SBM_sampler",
  inherit = SBM,
  ## fields for internal use (refering to mathematical notations)
  private = list(
    Z     = NULL, # the sampled indicator of blocks
    X     = NULL, # the sampled adjacency matrix
    r_law = NULL  # random generation for the emission law of the edges
  ),
  public = list(
    initialize = function(family = "Bernoulli", directed = FALSE, nNodes=NA, mixtureParam=NA, connectParam=NA) {
      super$initialize(family, directed, nNodes, mixtureParam, connectParam)
      private$r_law <- switch(family,
        "Bernoulli" = function(n, prob) {rbinom(n, 1, prob)},
        "Poisson"   = function(n, prob) {rpois( n,    prob)}
        )
    },
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
  ## accessor to private fields
  active = list(
    blocks = function(value) {private$Z},
    memberships = function(value) {apply(private$Z, 1, which.max)},
    adjacencyMatrix = function(value) {private$X}
  )
)
