#' R6 Class definition of a sampler for a Stochastic Block Model
#'
#' this class is used to sample from an SBM, thus has some additional fields and methods related to the blocks and the adjacency matrix (Z and Y)
#'
#' @include SBM-Class.R
SBM_sampler <-
  R6::R6Class(classname = "SBM_sampler",
  inherit = SBM,
  ## fields for internal use (refering to mathematical notations)
  private = list(
    Z     = NULL, # the sampled indicator of blocks
    Y     = NULL  # the sampled adjacency matrix
  ),
  public = list(
    initialize = function(directed = FALSE, nNodes=NA, mixtureParam=NA, connectParam=NA, covarParam=NULL, covarArray=NULL) {
      super$initialize(directed, nNodes, mixtureParam, connectParam, covarParam, covarArray)
    },
    ## constructor is the same as the above, so no need to specify initialize
    ## a method to generate a vector of clusters indicators
    rBlocks = function() {
      private$Z <- t(rmultinom(private$N, size = 1, prob = private$alpha))
    },
    ## a method to sample an adjacency matrix for the current SBM
    rAdjMatrix = function() {
      Y <- matrix(rbinom(private$N^2, 1, self$connectProb), private$N)
      if (!private$directed) Y <- Y * lower.tri(Y) + t(Y * lower.tri(Y))
      private$Y <- Y
    }
  ),
  ## accessor to private fields
  active = list(
    blocks = function(value) {private$Z},
    memberships = function(value) {apply(private$Z, 1, which.max)},
    adjMatrix   = function(value) {private$Y},
    connectProb = function(value) {
      PI <- private$Z %*% private$pi %*% t(private$Z)
      if (self$hasCovariates) {
        PI <- logistic(PI + roundProduct(private$phi, private$beta))
      }
      PI
    }
  )
)

## ----------------------------------------------------------------------
## PUBLIC METHODS FOR THE USERS
## ----------------------------------------------------------------------

SBM_sampler$set("public", "show",
function(model = "Sampler for Stochastic Block Model\n") {
  super$show(model)
  cat("* Sampling methods \n")
  cat("  $rBlocks(), $rAdjmatrix()\n")
  cat("* Additional fields \n")
  cat("  $blocks, $memberships, $adjMatrix, $connectProb\n")
  cat("* S3 methods \n")
  cat("  $plot\n")
})

SBM_sampler$set("public", "print", function() self$show())

SBM_sampler$set("public", "plot",
  function(type = c("network", "connectivity")) {
    type <- match.arg(type)
    if (type == "network") {
      image_NA(self$adjMatrix[order(self$memberships), order(self$memberships)])
    }
    if (type == "connectivity") {
      plot(
        graph_from_adjacency_matrix(
          private$pi,
          mode = ifelse(private$directed, "directed", "undirected"),
          weighted = TRUE, diag = TRUE
        ), main = "SBM connectivity matrix"
      )
    }
  }
)
