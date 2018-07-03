#' R6 Class definition of a sampler for a Stochastics Bloc Model
#'
#' this class is used to sample from an SBM, thus has some additional fields and methods related to the blocks and the adjancecy matrix (Z and Y)
#'
#' @include SBM-Class.R
#' @export
SBM_sampler <-
  R6::R6Class(classname = "SBM_sampler",
  inherit = SBM,
  ## fields for internal use (refering to mathematical notations)
  private = list(
    Z     = NULL, # the sampled indicator of blocks
    Y     = NULL  # the sampled adjacency matrix
  ),
  public = list(
    initialize = function(directed = FALSE, nNodes=NA, mixtureParam=NA, connectParam=NA, covariates=NULL, covarParam=NULL, covarSimilarity=NULL) {
      super$initialize(directed, nNodes, mixtureParam, connectParam, covariates, covarParam, covarSimilarity)
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
    adjacencyMatrix = function(value) {private$Y},
    connectProb = function(value) {
      PI <- private$Z %*% private$pi %*% t(private$Z)
      if (self$has_covariates) {
        PI <- logistic(PI + roundProduct(private$cov, private$beta))
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
  cat("  $blocks, $memberships, $adjacencyMatrix, $connectProb\n")

})

SBM_sampler$set("public", "print", function() self$show())

SBM_sampler$set("public", "plot",
  function(type = c("network", "connectivity")) {
    type <- match.arg(type)
    if (type == "network") {
      g <- gg_image_NA(self$adjacencyMatrix, self$memberships)
      print(g)
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
