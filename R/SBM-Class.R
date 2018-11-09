#' @import R6
SBM <- # this virtual class is the mother of all subtypes of SBM (either sample or fit)
R6::R6Class(classname = "SBM",
  ## fields for internal use (refering to mathematical notations)
  private = list(
    N        = NULL, # number of nodes
    Q        = NULL, # number of blocks
    M        = NULL, # number of covariates
    alpha    = NULL, # vector of block parameters
    pi       = NULL, # matrix of connectivity
    beta     = NULL, # vector of covariates parameters
    directed = NULL, # directed or undirected network
    phi      = NULL  # the similarity array (N x N x M)
  ),
  public = list(
    ## constructor
    initialize = function(directed = FALSE, nNodes=NA, mixtureParam=NA, connectParam=NA, covariates=NULL, covarParam=NULL) {
      private$directed <- directed
      private$N        <- nNodes
      private$Q        <- nrow(connectParam)
      private$M        <- ifelse(is.null(covariates), 0, length(covarParam))
      private$alpha    <- mixtureParam
      if (!is.null(covariates)) {
        stopifnot(dim(covariates)[1] == nNodes)
        stopifnot(dim(covariates)[2] == nNodes)
        stopifnot(dim(covariates)[3] == length(covarParam))
        private$phi  <- covariates
        private$beta <- covarParam
        private$pi   <- logit(connectParam)
      } else {
        private$pi   <- connectParam
      }
    }
  ),
  active = list(
    ## active binding to access fields outside the class
    nNodes      = function(value) {private$N}, # number of nodes
    nBlocks     = function(value) {private$Q}, # number of blocks
    nCovariates = function(value) {private$M}, # number of covariates
    nDyads      = function(value) {ifelse(private$directed, self$nNodes*(self$nNodes - 1), self$nNodes*(self$nNodes - 1)/2)},
    direction   = function(value) {if (private$directed) "directed" else "undirected"}, # directed network or not
    has_covariates = function(value) {ifelse(self$nCovariates > 0 , TRUE, FALSE)}, # with or without covariates
    ## the following fields may change if a SBM is fitted
    mixtureParam = function(value) { # vector of block parameters (a.k.a. alpha)
      if (missing(value)) return(private$alpha) else private$alpha <- value
    },
    connectParam = function(value) { # matrix of connectivity (a.k.a. pi)
      if (missing(value)) {
        if (self$has_covariates)
          return(logistic(private$pi))
        else
          return(private$pi)
      } else {
        private$pi <- value
      }
    },
    covarParam = function(value) { # vector of covariates (a.k.a. beta)
      if (missing(value)) return(private$beta) else private$beta <- value
    },
    covariates = function(value) {private$phi},
    df_mixtureParams = function(value) {self$nBlocks - 1},
    df_connectParams = function(value) {ifelse(private$directed, self$nBlocks^2, self$nBlocks*(self$nBlocks + 1)/2)},
    df_covarParams   = function(value) {self$nCovariates}
  )
)

## ----------------------------------------------------------------------
## PUBLIC METHODS FOR THE USERS
## ----------------------------------------------------------------------

SBM$set("public", "show",
function(model = "Stochastic Block Model\n") {
  cat(model)
  cat("==================================================================\n")
  cat("Model", self$direction, "with",
        self$nNodes,"nodes,",
        self$nBlocks, "blocks and",
        ifelse(self$has_covariates, self$nCovariates, "no"), "covariate(s).\n")
  cat("==================================================================\n")
  cat("* Useful fields \n")
  cat("  $nNodes, $nBlocks, $nCovariates, $nDyads\n", " $mixtureParam, $connectParam, $covarParam, $covariates\n")
})
SBM$set("public", "print", function() self$show())

