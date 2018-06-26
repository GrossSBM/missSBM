#' @import R6
#' @export
SBM <- # this virtual class is the mother of all subtypes of SBM (either sample or fit)
R6::R6Class(classname = "SBM",
  ## fields for internal use (refering to mathematical notations)
  private = list(
    N        = NULL, # number of nodes
    Q        = NULL, # number of blocks
    D        = NULL, # the number of covariates
    alpha    = NULL, # vector of block parameters
    pi       = NULL, # matrix of connectivity
    beta     = NULL, # vector of covariates parameters
    X        = NULL, # the matrix of covariates (N x D)
    directed = NULL  # directed or undirected network
  ),
  public = list(
    ## constructor
    initialize = function(directed = FALSE, nNodes=NA, mixtureParam=NA, connectParam=NA, covariates=NULL, covarParam=NULL) {
      private$N        <- nNodes
      private$Q        <- nrow(as.matrix(connectParam))
      private$D        <- ifelse(is.null(covariates), 0, ncol(X))
      private$alpha    <- mixtureParam
      private$pi       <- connectParam
      private$beta     <- covarParam
      private$directed <- directed
      private$X        <- covariates
    }
  ),
  active = list(
    ## active binding to access fields outside the class
    nNodes      = function(value) {private$N}, # number of nodes
    nBlocks     = function(value) {private$Q}, # number of blocks
    nCovariates = function(value) {private$D}, # number of covariates
    nDyads      = function(value) {ifelse(private$directed, self$nNodes*(self$nNodes - 1), self$nNodes*(self$nNodes - 1)/2)},
    direction   = function(value) {if (private$directed) "directed" else "undirected"} , # directed network or not
    has_covariates = function(value) {ifelse(private$D > 0, TRUE, FALSE)}, # covariates or not
    ## the following fields may change if a SBM is fitted
    mixtureParam = function(value) { # vector of block parameters (a.k.a. alpha)
      if (missing(value)) return(private$alpha) else private$alpha <- value
    },
    connectParam = function(value) { # matrix of connectivity (a.k.a. pi)
      if (missing(value)) return(private$pi) else private$pi <- value
    },
    covarParam = function(value) { # vector of covariates (a.k.a. beta)
      if (missing(value)) return(private$beta) else private$beta <- value
    },
    df_mixtureParams = function(value) {self$nBlocks - 1},
    df_connectParams = function(value) {ifelse(private$directed, self$nBlocks^2, self$nBlocks*(self$nBlocks + 1)/2)}
  )
)

#' @import R6
#' @export
SBM_cov_model1 <- # virtual class,  mother of all subtypes of SBM with model 1 of covariates (either sample or fit)
R6::R6Class(classname = "SBM_cov_model1",
  inherit = SBM,
  active = list(
    df_covarParams = function(value) {self$nCovariates * (self$nBlocks - 1)}
  )
)

#' @import R6
#' @export
SBM_cov_model2 <- # virtual class, mother of all subtypes of SBM with model 2 of covariates (either sample or fit)
R6::R6Class(classname = "SBM_cov_model2",
  inherit = SBM,
  ## fields for internal use (refering to mathematical notations)
  private = list(
    phi = NULL,  # a similarity function N x N -> D
    sim = NULL,  # a similarity function N x N -> D
  ),
  public = list(
    ## constructor
    initialize = function(directed = FALSE, nNodes=NA, mixtureParam=NA, connectParam=NA, covariates=NA, covarParam=NA, similarity=NA) {
      super$initialize(directed, nNodes, mixtureParam, connectParam, covariates, covarParam)
      private$phi  <- similarity
    }
  ),
  active = list(
      df_covarParams = function(value) {self$nCovariates}
    )
  )

