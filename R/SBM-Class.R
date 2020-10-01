#'  An R6 Class to represent a SBM
#'
#' @import R6
SBM <- # this virtual class is the mother of all subtypes of SBM (either sample or fit)
R6::R6Class(classname = "SBM",
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## fields for internal use (referring to mathematical notations)
  private = list(
    model   = NULL, # characters, the model name: distribution of the edges (bernoulli, poisson, gaussian)
    link    = NULL, # the link function (GLM-like)
    invlink = NULL, # the inverse link function (GLM-like)
    directed = NULL, # directed or undirected network
    N        = NULL, # number of nodes
    alpha    = NULL, # vector of block parameters
    pi       = NULL, # matrix of connectivity
    beta     = NULL, # vector of covariates parameters
    Y        = NULL, # the adjacency matrix
    X        = NULL  # list of covariates (list of dim[1] x dim[2] matrices)
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Initialize an [`SBM`] object
    #' @param directed logical indicating if the network is directed or not
    #' @param nbNodes The number of nodes
    #' @param blockProp the vector of mixture parameters
    #' @param connectParam the matrix of connectivity: inter/intra probabilities of connection when the network does not have covariates, or a logit scaled version of it.
    #' @param covarParam the vector of parameters associated with the covariates
    #' @param covarList A list with M entries (the M covariates). Each entry of the list must be an N x N matrix
    initialize = function(directed=FALSE, nbNodes=NA, blockProp=NA, connectParam=NA, covarParam=numeric(length(covarList)), covarList=list()) {

      stopifnot(all.equal(length(covarParam), length(covarList)))
      stopifnot(all(sapply(covarList, nrow) == nbNodes))
      stopifnot(all(sapply(covarList, ncol) == nbNodes))

      ## MODEL & PARAMETERS

      private$directed <- directed
      private$N        <- nbNodes
      private$alpha    <- blockProp
      private$pi       <- connectParam

      private$model   <- "bernoulli"
      private$X       <- covarList
      private$beta    <- covarParam
      private$link    <- .logit
      private$invlink <- .logistic
    },
    #' @description show method
    #' @param type character used to specify the type of SBM
    show = function(type = "Stochastic Block Model\n") {
      cat(type)
      cat("==================================================================\n")
      cat("Model", self$direction, "with",
          self$nbNodes,"nodes,",
          self$nbBlocks, "blocks and",
          ifelse(self$nbCovariates > 0, self$nbCovariates, "no"), "covariate(s).\n")
      cat("==================================================================\n")
      cat("* Useful fields \n")
      cat("  $nbNodes, $nbBlocks, $nbCovariates, $nbDyads\n", " $blockProp, $connectParam\n", "$covarParam, $covarArray \n")
    },
    #' @description User friendly print method
    print = function() { self$show() },
    #' @description basic matrix plot method for SBM object
    #' @param type character for the type of plot: either 'network' (true connection) or 'connectivity' (fitted connection). Default to 'network'.
    #' @importFrom corrplot corrplot
    plot = function(type = c("network", "connectivity")) {
      type <- match.arg(type)
      if (type == "network") {
        Z <- missSBM:::clustering_indicator(as.factor(self$memberships))
        colors <- matrix(-ncol(Z), ncol(Z), ncol(Z)); diag(colors) <- floor(ncol(Z)/2) + (1:ncol(Z)) # discriminate intra/inter cols
        colorMat <- Z %*% colors %*% t(Z)
        colorMap <- colorMat[order(self$memberships),order(self$memberships)]
        adjMatrix <- self$adjacencyMatrix[order(self$memberships), order(self$memberships)] * colorMap
        corrplot(adjMatrix, is.corr = F, tl.pos = "n", method = "color", cl.pos = "n", mar = c(0,0,1,0))
      }
      if (type == "connectivity") {
        corrplot(self$expectation[order(self$memberships), order(self$memberships)],
                 tl.pos = "n", method = "color", is.corr = FALSE, main = "")
      }
    }
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## ACTIVE BINDING
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field nbNodes The number of nodes
    nbNodes        = function(value) {private$N}, # number of nodes
    #' @field nbBlocks The number of blocks
    nbBlocks       = function(value) {length(self$blockProp)}, # number of blocks
    #' @field nbCovariates The number of covariates
    nbCovariates   = function(value) {length(private$X)}, # number of covariates
    #' @field nbDyads The number of dyads
    nbDyads        = function(value) {ifelse(private$directed, self$nbNodes*(self$nbNodes - 1), self$nbNodes*(self$nbNodes - 1)/2)},
    #' @field direction character indicating if the network is directed or not
    direction     = function(value) {if (private$directed) "directed" else "undirected"}, # directed network or not
    #' @field blockProp the vector of mixture parameters (block proportions)
    blockProp = function(value) {if (missing(value)) return(private$alpha) else private$alpha <- value},
    #' @field connectParam the matrix of connectivity: inter/intra probabilities of connection when the network does not have covariates, or a logit scaled version of it.
    connectParam     = function(value) {if (missing(value)) return(private$pi) else private$pi <- values},
    #' @field adjacencyMatrix  The adjacency matrix of the network
    adjacencyMatrix  = function(value) {if (missing(value)) return(private$Y) else private$Y <- value},
    #' @field covarParam the vector of parameters associated with the covariates
    covarParam       = function(value) {if (missing(value)) return(private$beta) else private$beta <- value},
    #' @field covarArray the array of covariates
    covarArray       = function(value) {simplify2array(private$X)},
    #' @field covarList list of matrices of covariates
    covarList    = function(value) {if (missing(value)) return(private$X) else private$X <- value},
    #' @field covarEffect effect of covariates
    covarEffect  = function(value) {if (self$nbCovariates > 0) return(roundProduct(self$covarArray, private$beta)) else return(numeric(0))},
    #' @field df_blockProps degrees of freedoms for the mixture parameters
    df_blockProps = function(value) {self$nbBlocks - 1},
    #' @field df_connectParams degrees of freedoms for the connectivity parameters
    df_connectParams = function(value) {ifelse(private$directed, self$nbBlocks^2, self$nbBlocks*(self$nbBlocks + 1)/2)},
    #' @field df_covarParams degrees of freedoms for the covariate parameters
    df_covarParams   = function(value) {self$nbCovariates}
  )
)

# SBM$set("public", "print"  , function() self$show())

## ----------------------------------------------------------------------
##
## PUBLIC S3 METHODS FOR SBM
##
## =========================================================================================

## Auxiliary functions to check the given class of an objet
is_SBM <- function(Robject) {inherits(Robject, "SBM")}

#' Extract model coefficients
#'
#' Extracts model coefficients from objects with class \code{\link[=SBM]{SBM}} and children
#'
#' @param object an R6 object inheriting from class SBM_fit (like SimpleSBM_fit or BipartiteSBM_fit)
#' @param type type of parameter that should be extracted. Either 'block' for \deqn{\pi}, 'connectivity' for \deqn{\theta},
#'  or "covariates" for \deqn{\beta}. Default is 'connectivity'.
#' @param ... additional parameters for S3 compatibility. Not used
#' @return vector or list of parameters.
#' @export
coef.SBM <- function(object, type = c( 'connectivity', 'block', 'covariates'), ...) {
  stopifnot(is_SBM(object))
  switch(match.arg(type),
         block        = object$blockProp,
         connectivity = object$connectParam,
         covariates   = object$covarParam)
}

#' @importFrom stats fitted
#' @export
fitted.SBM <- function(object, ...) {
  stopifnot(is_SBM(object))
  object$expectation
}

#' @export
summary.SBM <- function(object, ...) {
  stopifnot(is_SBM(object))
  object$show()
}
