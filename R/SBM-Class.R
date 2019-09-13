#' @import R6
SBM <- # this virtual class is the mother of all subtypes of SBM (either sample or fit)
R6::R6Class(classname = "SBM",
  ## fields for internal use (refering to mathematical notations)
  private = list(
    directed = NULL, # directed or undirected network
    N        = NULL, # number of nodes
    Q        = NULL, # number of blocks
    M        = NULL, # number of covariates
    alpha    = NULL, # vector of block parameters
    pi       = NULL, # matrix of connectivity
    beta     = NULL, # vector of covariates parameters
    Y        = NULL, # the adjacency matrix
    X        = NULL  # the array of covariates (N x N x M)
  ),
  public = list(
    ## constructor
    initialize = function(directed=FALSE, nNodes=NA, mixtureParam=NA, connectParam=NA, covarParam=NULL, covarArray=NULL) {
      private$directed <- directed
      private$N        <- nNodes
      private$Q        <- length(mixtureParam)
      private$M        <- ifelse(is.null(covarArray), 0, length(covarParam))
      private$alpha    <- mixtureParam
      private$pi       <- connectParam
      if (!is.null(covarArray)) {
        stopifnot(all.equal(dim(covarArray), c(private$N, private$N, private$M)))
        private$X    <- covarArray
        private$beta <- covarParam
      }
    }
  ),
  active = list(
    ## active binding to access fields outside the class
    nNodes        = function(value) {private$N}, # number of nodes
    nBlocks       = function(value) {private$Q}, # number of blocks
    nCovariates   = function(value) {private$M}, # number of covariates
    nDyads        = function(value) {ifelse(private$directed, self$nNodes*(self$nNodes - 1), self$nNodes*(self$nNodes - 1)/2)},
    direction     = function(value) {if (private$directed) "directed" else "undirected"}, # directed network or not
    hasCovariates = function(value) {ifelse(self$nCovariates > 0 , TRUE, FALSE)}, # with or without covariates
    ## the following fields may change if a SBM is fitted
    mixtureParam = function(value) { # vector of block parameters (a.k.a. alpha)
      if (missing(value)) return(private$alpha) else private$alpha <- value
    },
    connectParam     = function(value) {if (missing(value)) return(private$pi) else private$pi <- values},
    adjacencyMatrix  = function(value) {if (missing(value)) return(private$Y) else private$Y <- value},
    covarParam       = function(value) {if (missing(value)) return(private$beta) else private$beta <- value},
    covarArray       = function(value) {private$X},
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
        ifelse(self$hasCovariates, self$nCovariates, "no"), "covariate(s).\n")
  cat("==================================================================\n")
  cat("* Useful fields \n")
  cat("  $nNodes, $nBlocks, $nCovariates, $nDyads\n", " $mixtureParam, $connectParam\n", "$covarParam, $covarArray \n")
})
SBM$set("public", "print"  , function() self$show())

SBM$set("public", "plot",
  function(type = c("network", "connectivity")) {
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
      corrplot(self$connectProb[order(self$memberships), order(self$memberships)],
               tl.pos = "n", method = "color", is.corr = FALSE, main = "")
    }
  }
)

## PUBLIC S3 METHODS FOR SBM
## =========================================================================================

## Auxiliary functions to check the given class of an objet
is_SBM <- function(Robject) {inherits(Robject, "SBM")}

#' @export
coef.SBM <- function(object, type = c("mixture", "connectivity", "covariates"), ...) {
  stopifnot(is_SBM(object))
  switch(match.arg(type),
         mixture      = object$mixtureParam,
         connectivity = object$connectParam,
         covariates   = object$covarParam)
}

#' @importFrom stats fitted
#' @export
fitted.SBM <- function(object, ...) {
  stopifnot(is_SBM(object))
  object$connectProb
}

#' @export
summary.SBM <- function(object, ...) {
  stopifnot(is_SBM(object))
  object$show()
}
