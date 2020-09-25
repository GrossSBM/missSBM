#' An R6 Class used for internal representation of sampled network data
#'
#' All fields of this class are only accessible for reading.
#'
#' @importFrom R6 R6Class
sampledNetwork <-
  R6::R6Class(classname = "sampledNetwork",
  ## FIELDS : encode network with missing edges
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    Y        = NULL, # adjacency matrix
    X        = NULL, # the covariates matrix
    phi      = NULL, # the covariates array
    directed = NULL, # directed network of not
    D        = NULL, # list of potential dyads in the network
    nas      = NULL, # all NA in Y
    D_obs    = NULL, # array indices of missing dyads
    D_miss   = NULL, # array indices of observed dyads
    R        = NULL, # matrix of observed and non-observed edges
    S        = NULL  # vector of observed and non-observed nodes
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## ACTIVE BINDING
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field samplingRate The percentage of observed dyads
    samplingRate = function(value) {length(private$D_obs)/self$nDyads},
    #' @field nNodes The number of nodes
    nNodes = function(value) {ncol(private$Y)},
    #' @field nDyads The number of dyads
    nDyads = function(value) {
      ifelse(private$directed, self$nNodes*(self$nNodes - 1), self$nNodes*(self$nNodes - 1)/2)
    },
    #' @field is_directed logical indicating if the network is directed or not
    is_directed = function(value) {private$directed},
    #' @field adjacencyMatrix  The adjacency matrix of the network
    adjacencyMatrix = function(value) {private$Y},
    #' @field covarArray the array of covariates
    covarArray = function(value) {private$phi},
    #' @field covarMatrix the matrix of covariates
    covarMatrix = function(value) {if (missing(value)) return(private$X) else  private$X <- value},
    #' @field dyads a list of potential dyads in the network
    dyads           = function(value) {private$D},
    #' @field missingDyads array indices of missing dyads
    missingDyads    = function(value) {private$D_miss},
    #' @field observedDyads array indices of observed dyads
    observedDyads   = function(value) {private$D_obs},
    #' @field samplingMatrix matrix of observed and non-observed edges
    samplingMatrix  = function(value) {private$R},
    #' @field observedNodes a vector of observed and non-observed nodes
    observedNodes   = function(value) {private$S},
    #' @field NAs boolean for NA entries in the adjacencyMatrix
    NAs             = function(value) {private$nas}
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description constructor
    #' @param adjacencyMatrix The adjacency matrix of the network
    #' @param covariates A list with M entries (the M covariates), each of whom being either a size-N vector or N x N matrix.
    #' @param similarity An R x R -> R function to compute similarities between node covariates. Default is \code{l1_similarity}, that is, -abs(x-y).
    initialize = function(adjacencyMatrix, covariates = NULL, similarity = l1_similarity) {

      ## adjacency matrix
      stopifnot(is.matrix(adjacencyMatrix))
      ## Only binary graph supported
      stopifnot(all.equal(unique(as.numeric(adjacencyMatrix[!is.na(adjacencyMatrix)])), c(0,1)))

      if (isSymmetric(adjacencyMatrix)) private$directed <- FALSE else private$directed <- TRUE
      private$Y  <- adjacencyMatrix

      ## covariates
      covar <- format_covariates(covariates, similarity)
      private$X   <- covar$Matrix
      private$phi <- covar$Array

      ## sets of observed / unobserved dyads
      private$nas <- is.na(adjacencyMatrix)
      if (private$directed) {
        ## remove diagonal (no loops)
        private$D <- which(upper.tri(adjacencyMatrix) | lower.tri(adjacencyMatrix))
      } else {
        private$D <- which(upper.tri(adjacencyMatrix))
      }
      private$D_miss <- intersect(which( private$nas), private$D )
      private$D_obs  <- intersect(which(!private$nas), private$D )

      ## sets of observed / unobserved nodes
      S <- rep(FALSE, self$nNodes)
      S[!is.na(rowSums(adjacencyMatrix))] <- TRUE
      private$S <- S

      ## sampling matrix (indicating who is observed) : USELESS ??
      R <- matrix(0, self$nNodes, self$nNodes)
      R[private$D_obs] <- 1
      if (!private$directed)  R <- t(R) | R
      private$R <- R
    },
    #' @description plot method for sampledNetwork
    #' @param clustering an optional vector of clustering memberships, default to \code{NULL}.
    #' @param main a character for the title of the plot
    #' @importFrom corrplot corrplot
    plot = function(clustering = NULL, main = paste("Network with sampling rate:", signif(self$samplingRate,3))) {
      if (is.null(clustering)) {
        adjMatrix <- self$adjacencyMatrix
      } else {
        Z <- missSBM:::clustering_indicator(as.factor(clustering))
        colors <- matrix(- ncol(Z), ncol(Z), ncol(Z)); diag(colors) <- floor(ncol(Z)/2) + (1:ncol(Z)) # discriminate intra/inter cols
        colorMat <- Z %*% colors %*% t(Z)
        colorMap <- colorMat[order(clustering),order(clustering)]
        adjMatrix <- self$adjacencyMatrix[order(clustering), order(clustering)] * colorMap
      }
      corrplot(adjMatrix, is.corr = F, tl.pos = "n", method = "color", cl.pos = "n", na.label.col = "grey", main = main, mar = c(0,0,1,0))
    },
    #' @description show method
    show = function() {
      cat("Sampled Network\n")
      cat("==================================================================\n")
      cat("Structure for storing a sampled network in missSBM\n")
      cat("==================================================================\n")
      cat("* Useful fields \n")
      cat("  $nNodes, $nDyads, $is_directed\n", "  $adjacencyMatrix, $covarMatrix, $covarArray\n",
          "  $dyads, $missingDyads, $observedDyads, $observedNodes\n",  "  $samplingRate, $samplingMatrix, $NAs\n")
      cat("* Useful method: plot(), summary() , print()  \n")
    },
    #' @description User friendly print method
    print = function() { self$show() }
  )
)



## PUBLIC S3 METHODS FOR sampledNetwork
## =========================================================================================

## Auxiliary functions to check the given class of an objet
is_sampledNetwork <- function(Robject) {inherits(Robject, "sampledNetwork")}
