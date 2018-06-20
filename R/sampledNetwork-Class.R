#' @import R6
#' @export
sampledNetwork <-
  R6::R6Class(classname = "sampledNetwork",
  ## FIELDS : encode network with missing edges
  private = list(
    X        = NULL, # adjacency matrix
    directed = NULL, # directed network of not
    D        = NULL, # list of potential dyads in the network
    nas      = NULL, # all NA in X
    D_obs    = NULL, # array indices of missing dyads
    D_miss   = NULL, # array indices of observed dyads
    R        = NULL, # matrix of observed and non-observed edges
    S        = NULL  # vector of observed and non-observed nodes
  ),
  ## Basically getters and setters for private fields
  active = list(
    ## percentage of observed dyads
    samplingRate = function(value) {length(private$D_obs)/self$nDyads},
    # number of nodes
    nNodes = function(value) {ncol(private$X)},
    # number of dyads
    nDyads = function(value) {#' An R package for adjusting Stochastic Block Models from networks data sampled under various missing data conditions
      #'
      #' The SBM package provides three categories of important functions:
      #' \code{\link{simulateSBM}}, \code{\link{samplingSBM}}, \code{\link{inferSBM}}.
      #'
      #' @details
      #' Package: missSBM
      #'
      #' Type:    Package
      #'
      #' Version: 0.1.0
      #'
      #' Date:    2017-12-01
      #'
      #' License: GPL-2 | GPL-3
      #'
      #' @docType package
      #' @author Timothée Tabouy, Pierre Barbillon, Julien Chiquet
      ifelse(private$directed, self$nNodes*(self$nNodes - 1), self$nNodes*(self$nNodes - 1)/2)
    },
    # direction
    is_directed = function(value) {private$directed},
    # adjacency matrix
    adjacencyMatrix = function(value) {private$X},
    # list of potential dyads in the network
    dyads           = function(value) {private$D},
    # array indices of missing dyads
    missingDyads    = function(value) {private$D_miss},
    # array indices of observed dyads
    observedDyads   = function(value) {private$D_obs},
    # matrix of observed and non-observed edges
    samplingMatrix  = function(value) {private$R},
    # vector of observed and non-observed nodes
    observedNodes   = function(value) {private$S},
    # boolean for NA entries in the adjacencyMatrix
    NAs             = function(value) {private$nas}
  ),
  ## Constructor
  public = list(
    initialize = function(adjacencyMatrix) {

      ## adjacency matrix
      stopifnot(is.matrix(adjacencyMatrix))
      if (isSymmetric(adjacencyMatrix)) private$directed <- FALSE else private$directed <- TRUE
      private$X  <- adjacencyMatrix

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
    plot = function(title = "Network sampling") {
      par(mfrow = c(1,2))
      image_NA(self$samplingMatrix , main = "sampling matrix")
      image_NA(self$adjacencyMatrix, main = "adjacency matrix")
      title(main = paste("\n",title,"with sampling rate:", signif(self$samplingRate,3)), outer = TRUE)
    }
  )
)
