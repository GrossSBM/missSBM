#' @import R6
#' @export
sampledNetwork <-
R6::R6Class(classname = "sampledNetwork",
  ## FIELDS : encode network with missing edges
  private = list(
    X        = NULL, # adjacency matrix
    directed = NULL, # directed network of not
    D        = NULL, # list of potential dyads in the network
    D_obs    = NULL, # array indices of missing dyads
    D_miss   = NULL, # array indices of observed dyads
    R        = NULL, # matrix of observed and non-observed edges
    S        = NULL  # vector of observed and non-observed nodes
  ),
  ## Basically getters and setters for private fields
  active = list(
    ## percentage of observed dyads
    samplingRate = function(value) {
      if (!private$directed) {
        res <- (length(private$D_obs) - self$nNodes)/(2*self$nDyads)
      } else {
        res <- (length(private$D_obs) - self$nNodes)/self$nDyads
      }
      res
    },
    # number of nodes
    nNodes = function(value) {ncol(private$X)},
    # number of dyads
    nDyads = function(value) {
      ifelse(directed, self$nNodes*(self$nNodes - 1), self$nNodes*(self$nNodes - 1)/2)
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
    observedNodes   = function(value) {private$S}
  ),
  ## Constructor
  public = list(
    initialize = function(adjacencyMatrix) {

      ## adjacency matrix
      if (isSymmetric(adjacencyMatrix)) private$directed <- FALSE else private$directed <- TRUE
      private$X  <- adjacencyMatrix

      ## sets of observed / unobserved dyads
      if (private$directed) {
        ## remove diagonal( no loops)
        private$D_miss <- which( is.na(adjacencyMatrix) & (upper.tri(adjacencyMatrix) | lower.tri(adjacencyMatrix)) )
        private$D_obs  <- which(!is.na(adjacencyMatrix) & (upper.tri(adjacencyMatrix) | lower.tri(adjacencyMatrix)) )
      } else {
        private$D_miss <- which( is.na(adjacencyMatrix) & upper.tri(adjacencyMatrix))
        private$D_obs  <- which(!is.na(adjacencyMatrix) & upper.tri(adjacencyMatrix))
      }

      ## sets of observed / unobserved nodes
      S <- rep(FALSE, self$nNodes)
      S[!is.na(rowSums(adjacencyMatrix))] <- TRUE
      private$S <- S

      ## sampling matrix (indicating who is observed) : USELESS ??
      R <- matrix(0, self$nNodes, self$nNodes)
      R[private$D_obs] <- 1
      private$R <- R
    }
  )
)
