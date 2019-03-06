#' @import R6
sampledNetwork <-
  R6::R6Class(classname = "sampledNetwork",
  ## FIELDS : encode network with missing edges
  private = list(
    Y        = NULL, # adjacency matrix
    X        = NULL, # M x N covariates matrix
    # s        = NULL, # the similarity function (N x N -> M)
    # phi      = NULL, # N x N x M covariates array
    directed = NULL, # directed network of not
    D        = NULL, # list of potential dyads in the network
    nas      = NULL, # all NA in Y
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
    nNodes = function(value) {ncol(private$Y)},
    # number of dyads
    nDyads = function(value) {
      ifelse(private$directed, self$nNodes*(self$nNodes - 1), self$nNodes*(self$nNodes - 1)/2)
    },
    # direction
    is_directed = function(value) {private$directed},
    # adjacency matrix
    adjMatrix = function(value) {private$Y},
    # covariates matrix
    covarMatrix = function(value) {if (missing(value)) return(private$X) else  private$X <- value},
    # covarArray  = function(value) {private$phi},
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
    initialize = function(adjMatrix, covarMatrix = NULL) {

      ## adjacency matrix
      stopifnot(is.matrix(adjMatrix))
      if (isSymmetric(adjMatrix)) private$directed <- FALSE else private$directed <- TRUE
      private$Y  <- adjMatrix

      ## array of covariates
      if (!is.null(covarMatrix)) {
        private$X <- covarMatrix
        # N <- nrow(covariatesMatrix)
        # M <- ncol(covariatesMatrix)
        # covariates <- array(dim = c(N, N, M))
        # for (i in 1:N)
        #   for (j in 1:N)
        #     covariates[i,j,] <- -abs(covariatesMatrix[i, ] - covariatesMatrix[j, ])
        # private$phi <- covariates
      }

      ## sets of observed / unobserved dyads
      private$nas <- is.na(adjMatrix)
      if (private$directed) {
        ## remove diagonal (no loops)
        private$D <- which(upper.tri(adjMatrix) | lower.tri(adjMatrix))
      } else {
        private$D <- which(upper.tri(adjMatrix))
      }
      private$D_miss <- intersect(which( private$nas), private$D )
      private$D_obs  <- intersect(which(!private$nas), private$D )

      ## sets of observed / unobserved nodes
      S <- rep(FALSE, self$nNodes)
      S[!is.na(rowSums(adjMatrix))] <- TRUE
      private$S <- S

      ## sampling matrix (indicating who is observed) : USELESS ??
      R <- matrix(0, self$nNodes, self$nNodes)
      R[private$D_obs] <- 1
      if (!private$directed)  R <- t(R) | R
      private$R <- R
    }
  )
)

#' @export
sampledNetwork$set("public", "plot",
function(title = "Network sampling") {
  par(mfrow = c(1,2))
  image_NA(self$samplingMatrix , main = "sampling matrix")
  image_NA(self$adjMatrix, main = "adjacency matrix")
  title(main = paste("\n",title,"with sampling rate:", signif(self$samplingRate,3)), outer = TRUE)
  par(mfrow=c(1,1))
})

#' @export
sampledNetwork$set("public", "show",
function(model = "Sampled Network\n") {
  cat(model)
  cat("==================================================================\n")
  cat("Structure for storing a sampled network in missSBM.\n")
  cat("==================================================================\n")
  cat("* Useful fields \n")
  cat("  $nNodes, $nDyads, $is_directed\n", "  $adjMatrix, $covarMatrix\n",
      "  $dyads, $missingDyads, $observedDyads, $observedNodes\n",  "  $samplingRate, $samplingMatrix, $NAs\n")
  cat("* Useful method: plot() \n")
})
sampledNetwork$set("public", "print", function() self$show())
