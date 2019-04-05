#' An R6 Class to represent a sampler for a SBM
#'
#' The function \code{\link{simulate}} produces an instance of an object with class \code{SBM_sampler}.
#'
#' All fields of this class are only accessible for reading. This class comes with a basic plot and print methods
#'
#' @field samplingRate percentage of observed dyads
#' @field nNodes number of nodes
#' @field nDyads number of dyads
#' @field is_directed direction
#' @field adjacencyMatrix adjacency matrix (with NA)
#' @field covarMatrix the matrix of covariates (if applicable)
#' @field covarArray the array of covariates (if applicable)
#' @field dyads list of potential dyads in the network
#' @field missingDyads array indices of missing dyads
#' @field observedDyads array indices of observed dyads
#' @field samplingMatrix matrix of observed and non-observed edges
#' @field observedNodes vector of observed and non-observed nodes
#' @field NAs boolean for NA entries in the adjacencyMatrix
#'
#' @importFrom R6 R6Class
#' @examples
#' ## SBM parameters
#' directed <- FALSE
#' N <- 300 # number of nodes
#' Q <- 3   # number of clusters
#' alpha <- rep(1,Q)/Q     # mixture parameter
#' pi <- diag(.45,Q) + .05 # connectivity matrix
#'
#' ## simulate a SBM without covariates
#' sbm <- missSBM::simulate(N, alpha, pi, directed)
#'
#' ## Sample network data
#' sampled_network <-
#'      missSBM::sample(
#'        adjacencyMatrix = sbm$adjacencyMatrix,
#'        sampling        = "double-standard",
#'        parameters      = c(0.4, 0.8)
#'      )
#'
#' print(sampled_network)
#' plot(sampled_network)
#'@export
sampledNetwork <-
  R6::R6Class(classname = "sampledNetwork",
  ## FIELDS : encode network with missing edges
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
    adjacencyMatrix = function(value) {private$Y},
    # covariates array
    covarArray = function(value) {private$phi},
    # covariates matrix
    covarMatrix = function(value) {if (missing(value)) return(private$X) else  private$X <- value},
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
    initialize = function(adjacencyMatrix, covarMatrix = NULL, covarArray = NULL) {

      ## adjacency matrix
      stopifnot(is.matrix(adjacencyMatrix))
      ## Only binary graph supported
      stopifnot(all.equal(unique(as.numeric(adjacencyMatrix[!is.na(adjacencyMatrix)])), c(0,1)))

      if (isSymmetric(adjacencyMatrix)) private$directed <- FALSE else private$directed <- TRUE
      private$Y  <- adjacencyMatrix

      ## covariates
      if (!is.null(covarMatrix)) {
        private$X <- covarMatrix
      }
      if (!is.null(covarArray)) {
        private$phi <- covarArray
      }

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
    }
  )
)

sampledNetwork$set("public", "plot",
function(main = paste("Network sampling with sampling rate:", signif(self$samplingRate,3))) {
  image_NA(self$adjacencyMatrix, main = main)
})

sampledNetwork$set("public", "show",
function(model = "Sampled Network\n") {
  cat(model)
  cat("==================================================================\n")
  cat("Structure for storing a sampled network in missSBM.\n")
  cat("==================================================================\n")
  cat("* Useful fields \n")
  cat("  $nNodes, $nDyads, $is_directed\n", "  $adjacencyMatrix, $covarMatrix, $covarArray\n",
      "  $dyads, $missingDyads, $observedDyads, $observedNodes\n",  "  $samplingRate, $samplingMatrix, $NAs\n")
  cat("* Useful method: plot() \n")
})
sampledNetwork$set("public", "print", function() self$show())
