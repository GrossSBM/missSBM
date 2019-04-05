#' An R6 Class to represent a sampler for a SBM
#'
#' The function \code{\link{simulate}} produces an instance of an object with class \code{SBM_sampler}.
#'
#' All fields of this class are only accessible for reading. This class comes with a set of methods,
#' some of them being useful for the user (see examples)
#' \itemize{
#' \item{R6 methods:}{$rBlocks(), $rAdjancencyMatrix()}
#' \item{S3 methods}{print(), plot()}
#' }
#'
#' @field nNodes The number of nodes
#' @field nBlocks The number of blocks
#' @field nCovariates  The number of covariates
#' @field nDyads The number of possible dyad in the network (depends on the direction)
#' @field direction A character indicating if the network is directed or undirected
#' @field hasCovariates a boolean indicating if the model has covariates
#' @field mixtureParam the vector of mixture parameters
#' @field connectParam the matrix of connectivity: inter/intra probabilities of connection when the network does not have covariates, or a logit scaled version of it.
#' @field covarParam the vector of parameters associated with the covariates
#' @field covarArray the array of covariates
#'
#' @include SBM-Class.R
#' @importFrom R6 R6Class
#' @examples
#' ## SBM parameters
#' directed <- FALSE
#' N <- 300 # number of nodes
#' Q <- 3   # number of clusters
#' alpha <- rep(1,Q)/Q     # mixture parameter
#' pi <- diag(.45,Q) + .05 # connectivity matrix
#' gamma <- log(pi/(1-pi)) # logit transform fo the model with covariates
#' M <- 2 # two Gaussian covariates
#' covarMatrix <- matrix(rnorm(N*M,mean = 0, sd = 1), N, M)
#' covarParam  <- rnorm(M, -1, 1)
#'
#' ## draw a SBM without covariates through simulateSBM
#' sbm <- missSBM::simulate(N, alpha, pi, directed)
#'
#' ## equivalent construction from the SBM_sampler class itslef
#' sbm_s <- SBM_sampler$new(directed, N, alpha, pi)
#' sbm_s$rBlocks() # draw some blocks
#' sbm_s$rAdjMatrix() # draw some edges
#'
#' @seealso The function \code{\link{simulate}}.
#' @export
SBM_sampler <-
  R6::R6Class(classname = "SBM_sampler",
  inherit = SBM,
  ## fields for internal use (refering to mathematical notations)
  private = list(
    Z     = NULL # the sampled indicator of blocks
  ),
  public = list(
    initialize = function(directed = FALSE, nNodes=NA, mixtureParam=NA, connectParam=NA, covarParam=NULL, covarArray=NULL) {
      super$initialize(directed, nNodes, mixtureParam, connectParam, covarParam, covarArray)
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
    connectProb = function(value) {
      PI <- private$Z %*% private$pi %*% t(private$Z)
      if (self$hasCovariates) {
        PI <- logistic(PI + roundProduct(private$X, private$beta))
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
  cat("  $blocks, $memberships, $adjMatrix, $connectProb\n")
  cat("* S3 methods \n")
  cat("  plot, print\n")
})

SBM_sampler$set("public", "print", function() self$show())

SBM_sampler$set("public", "plot",
  function(type = c("network", "connectivity")) {
    type <- match.arg(type)
    if (type == "network") {
      image_NA(self$adjacencyMatrix[order(self$memberships), order(self$memberships)])
    }
    if (type == "connectivity") {
      plot(
        igraph::graph_from_adjacency_matrix(
          private$pi,
          mode = ifelse(private$directed, "directed", "undirected"),
          weighted = TRUE, diag = TRUE
        ), main = "SBM connectivity matrix"
      )
    }
  }
)
