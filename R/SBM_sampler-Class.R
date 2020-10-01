#' An R6 Class to represent a sampler for a SBM
#'
#' The function [`simulate`] produces an instance of an object with class [`SBM_sampler`].
#'
#' All fields of this class are only accessible for reading. This class comes with a set of methods,
#' some of them being useful for the user (see examples)
#' \itemize{
#' \item{R6 methods:}{$rMemberships(), $rAdjacency()}
#' \item{S3 methods}{print(), plot()}
#' }
#'
#' Fields are accessed via active binding and cannot be changed by the user.
#'
#' @inherit SBM details
#'
#' @include SBM-Class.R
#' @importFrom R6 R6Class
#' @examples
#' ## SBM parameters
#' directed <- FALSE
#' N <- 300 # number of nodes
#' Q <- 3   # number of clusters
#' pi <- rep(1,Q)/Q     # mixture parameter
#' theta <- diag(.45,Q) + .05 # connectivity matrix
#'
#' ## draw a SBM without covariates through simulateSBM
#' sbm <- missSBM::simulate(N, pi, theta, directed)
#'
#' ## equivalent construction from the SBM_sampler class itslef
#' sbm_s <- SBM_sampler$new(directed, N, pi, theta)
#' sbm_s$rMemberships() # draw some blocks
#' sbm_s$rAdjacency() # draw some edges
#'
#' coef(sbm_s, "block")
#' coef(sbm_s, "connectivity")
#' summary(sbm_s)
#' @seealso The function \code{\link{simulate}}.
#' @export
SBM_sampler <-
  R6::R6Class(classname = "SBM_sampler",
  inherit = SBM,
  ## fields for internal use (refering to mathematical notations)
  private = list(
    Z     = NULL # the sampled indicator of blocks
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Initialize a [`SBM_sampler`] object
    #' @param directed logical indicating if the network is directed or not
    #' @param nbNodes The number of nodes
    #' @param blockProp the vector of mixture parameters
    #' @param connectParam the matrix of connectivity: inter/intra probabilities of connection when the network does not have covariates, or a logit scaled version of it.
    #' @param covarParam the vector of parameters associated with the covariates
    #' @param covarList A list with M entries (the M covariates). Each entry of the list must be an N x N matrix
    initialize = function(directed = FALSE, nbNodes=NA, blockProp=NA, connectParam=NA, covarParam=numeric(0), covarList=list()) {
      super$initialize(directed, nbNodes, blockProp, connectParam, covarParam, covarList)
    },
    #' @description a method to generate a vector of clusters indicators
    rMemberships = function() {
      private$Z <- t(rmultinom(private$N, size = 1, prob = private$pi))
    },
    #' @description a method to sample an adjacency matrix for the current SBM
    rAdjacency = function() {
      Y <- matrix(rbinom(private$N^2, 1, self$expectation), private$N)
      if (!private$directed_) Y <- Y * lower.tri(Y) + t(Y * lower.tri(Y))
      private$Y <- Y
    },
    #' @description show method
    show = function() {
      super$show("Sampler for Stochastic Block Model\n")
      cat("* Sampling methods \n")
      cat("  $rMemberships(), $rAdjacency()\n")
      cat("* Additional fields \n")
      cat("  $indMemberships, $memberships, $adjMatrix, $expectation\n")
      cat("* S3 methods \n")
      cat("  plot, print, summary, coef\n")
    },
    #' @description User friendly print method
    print = function() { self$show() }
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## ACTIVE BINDING
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field indMemberships matrix for clustering memberships
    indMemberships = function(value) {private$Z},
    #' @field memberships vector of clustering
    memberships = function(value) {apply(private$Z, 1, which.max)},
    #' @field expectation expected values of connection under the current model
    expectation = function(value) {
      PI <- private$Z %*% private$theta %*% t(private$Z)
      if (self$nbCovariates > 0) {
        PI <- .logistic(PI + self$covarEffect)
      }
      PI
    }
  )
)
