#' Definition of R6 Class 'networkSampling_sampler'
#'
#' This class is use to define a sampling model for a network. Inherits from 'networkSampling'
#' it has a rSampling method which takes an adjacency matrix as an input and send back an object with class sampledNetwork.
#'
#' @include networkSampling-Class.R
#' @include utils.R
#' @import R6
#' @export
networkSampler <-
  R6::R6Class(classname = "networkSampler",
  inherit = networkSampling,
  private = list(
    N        = NULL, # number of nodes
    directed = NULL, # directed or undirected network
    R        = NULL  # the sampling matrix, indicating observed dyads
  ),
  public = list(
    initialize = function(type=NA, parameters=NA, nNodes=NA, directed=FALSE) {
      super$initialize(type, parameters)
      private$N <- nNodes
      private$directed <- directed
    },
    rSamplingMatrix = function() {
      D_obs <- private$rObservedDyads()
      R <- matrix(0, private$N, private$N)
      R[D_obs] <- 1; diag(R) <- 1
      if (!private$directed)  R <- t(R) | R
      private$R <- R
    }
  ),
  active = list(
    samplingMatrix = function(value) {private$R}
  )
)

# ================================================================================
# DYAD-CENTERED SAMPLINGS
#
#' Virtual class for all dyad-centered samplers
dyadSampler <-
R6::R6Class(classname = "dyadSampler",
  inherit = networkSampler,
  public = list(
    initialize = function(type = NA, parameters = NA, nNodes = NA, directed = FALSE) {
      super$initialize(type, parameters, nNodes, directed)
      tmp_mat <- matrix(NA, private$N, private$N)
      if (!directed) {
        private$dyads <- which(lower.tri(tmp_mat))
      } else {
        private$dyads <- which(lower.tri(tmp_mat) | upper.tri(tmp_mat))
      }
    }
  ),
  private = list(
    dyads = NULL,
    rObservedDyads = function() {
      success <- rbinom(length(private$dyads), 1, prob = private$rho)
      private$dyads[success == 1]
    }
  )
)

# ================================================================================
# NODE-CENTERED SAMPLINGS
#
# Virtual class for all node-centerd samplers
nodeSampler <-
R6::R6Class(classname = "nodeSampler",
  inherit = networkSampler,
  private = list(
    rObservedDyads = function() {
      N_obs <- which(runif(private$N) < private$rho)
      candidates <- expand.grid(N_obs, 1:private$N)
      D_obs <- unique(rbind(as.matrix(candidates), as.matrix(rev(candidates))))
      D_obs
    }
  )
)

## ================================================================================
## DYAD-CENTERED SAMPLINGS
##
## - DYAD (MAR)
## - DOUBLE-STANDARD (NMAR)
## - BLOCK-DYAD (NMAR)
simpleDyadSampler <-
R6::R6Class(classname = "simpleDyadSampler",
  inherit = dyadSampler,
  public = list(
    initialize = function(parameters = NA, nNodes = NA, directed = FALSE, covarArray = NULL) {
      super$initialize("dyad", parameters, nNodes, directed)
      if (is.null(covarArray)) {
        stopifnot(length(parameters) == 1, all(parameters >= 0), all(parameters <= 1))
        sampling_prob <- rep(parameters, length(private$dyads))
      } else {
        stopifnot(length(parameters) == dim(covarArray)[3])
        sampling_prob <- logistic(roundProduct(covarArray, parameters))[private$dyads]
      }
      private$rho <- sampling_prob
    }
  )
)

doubleStandardSampler <-
R6::R6Class(classname = "doubleStandardSampler",
  inherit = dyadSampler,
  public = list(
    initialize = function(parameters = NA, adjMatrix = NA, directed = FALSE) {
      stopifnot(length(parameters) == 2, all(parameters >= 0), all(parameters <= 1))
      super$initialize("double-standard", parameters, ncol(adjMatrix), directed)
      if (!directed) {
        D_0 <- which((adjMatrix == 0) & lower.tri(adjMatrix))
        D_1 <- which((adjMatrix == 1) & lower.tri(adjMatrix))
      } else {
        D_0 <- which((adjMatrix == 0) & (upper.tri(adjMatrix) | lower.tri(adjMatrix)))
        D_1 <- which((adjMatrix == 1) & (upper.tri(adjMatrix) | lower.tri(adjMatrix)))
      }
      sampling_prob <- rep(NA, length(private$dyads))
      sampling_prob[private$dyads %in% D_0] <- parameters[1]
      sampling_prob[private$dyads %in% D_1] <- parameters[2]
      private$rho <- sampling_prob
    }
  )
)

blockDyadSampler <-
R6::R6Class(classname = "blockDyadSampler",
  inherit = dyadSampler,
  private = list(Q = NULL),
  public = list(
    initialize = function(parameters = NA, nNodes = NA, directed = FALSE, clusters = NA) {
      Q <- length(unique(clusters))
      stopifnot(
        all(parameters >= 0),
        all(parameters <= 1),
        length(unique(clusters)) == Q,
        dim(parameters) == c(Q,Q)
      )
      super$initialize("block-dyad", parameters, nNodes, directed)
      Z <- matrix(0, nNodes, Q)
      Z[cbind(1:private$N, private$clusters)] <- 1
      private$rho <- (Z %*% self$parameters %*% t(Z))[private$dyads]
      private$Q <- Q
    }
  ),
  active = list(
    df = function(value) {ifelse(private$directed, private$Q^2, private$Q * (private$Q + 1) / 2) }
  )
)

## ================================================================================
## NODE-CENTERED SAMPLINGS
## - SIMPLE NODE (MAR)
## - BLOCK-NODE (NMAR)
## - DEGREE (NMAR)

simpleNodeSampler <-
R6::R6Class(classname = "simpleNodeSampler",
  inherit = nodeSampler,
  public = list(
    initialize = function(parameters = NA, nNodes = NA, directed = FALSE, covarMatrix = NULL) {
      ## w/o covariates
      if (is.null(covarMatrix)) {
        stopifnot(length(parameters) == 1, all(parameters >= 0), all(parameters <= 1))
        sampling_prob <- rep(parameters, nNodes)
      } else {
        stopifnot(length(parameters) == ncol(covarMatrix))
        sampling_prob <- logistic(covarMatrix %*% parameters)
      }
      super$initialize("node", parameters, nNodes, directed)
      private$rho <- sampling_prob
    }
  )
)

blockNodeSampler <-
R6::R6Class(classname = "blockNodeSampler",
  inherit = nodeSampler,
  public = list(
    initialize = function(parameters = NA, nNodes = NA, directed = FALSE, clusters = NA) {
      stopifnot(all(parameters >= 0), all(parameters <= 1))
      stopifnot(length(clusters) == nNodes, length(parameters) == length(unique(clusters)))
      super$initialize("block-node", parameters, nNodes, directed)
      private$rho <- self$parameters[clusters]
    }
  )
)

degreeSampler <-
R6::R6Class(classname = "degreeSampler",
  inherit = nodeSampler,
  public = list(
    initialize = function(parameters = NA, degrees = NA, directed = FALSE) {
      stopifnot(length(parameters) == 2)
      super$initialize("degree", parameters, length(degrees), directed)
      private$rho <- logistic(parameters[1] + self$parameters[2]*degrees)
    }
  )
)

### TODO: SNOWBALL SAMPLING add a parameter for the number of waves
# "snowball" = function(adjMatrix, ...) {
#   # initial set
#   N <- nrow(adjMatrix)
#   wave1 <- which(runif(N) < self$parameters)
#   # first wave
#   wave2 <- unique(unlist(adjMatrix[wave1, ] != 0, 1, function(x) which(x != 0)))
#   N_obs <- union(wave1, wave2)
#   D_obs <- unique(rbind(as.matrix(expand.grid(N_obs, 1:N)), as.matrix(rev(expand.grid(N_obs, 1:N)))))
# }
