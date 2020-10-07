#' Definition of R6 Class 'networkSampling_sampler'
#'
#' This class is use to define a sampling model for a network. Inherits from 'networkSampling'.
#' Owns a rSampling method which takes an adjacency matrix as an input and send back an object
#' with class partlyObservedNetwork.
#'
#' @include R6Class-networkSampling.R
#' @include utils_missSBM.R
#' @seealso \code{\link{partlyObservedNetwork}}
#' @import R6
networkSampler <-
  R6::R6Class(classname = "networkSampler",
  inherit = networkSampling,
  private = list(
    N        = NULL, # number of nodes
    directed = NULL, # directed or undirected network
    R        = NULL  # the sampling matrix, indicating observed dyads
  ),
  public = list(
    #' @description constructor for networkSampling
    #' @param type character for the type of sampling. must be in ("dyad", "covar-dyad", "node", "covar-node", "block-node", "block-dyad", "double-standard", "degree")
    #' @param parameters the vector of parameters associated to the sampling at play
    #' @param nbNodes number of nodes in the network
    #' @param directed logical, directed network of not
    initialize = function(type=NA, parameters=NA, nbNodes=NA, directed=FALSE) {
      super$initialize(type, parameters)
      private$N <- nbNodes
      private$directed <- directed
    },
    #' @description a method for drawing a sampling matrix according to the current sampling design
    rSamplingMatrix = function() {
      D_obs <- private$rObservedDyads()
      R <- matrix(0, private$N, private$N)
      R[D_obs] <- 1; diag(R) <- 1
      if (!private$directed)  R <- t(R) | R
      private$R <- R
    }
  ),
  active = list(
    #' @field samplingMatrix a matrix of booleans indicating sampled entries
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
    #' @description constructor for networkSampling
    #' @param type character for the type of sampling. must be in ("dyad", "covar-dyad", "node", "covar-node", "block-node", "block-dyad", "double-standard", "degree")
    #' @param parameters the vector of parameters associated to the sampling at play
    #' @param nbNodes number of nodes in the network
    #' @param directed logical, directed network of not
    initialize = function(type = NA, parameters = NA, nbNodes = NA, directed = FALSE) {
      super$initialize(type, parameters, nbNodes, directed)
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

#' Virtual class for all node-centered samplers
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

#' Class for defining a simple dyad sampler
simpleDyadSampler <-
R6::R6Class(classname = "simpleDyadSampler",
  inherit = dyadSampler,
  public = list(
    #' @description constructor for networkSampling
    #' @param parameters the vector of parameters associated to the sampling at play
    #' @param nbNodes number of nodes in the network
    #' @param directed logical, directed network of not
    #' @param covarArray an array of covariates used
    #' @param intercept double, intercept term used to compute the probability of sampling in the presence of covariates. Default 0.
    initialize = function(parameters = NA, nbNodes = NA, directed = FALSE, covarArray = NULL, intercept = 0) {
      super$initialize("dyad", parameters, nbNodes, directed)
      if (is.null(covarArray)) {
        stopifnot(length(parameters) == 1, all(parameters >= 0), all(parameters <= 1))
        sampling_prob <- rep(parameters, length(private$dyads))
      } else {
        stopifnot(length(parameters) == dim(covarArray)[3])
        sampling_prob <- .logistic(intercept + roundProduct(array2list(covarArray), parameters))[private$dyads]
      }
      private$rho <- sampling_prob
    }
  )
)

#' Class for defining a double-standard sampler
doubleStandardSampler <-
R6::R6Class(classname = "doubleStandardSampler",
  inherit = dyadSampler,
  public = list(
    #' @description constructor for networkSampling
    #' @param parameters the vector of parameters associated to the sampling at play
    #' @param adjMatrix matrix of adjacency
    #' @param directed logical, directed network of not
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

#' Class for defining a block dyad sampler
blockDyadSampler <-
R6::R6Class(classname = "blockDyadSampler",
  inherit = dyadSampler,
  private = list(Q = NULL),
  public = list(
    #' @description constructor for networkSampling
    #' @param parameters the vector of parameters associated to the sampling at play
    #' @param nbNodes number of nodes in the network
    #' @param directed logical, directed network of not
    #' @param clusters a vector of class memberships
    initialize = function(parameters = NA, nbNodes = NA, directed = FALSE, clusters = NA) {
      Q <- length(unique(clusters))
      stopifnot(
        all(parameters >= 0),
        all(parameters <= 1),
        dim(parameters) == c(Q,Q)
      )
      super$initialize("block-dyad", parameters, nbNodes, directed)
      Z <- matrix(0, nbNodes, Q)
      Z[cbind(1:private$N, private$clusters)] <- 1
      private$rho <- (Z %*% self$parameters %*% t(Z))[private$dyads]
      private$Q <- Q
    }
  ),
  active = list(
    #' @field df the number of parameters of this sampling
    df = function(value) {ifelse(private$directed, private$Q^2, private$Q * (private$Q + 1) / 2) }
  )
)

## ================================================================================
## NODE-CENTERED SAMPLINGS
## - SIMPLE NODE (MAR)
## - BLOCK-NODE (NMAR)
## - DEGREE (NMAR)

#' Class for defining a simple node sampler
simpleNodeSampler <-
R6::R6Class(classname = "simpleNodeSampler",
  inherit = nodeSampler,
  public = list(
    #' @description constructor for networkSampling
    #' @param parameters the vector of parameters associated to the sampling at play
    #' @param nbNodes number of nodes in the network
    #' @param directed logical, directed network of not
    #' @param covarMatrix a matrix of covariates used
    #' @param intercept double, intercept term used to compute the probability of sampling in the presence of covariates. Default 0.
    initialize = function(parameters = NA, nbNodes = NA, directed = FALSE, covarMatrix = NULL, intercept = 0) {
      ## w/o covariates
      if (is.null(covarMatrix)) {
        stopifnot(length(parameters) == 1, all(parameters >= 0), all(parameters <= 1))
        sampling_prob <- rep(parameters, nbNodes)
      } else {
        stopifnot(length(parameters) == ncol(covarMatrix))
        sampling_prob <- .logistic(intercept + covarMatrix %*% parameters)
      }
      super$initialize("node", parameters, nbNodes, directed)
      private$rho <- sampling_prob
    }
  )
)

#' Class for defining a block node sampler
blockNodeSampler <-
R6::R6Class(classname = "blockNodeSampler",
  inherit = nodeSampler,
  public = list(
    #' @description constructor for networkSampling
    #' @param parameters the vector of parameters associated to the sampling at play
    #' @param nbNodes number of nodes in the network
    #' @param directed logical, directed network of not
    #' @param clusters a vector of class memberships
    initialize = function(parameters = NA, nbNodes = NA, directed = FALSE, clusters = NA) {
      stopifnot(all(parameters >= 0), all(parameters <= 1))
      stopifnot(length(clusters) == nbNodes, length(parameters) == length(unique(clusters)))
      super$initialize("block-node", parameters, nbNodes, directed)
      private$rho <- self$parameters[clusters]
    }
  )
)

#' Class for defining a degree sampler
degreeSampler <-
R6::R6Class(classname = "degreeSampler",
  inherit = nodeSampler,
  public = list(
    #' @description constructor for networkSampling
    #' @param parameters the vector of parameters associated to the sampling at play
    #' @param degrees vector of nodes' degrees
    #' @param directed logical, directed network of not
    initialize = function(parameters = NA, degrees = NA, directed = FALSE) {
      stopifnot(length(parameters) == 2)
      super$initialize("degree", parameters, length(degrees), directed)
      private$rho <- .logistic(parameters[1] + self$parameters[2]*degrees)
    }
  )
)

#' Class for defining a snowball sampler
snowballSampler <-
  R6::R6Class(classname = "snowballSampler",
    inherit = nodeSampler,
    public = list(
      #' @description constructor for networkSampling
      #' @param parameters the vector of parameters associated to the sampling at play
      #' @param adjacencyMatrix the adjacency matrix of the network
      #' @param directed logical, directed network of not
      initialize = function(parameters = NA, adjacencyMatrix=NA, directed = FALSE) {
        stopifnot(length(parameters) == 2)
        n <- nrow(adjacencyMatrix)
        nWaves <- parameters[1] # number of waves
        pfirstwave <- parameters[2] # proportion of nodes seen in the first wave
        # wave 1
        observedNodes <- (runif(n) < pfirstwave)*1
        nRemainingWaves <- nWaves - 1
        while (nRemainingWaves>0 & sum(observedNodes)<n)
        {
          link <- adjacencyMatrix[which(observedNodes==1),,drop=FALSE]
          observedNodes[which(colSums(link)>0)] <- 1 # link tracing in giver to receiver
          nRemainingWaves <- nRemainingWaves - 1
        }
        super$initialize("snowball", parameters, n, directed)
        private$rho <- observedNodes # 0-1 probabilities
      }
    )
  )
