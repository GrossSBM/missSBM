#' Definition of R6 Class 'networkSampling_sampler'
#'
#' This class is use to define a sampling model for a network. Inherits from 'networkSampling'
#' it has a rSampling method which takes an adjacency matrix as an input and send back an object with class sampledNetwork.
#'
#' @include networkSampling-Class.R
#' @include sampledNetwork-Class.R
#' @include utils.R
#' @import R6
#' @export
networkSampling_sampler <-
  R6::R6Class(classname = "networkSampling_sampler",
  inherit = networkSampling,
  private = list(
    ## a private function to generate a collection of observed dyads for different sampling process
    rObservedDyads = NULL
  ),
  public = list(
    ## methods
    initialize = function(type = NA, parameters = NA) {
      super$initialize(type = type)

      if (!switch(type,
                  "double_standard" = ifelse(length(parameters) == 2, TRUE, FALSE),
                  "block_dyad"      = TRUE, # handled in rSampling once the vector 'clusters' is known
                  "degree"          = ifelse(length(parameters) == 2, TRUE, FALSE),
                  "dyad"            = TRUE, # handled in rSampling once the the number of dyads is known
                  "block"           = TRUE, # handled in rSampling once the vector 'clusters' is known
                  "node"            = TRUE, # handled in rSampling once the the number of nodes is known
                  ## "snowball"        = ifelse(length(parameters) == 1, TRUE, FALSE))
                  )) {
        stop("Sampling parameters does not have the required size.")
      }

      if (type != "degree") {
        if (any(parameters < 0) | any(parameters > 1)) {
          stop("Sampling parameters must be probabilities (i.e between 0 and 1)")
        }
      }
      private$psi  <- parameters

      ## a private function to generate a collection of
      ## observed dyads for different sampling process
      private$rObservedDyads <- switch(type,
      "dyad" = function(adjMatrix, ...) {
        if (isSymmetric(adjMatrix)) {
          D <- which(lower.tri(adjMatrix))
        } else {
          D <- which(lower.tri(adjMatrix) | upper.tri(adjMatrix))
        }
        if (length(self$parameters) == 1) {
          D_obs <- D[rbinom(length(D), 1, prob = self$parameters) == 1]
        } else {
          if (length(self$parameters) != (ncol(adjMatrix)^2) ) {
            stop("Sampling parameters does not have the required size.")
          }
          D_obs <- D[rbinom(length(D), 1, prob = self$parameters[D]) == 1]
        }
      },
      "double_standard" = function(adjMatrix, ...) {
        if (isSymmetric(adjMatrix)) {
          D_0 <- which((adjMatrix == 0) & upper.tri(adjMatrix))
          D_1 <- which((adjMatrix == 1) & upper.tri(adjMatrix))
        } else {
          D_0 <- which((adjMatrix == 0) & (upper.tri(adjMatrix) | lower.tri(adjMatrix)))
          D_1 <- which((adjMatrix == 1) & (upper.tri(adjMatrix) | lower.tri(adjMatrix)))
        }
        D_obs_0 <- D_0[runif(length(D_0)) < self$parameters[1]]
        D_obs_1 <- D_1[runif(length(D_1)) < self$parameters[2]]
        D_obs <- union(D_obs_0, D_obs_1)
      },
      "block_dyad" = function(adjMatrix, clusters) {
        N <- ncol(adjMatrix)
        Q <- length(unique(clusters))
        if (length(self$parameters) != Q*Q)
          stop("Sampling parameters does not have the required size.")
        Z <- matrix(0, N, Q); Z[cbind(1:N,clusters)] <- 1
        D_obs <- matrix(rbinom(N^2, 1, Z %*% self$parameters %*% t(Z)), N)
        if (isSymmetric(adjMatrix)){
          D_obs <- D_obs * lower.tri(D_obs) + t(D_obs * lower.tri(D_obs))
          D_obs <- which(D_obs == 1)
        } else {
          D_obs <- D_obs * lower.tri(D_obs) + t(D_obs * upper.tri(D_obs))
          D_obs <- which(D_obs == 1)
        }
      },
      "node" = function(adjMatrix, ...) {
        N <- ncol(adjMatrix)
        if (length(self$parameters) != 1 & length(self$parameters) != N)
            stop("Sampling parameters does not have the required size.")
  ##      sampling_prob <- logistic(private$X %*% self$parameters)
        N_obs <- which(runif(N) < self$parameters)
        D_obs <- unique(rbind(as.matrix(expand.grid(N_obs, 1:N)), as.matrix(rev(expand.grid(N_obs, 1:N)))))

      },
      "block" = function(adjMatrix, clusters) {
        N <- nrow(adjMatrix)
        if (length(clusters) != N)
          stop(paste("The vector 'clusters' must have ", N," entries!"))

        if (length(self$parameters) != length(unique(clusters)))
          stop("Sampling parameters does not have the required size.")
        N_obs <- which(runif(N) < self$parameters[clusters])
        D_obs <- unique(rbind(as.matrix(expand.grid(N_obs, 1:N)), as.matrix(rev(expand.grid(N_obs, 1:N)))))
      },
      "degree" = function(adjMatrix, ...) {
        N <- nrow(adjMatrix)
        N_obs <- which(runif(N) < logistic(self$parameters[1] + self$parameters[2]*rowSums(adjMatrix)))
        D_obs <- unique(rbind(as.matrix(expand.grid(N_obs, 1:N)), as.matrix(rev(expand.grid(N_obs, 1:N)))))
      } #,
      # ### TODO: add a parameter for the number of waves
      # "snowball" = function(adjMatrix, ...) {
      #   # initial set
      #   N <- nrow(adjMatrix)
      #   wave1 <- which(runif(N) < self$parameters)
      #   # first wave
      #   wave2 <- unique(unlist(adjMatrix[wave1, ] != 0, 1, function(x) which(x != 0)))
      #   N_obs <- union(wave1, wave2)
      #   D_obs <- unique(rbind(as.matrix(expand.grid(N_obs, 1:N)), as.matrix(rev(expand.grid(N_obs, 1:N)))))
      # }
      )
    },
    rSampling = function(adjMatrix, ...) {
      D_obs <- private$rObservedDyads(adjMatrix, ...)
      R <- matrix(0,ncol(adjMatrix), ncol(adjMatrix))
      R[D_obs] <- 1; diag(R) <- 1
      if (isSymmetric(adjMatrix))  R <- t(R) | R
      adjMatrix[R == 0] <- NA
      sampledNetwork$new(adjMatrix)
    }
  )
)

