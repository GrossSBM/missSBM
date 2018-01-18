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
R6Class(classname = "networkSampling_sampler",
  inherit = networkSampling,
  public = list(
    ## methods
    initialize = function(type = NA, parameters = NA) {
      super$initialize(type = type)

      if (!switch(type,
                  "double_standard" = ifelse(length(parameters) == 2, TRUE, FALSE),
                  "degree"          = ifelse(length(parameters) == 2, TRUE, FALSE),
                  "dyad"            = ifelse(length(parameters) == 1, TRUE, FALSE),
                  "block"           = TRUE, ## handled in rSampling once the vector 'clusters' is known
                  "node"            = ifelse(length(parameters) == 1, TRUE, FALSE),
                  "snowball"        = ifelse(length(parameters) == 1, TRUE, FALSE))) {
        stop("Sampling parameters does not have the required size.")
      }

      if (type != "degree") {
        if (any(parameters < 0) | any(parameters > 1)) {
          stop("Sampling parameters must be probabilities (i.e between 0 and 1)")
        }
      }
      private$psi  <- parameters

      self$rSampling <- switch(type,
      "dyad" = function(adjMatrix, ...) {
        N <- ncol(adjMatrix)
        R <- diag(N)

        if (isSymmetric(adjMatrix)) {
          D <- which(lower.tri(adjMatrix))
        } else {
          D <- which(lower.tri(adjMatrix) | upper.tri(adjMatrix))
        }
        D_obs <- D[rbinom(length(D), 1, prob = self$parameters) == 1]
        R[D_obs] <- 1

        if (isSymmetric(adjMatrix))  R <- t(R) | R

        adjMatrix[R == 0] <- NA
        sampledNetwork$new(adjMatrix)
      },
      "double_standard" = function(adjMatrix, ...) {
        N <- ncol(adjMatrix)
        R <- diag(N)

        if (isSymmetric(adjMatrix)) {
          D_obs_0 <- which(runif(sum((adjMatrix == 0) & upper.tri(adjMatrix))) < self$parameters[1])
          D_obs_1 <- which(runif(sum((adjMatrix == 1) & upper.tri(adjMatrix))) < self$parameters[2])
        } else {
          D_obs_0 <- which(runif(sum(adjMatrix == 0)) < self$parameters[1])
          D_obs_1 <- which(runif(sum(adjMatrix == 1)) < self$parameters[2])
        }
## Check THIS...
        R[union(D_obs_0, D_obs_1)] <- 1
        if (isSymmetric(adjMatrix))  R <- t(R) | R

        adjMatrix[R == 0] <- NA
        sampledNetwork$new(adjMatrix)
      },
      "node" = function(adjMatrix, ...) {
        N <- ncol(adjMatrix)
        R <- diag(N)

        N_obs <- which(runif(N) < self$parameters)

        N_obs <- sample(1:N, floor(N * self$parameters))
        R[N_obs,] <- 1
        R <- t(R) | R

        adjMatrix[R == 0] <- NA
        sampledNetwork$new(adjMatrix)
      },
      "block" = function(adjMatrix, clusters) {
        N <- nrow(adjMatrix)
        if (length(clusters) != N)
          stop(paste("The vector 'clusters' must have ", N," entries!"))

        if (length(self$parameters) != length(unique(clusters)))
          stop("Sampling parameters does not have the required size.")

        R <- diag(N)

        N_obs <- which(runif(N) < self$parameters[clusters])
        R[N_obs,] <- 1
        R <- t(R) | R

        adjMatrix[R == 0] <- NA
        sampledNetwork$new(adjMatrix)
      },
      "degree" = function(adjMatrix, ...) {
        N <- nrow(adjMatrix)
        R <- diag(N)

        N_obs <- which(runif(N) < logistic(self$parameters[1] + self$parameters[2]*rowSums(adjMatrix)))
        R[N_obs,] <- 1
        R <- t(R) | R

        adjMatrix[R == 0] <- NA
        sampledNetwork$new(adjMatrix)
      },
      ### TODO: add a parameter for the number of waves
      "snowball" = function(adjMatrix, ...) {
        N <- nrow(adjMatrix)
        R <- diag(N)

        # First wave
        wave1 <- which(runif(N) < self$parameters)
        # Second wave
        wave2 <- unique(unlist(adjMatrix[wave1, ] != 0, 1, function(x) which(x != 0)))
        N_obs <- union(wave1, wave2)

        R[N_obs,] <- 1; R[,N_obs] <- 1

        adjMatrix[R == 0] <- NA
        sampledNetwork$new(adjMatrix)
      })
    },
    rSampling = NULL ## the sampling function (initialize with in the constructor)
  )
)
