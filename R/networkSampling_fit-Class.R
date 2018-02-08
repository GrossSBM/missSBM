#' Virtual class use to define a family of networkSampling_fit
#' @include networkSampling-Class.R
#' @import R6
networkSampling_fit <-
R6Class(classname = "networkSampling_fit",
  inherit = networkSampling,
  private = list(
    NAs      = NULL, # where are the missing entries
    dyads    = NULL, # set of dyads (differ according to directed/undirected graphs)
    card_D_o = NULL, # stats required by the likelihood
    card_D   = NULL  # number of possible dyads in the network
  ),
  public = list(
    ## initialize estimation and imputation function
    ## by default, nothing to do (corresponds to MAR sampling)
    update_parameters = function(...) {invisible(NULL)},
    update_imputation = function(...) {NA}
  ),
  active = list(
    ## nDyads automatically handles the directed/undirected cases
    penalty = function(value) {log(private$card_D) * self$df}
  )
)

#' Definition of R6 Class 'networkSampling_fit'
#'
#' This class is use to define a fit for a networkSampling. Inherits from 'networkSampling'
#' @export
dyadSampling_fit <-
R6Class(classname = "dyadSampling_fit",
  inherit = networkSampling_fit,
  private = list(
    card_D_m = NULL  # number of nodes, number of observed, repectively missing dyads
  ),
  public = list(
    initialize = function(sampledNetwork, ...) {
      private$name <- "dyad"
      private$card_D_o <- length(sampledNetwork$observedDyads)
      private$card_D_m <- length(sampledNetwork$missingDyads )
      private$card_D   <- sampledNetwork$nDyads
      private$psi      <- private$card_D_o / private$card_D
    }
  ),
  active = list(
    logLik = function() {
      res <- private$card_D_o * logx(private$psi) + private$card_D_m * log1mx(private$psi)
      res
    },
    ## overloading the penalty in MAR case
    penalty = function(value) {log(private$card_D_o) * self$df}
  )
)

#' @export
nodeSampling_fit <-
R6Class(classname = "nodeSampling_fit",
  inherit = networkSampling_fit,
  private = list(
    card_N_o = NULL, # stats required by the likelihood
    card_N_m = NULL, # number of observed, repectively missing nodes
    N        = NULL  #
  ),
  public = list(
    initialize = function(sampledNetwork, ...) {
      private$name <- "node"
      private$card_N_o <- sum( sampledNetwork$observedNodes)
      private$card_N_m <- sum(!sampledNetwork$observedNodes)
      private$card_D_o <- length(sampledNetwork$observedDyads)
      private$card_D   <- sampledNetwork$nDyads
      private$psi <- private$card_N_o / (private$card_N_o + private$card_N_m)
    }
  ),
  active = list(
    logLik = function() {
      res <- private$card_N_o * logx(private$psi) + private$card_N_m * log1mx(private$psi)
      res
    },
    ## overloading the penalty in MAR case
    penalty = function(value) {log(private$card_N_o) * self$df}
  )
)

#' @export
doubleStandardSampling_fit <-
R6Class(classname = "doubleStandardSampling_fit",
  inherit = networkSampling_fit,
  private = list(
    So     = NULL, ## statistics only requiring the observed part of the network
    So.bar = NULL, ## can be computed once for all during the initialization
    Sm     = NULL, ## these ones will be updated during the optimization
    Sm.bar = NULL
  ),
  public = list(
    initialize = function(sampledNetwork, ...) {
      private$name <- "double_standard"
      private$card_D <- sampledNetwork$nDyads
      private$NAs    <- sampledNetwork$NAs
      private$dyads  <- sampledNetwork$dyads
      private$So     <- sum(    sampledNetwork$adjacencyMatrix[sampledNetwork$observedDyads])
      private$So.bar <- sum(1 - sampledNetwork$adjacencyMatrix[sampledNetwork$observedDyads])
      imputedNet     <- matrix(mean(sampledNetwork$adjacencyMatrix, na.rm = TRUE), sampledNetwork$nNodes, sampledNetwork$nNodes)
      self$update_parameters(imputedNet)
    },
    update_parameters = function(ImputedNet, ...) {
      private$Sm     <- sum(    ImputedNet[private$NAs & private$dyads])
      private$Sm.bar <- sum(1 - ImputedNet[private$NAs & private$dyads])
      private$psi    <- c(private$So.bar / (private$So.bar + private$Sm.bar), private$So / (private$So + private$Sm))
    },
    update_imputation = function(Z, pi) {
      nu <- logistic(log((1 - private$psi[2]) / (1 - private$psi[1])) + Z %*% log(pi/(1 - pi)) %*% t(Z))
      nu[private$NAs]
    }
  ),
  active = list(
    logLik = function(value) {
      res <- logx(private$psi[2]) * private$So + logx(private$psi[1]) * private$So.bar +
        log1mx(private$psi[2]) * private$Sm + log1mx(private$psi[1]) * private$Sm.bar
      res
    }
  )
)

#' @export
blockSampling_fit <-
R6Class(classname = "blockSampling_fit",
  inherit = networkSampling_fit,
  private = list(
    N_obs  = NULL, ## boolean for observed nodes
    So     = NULL, ## sum_(i in Nobs ) Z_iq
    Sm     = NULL  ## sum_(i in Nmiss) Z_iq
  ),
  public = list(
    initialize = function(sampledNetwork, blockInit) {
      private$name <- "block"
      private$card_D <- sampledNetwork$nDyads
      private$NAs    <- sampledNetwork$NAs
      private$N_obs  <- sampledNetwork$observedNodes
      self$update_parameters(NA, blockInit)
    },
    update_parameters = function(imputedNet, Z) {
      private$So <- colSums(Z[ private$N_obs, ])
      private$Sm <- colSums(Z[!private$N_obs, ])
      private$psi <- private$So / (private$So + private$Sm)
    },
    update_imputation = function(Z, pi) {
      nu <- logistic(Z %*% log(pi/(1 - pi)) %*% t(Z))
      nu[private$NAs]
    }
  ),
  active = list(
    logLik = function() {
      res <- c(crossprod(private$So, log(private$psi)) +  crossprod(private$Sm, log(1 - private$psi)))
      res
    }
  )
)

#' @export
degreeSampling_fit <-
R6Class(classname = "degreeSampling_fit",
  inherit = networkSampling_fit,
  private = list(
    Dij      = NULL, # stat required to perform imputation
    N_obs    = NULL, # observed nodes
    ksi      = NULL, # additional variational parameter for approximation in the logistic function
    degree_o = NULL  # estimation of the degrees on the observed part of the network
  ),
  public = list(
    initialize = function(sampledNetwork, ...) {

      private$name <- "degree"
      private$card_D <- sampledNetwork$nDyads
      private$NAs    <- sampledNetwork$NAs
      private$N_obs  <- sampledNetwork$observedNodes

      ## will remain the same
      private$degree_o <- rowSums(sampledNet$adjacencyMatrix, na.rm = TRUE)

      ## will fluctuate along the algorithm
      private$psi <- coefficients(glm(1*(sampledNet$observedNodes) ~ private$degree_o, family = binomial(link = "logit")))

      imputedNet     <- matrix(mean(sampledNetwork$adjacencyMatrix, na.rm = TRUE), sampledNetwork$nNodes, sampledNetwork$nNodes)
      self$update_parameters(imputedNet)
    },
    update_parameters = function(imputedNet, ...) {
      nu <- imputedNet
      nu[!private$NAs] <- NA
      D  <- rowSums(nu, na.rm = TRUE) + private$degree_o
      D2 <- rowSums(nu * (1 - nu), na.rm = TRUE) + D^2
      private$ksi <- sqrt( private$psi[1]^2 + private$psi[2]^2 * D2  + 2 * private$psi[1] * private$psi[2] * D)
      C <- .5 * nrow(imputedNet) - sum(!private$N_obs)
      s_hksi     <- sum(h(private$ksi))
      s_hksiD    <- sum(h(private$ksi) * D)
      s_hksiDhat <- sum(h(private$ksi) * D2)
      b <- (2 * C * s_hksiD  - (.5 * sum(D) - sum(D[!private$N_obs])) * s_hksi) / (2 * s_hksiDhat * s_hksi - (2 * s_hksiD)^2)
      a <- -(b * s_hksiD + C) / s_hksi
      private$psi    <- c(a, b)

      ## update stat required to perform imputation
      nu[!private$NAs] <- 0
      private$Dij <- matrix(D, nrow(imputedNet), ncol(imputedNet)) - nu
    },
    update_imputation = function(Z, pi) {
      C <- 2 * h(private$ksi) * (private$ksi[1] * private$ksi[2] + private$ksi[2]^2 * (1 + private$Dij))
      nu <- logistic(Z %*% log(pi/(1 - pi)) %*% t(Z) - private$psi[2] + C + t(C) )
      nu[private$NAs]
    }
  ),
  active = list(
    logLik = function() {
      prob <- logistic(private$psi[1] + private$psi[2] * private$D)
      ## ????
      res  <- log( prob^private$N_obs %*% (1 - prob)^(!private$N_obs) )
      res
    }
  )
)

## TODO: handle multiple waves
snowballSampling_fit <-
R6Class(classname = "snowballSampling_fit",
  inherit = networkSampling_fit,
  active = list(
    logLik = function(value) {
      NA
    }
  )
)

