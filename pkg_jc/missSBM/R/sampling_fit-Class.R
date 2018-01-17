#' Definition of R6 Class 'sampling_fit'
#'
#' This class is use to define a sampling fit. Inherits from 'sampling'
#'
#' @include sampling-Class.R
#' @import R6
#' @export
sampling_fit <-
R6Class(classname = "sampling_fit",
  inherit = sampling,
  public = list(
    initialize = function(adjMatrix) {
      private$net <- sampledNetwork$new(adjMatrix)
    }
  ),
  active = list(
    ## nDyads automatically handle the directed/undirected cases
    penalty = function(value) {log(private$net$nDyads) * self$df}
  )
)

#' @export
sampling_fit_dyad <-
R6Class(classname = "sampling_fit_dyad",
  inherit = sampling_fit,
  private = list(
    card_D_o = NULL, # stats required by the likelihood
    card_D_m = NULL  # number of observed, repectively missing dyads
  ),
  public = list(
    initialize = function(adjMatrix) {
      super$initialize(adjMatrix)
      private$name <- "dyad"
      private$card_D_o <- length(private$net$observedDyads)
      private$card_D_m <- length(private$net$missingDyads )
      private$psi <- private$card_D_o / private$net$nDyads
    }
  ),
  active = list(
    vLogLik = function() {
      res <- (private$card_D_o - private$net$nNodes) * logx(private$psi) + private$card_D_m * log1mx(private$psi)
      res
    }
  )
)

#' @export
sampling_fit_node <-
R6Class(classname = "sampling_fit_node",
  inherit = sampling_fit,
  private = list(
    card_N_o = NULL, # stats required by the likelihood
    card_N_m = NULL  # number of observed, repectively missing nodes
  ),
  public = list(
    initialize = function(adjMatrix) {
      super$initialize(adjMatrix)
      private$name <- "node"
      private$card_N_o <- sum( private$net$observedNodes)
      private$card_N_m <- sum(!private$net$observedNodes)
      private$psi <- private$card_N_o / (private$card_N_o + private$card_N_m)
    }
  ),
  active = list(
    vLogLik = function() {
### SHOULD'NT BE A SQUARE FOR N_obs SOMEWHERE?
### THE LOGLIK IS HIGH COMPARED TO THE ONE IN THE DYAD CASE...
      res <- private$card_N_o * logx(private$psi) + private$card_N_m * log1mx(private$psi)
      res
    }
  )
)

#' @export
sampling_fit_double_standard <-
R6Class(classname = "sampling_fit_double_standard",
  inherit = sampling_fit,
  private = list(
    So     = NULL, ## statistics only requiring the observed part of the network
    So.bar = NULL, ## can be computed once for all during the initialization
    Sm     = NULL, ## these ones will be updated during the algorithm
    Sm.bar = NULL
  ),
  public = list(
    initialize = function(adjMatrix, missingInit = NA) {
      super$initialize(adjMatrix)
      private$name <- "double_standard"
      private$So     <- sum(    private$net$adjacencyMatrix[private$net$observedDyads])
      private$So.bar <- sum(1 - private$net$adjacencyMatrix[private$net$observedDyads])
      ## SEE HOW TO "COMPLETE" THE NETWORK AT START-UP IN ORDER TO INITIALIZE PSI
      if (is.na(missingInit))
        missingInit <- rep(mean(private$net$adjacencyMatrix, na.rm = TRUE), length(private$net$missingDyads))
      self$update_missing(missingInit)
    },
    update_missing = function(nu) {
      private$Sm     <- sum(    nu)
      private$Sm.bar <- sum(1 - nu)
      private$psi    <- c(private$So.bar / (private$So.bar + private$Sm.bar), private$So / (private$So + private$Sm))
    }
  ),
  active = list(
    vLogLik = function(value) {
      res <- logx(private$psi[2]) * private$So + logx(private$psi[1]) * private$So.bar +
        log1mx(private$psi[2]) * private$Sm + log1mx(private$psi[1]) * private$Sm.bar
      res
    }
  )
)

#' @export
sampling_fit_block <-
R6Class(classname = "sampling_fit_block",
  inherit = sampling_fit,
  private = list(
    So = NULL, ## sum_(i in Nobs ) Z_iq
    Sm = NULL  ## sum_(i in Nmiss) Z_iq
  ),
  public = list(
    initialize = function(adjMatrix, blockInit) {
      super$initialize(adjMatrix)
      private$name <- "block"
      self$update_missing(blockInit)
    },
    update_missing = function(Z) {
      private$So <- colSums(Z[ private$net$observedNodes, ])
      private$Sm <- colSums(Z[!private$net$observedNodes, ])
      private$psi <- private$So / (private$So + private$Sm)
    }
  ),
  active = list(
    vLogLik = function() {
      res <- c(crossprod(private$So, log(private$psi)) +  crossprod(private$Sm, log(1 - private$psi)))
      res
    }
  )
)

#' @export
sampling_fit_degree <-
R6Class(classname = "sampling_fit_degree",
  inherit = sampling_fit,
  private = list(
    ksi      = NULL, # additional variational parameter for approximation in the logistic function
    degree_o = NULL  # estimation of the degrees on the observed part of the network
  ),
  public = list(
    initialize = function(adjMatrix) {
      super$initialize(adjMatrix)
      private$name <- "degree"

      ## will remain the same
      private$degree_o <- rowSums(private$net$adjacencyMatrix, na.rm = TRUE)

      ## will fluctuate along the algorithm
      private$psi <- coefficients(glm(1*(private$net$observedNodes) ~ private$degree_o, family = binomial(link = "logit")))
      nu <- matrix(NA, private$net$nNodes, private$net$nNodes)
      nu[private$net$missingDyads] <- mean(private$net$adjacencyMatrix, na.rm = TRUE)
    },
    update_missing = function(nu) {
      D  <- rowSums(nu, na.rm = TRUE) + private$degree_o
      D2 <- rowSums(nu * (1 - nu), na.rm = TRUE) + D^2
      private$ksi <- sqrt( private$psi[1]^2 + private$psi[2]^2 * D2  + 2 * private$psi[1] * private$psi[2] * D)
      C <- .5 * private$net$nNodes - sum(!private$net$observedNodes)
      s_hksi     <- sum(h(private$ksi))
      s_hksiD    <- sum(h(private$ksi) * D)
      s_hksiDhat <- sum(h(private$ksi) * D2)
      b <- (2 * C * s_hksiD  - (.5 * sum(D) - sum(D[!private$net$observedNodes])) * s_hksi) / (2 * s_hksiDhat * s_hksi - (2 * s_hksiD)^2)
      a <- -(b * s_hksiD + C) / s_hksi
      private$psi    <- c(a, b)
    }
  ),
  active = list(
    vLogLik = function() {
      prob <- logistic(private$psi[1] + private$psi[2] * private$D)
      res  <- log( prob^private$net$observedNodes %*% (1 - prob)^(!private$net$observedNodes) )
      res
    }
  )
)

sampling_fit_snowball <-
R6Class(classname = "sampling_fit_snowball",
  inherit = sampling_fit,
  active = list(
    vLogLik = function(value) {
      NA
    }
  )
)

