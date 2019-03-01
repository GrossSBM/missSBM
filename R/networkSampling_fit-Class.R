#' Virtual class used to define a family of networkSamplingDyads_fit
#' @include networkSampling-Class.R
#' @import R6
networkSamplingDyads_fit <-
  R6::R6Class(classname = "networkSamplingDyads_fit",
  inherit = networkSampling,
  private = list(
    card_D = NULL, # number of possible dyads in the network
    D_miss = NULL  # where are the missing dyads
  ),
  public = list(
    initialize = function(sampledNetwork, name) {
      private$name    <- name
      private$D_miss  <- sampledNetwork$missingDyads
      private$card_D  <- sampledNetwork$nDyads
    },
    ## initialize estimation and imputation function
    ## by default, nothing to do (corresponds to MAR sampling)
    update_parameters = function(...) {invisible(NULL)},
    update_imputation = function(Z, pi) { ## good for MCAR on node, dyads and NMAR with blocks
      nu <- check_boundaries(logistic(Z %*% (log(pi/(1 - pi))) %*% t(Z)))
      nu
    }
  ),
  active = list(
    ## nDyads automatically handles the directed/undirected cases
    penalty = function(value) {log(private$card_D) * self$df},
##    entropy = function(value) {-sum(xlogx(private$nu[private$NAs]) + xlogx(1 - private$nu[private$NAs]))}
    log_lambda = function(value) {0}
  )
)

#' Virtual class used to define a family of networkSamplingNodes_fit
#' @include networkSampling-Class.R
#' @import R6
networkSamplingNodes_fit <-
  R6::R6Class(classname = "networkSamplingNodes_fit",
  inherit = networkSampling,
  private = list(
    card_N = NULL, # number of nodes in the network
    N_obs  = NULL  # boolean for observed nodes
  ),
  public = list(
    initialize = function(sampledNetwork, name) {
      private$name   <- name
      private$N_obs  <- sampledNetwork$observedNodes
      private$card_N <- sampledNetwork$nNodes
    },
    ## initialize estimation and imputation function
    ## by default, nothing to do (corresponds to MAR sampling)
    update_parameters = function(...) {invisible(NULL)},
    update_imputation = function(Z, pi) { ## good for MCAR on node, dyads and NMAR with blocks (dyad and node)
      nu <- check_boundaries(logistic(Z %*% (log(pi/(1 - pi))) %*% t(Z)))
      nu
    }
  ),
  active = list(
    ## nDyads automatically handles the directed/undirected cases
    penalty = function(value) {log(private$card_N) * self$df},
    log_lambda = function(value) {0}
  )
)

dyadSampling_fit <-
  R6::R6Class(classname = "dyadSampling_fit",
  inherit = networkSamplingDyads_fit,
  private = list(
    card_D_o = NULL, # number of observed dyads
    card_D_m = NULL  # number of missing dyads
  ),
  public = list(
    initialize = function(sampledNetwork, ...) {
      super$initialize(sampledNetwork, "dyad")
      private$card_D_o <- length(sampledNetwork$observedDyads)
      private$card_D_m <- length(sampledNetwork$missingDyads )
      private$psi      <- check_boundaries(private$card_D_o / (private$card_D_m + private$card_D_o))
    }
  ),
  active = list(
    vExpec = function(value) {
      res <- private$card_D_o * log(private$psi) + private$card_D_m * log(1 - private$psi)
      res
    }
  )
)

dyadSampling_fit_covariates <-
  R6::R6Class(classname = "dyadSampling_fit_covariates",
  inherit = networkSamplingDyads_fit,
  private = list(
    D_obs = NULL, # observed dyads
    rho   = NULL # matrix of predicted probabilities of observation
  ),
  public = list(
    initialize = function(sampledNetwork, ...) {
      super$initialize(sampledNetwork, "dyad")
      X <- apply(sampledNetwork$covariatesArray, 3, as.vector)
      y <- 1*as.vector(!sampledNetwork$NAs)
      glm_out     <- glm.fit(X, y, family = binomial())
      private$psi <- coefficients(glm_out)
      private$rho <- fitted(glm_out)
      private$D_obs <- sampledNetwork$D_obs
    }
  ),
  active = list(
    prob_obs = function(value) {private$rho},
    vExpec = function(value) {
      res <- sum(log(1 - private$rho[private$D_miss])) + sum(log(private$rho[private$D_obs]))
      res
    }
  )
)

nodeSampling_fit <-
  R6::R6Class(classname = "nodeSampling_fit",
  inherit = networkSamplingNodes_fit,
  private = list(
    card_N_o = NULL, # number of observed nodes
    card_N_m = NULL  # number of missing nodes
  ),
  public = list(
    initialize = function(sampledNetwork, ...) {
      super$initialize(sampledNetwork, "node")
      private$card_N_o <- sum( sampledNetwork$observedNodes)
      private$card_N_m <- sum(!sampledNetwork$observedNodes)
      private$psi <- private$card_N_o / (private$card_N_o + private$card_N_m)
    }
  ),
  active = list(
    vExpec = function() {
      res <- private$card_N_o * log(private$psi) + private$card_N_m * log(1 - private$psi)
      res
    }
  )
)

nodeSampling_fit_covariates <-
  R6::R6Class(classname = "nodeSampling_fit_covariates",
  inherit = networkSamplingNodes_fit,
  private = list(
    rho = NULL # vector of predicted probabilities of observation
  ),
  public = list(
    initialize = function(sampledNetwork, ...) {
      super$initialize(sampledNetwork, "node")
      y <- 1 * (sampledNetwork$observedNodes)
      X <- sampledNetwork$covariatesMatrix
      glm_out     <- glm.fit(X, y, family = binomial())
      private$psi <- coefficients(glm_out)
      private$rho <- fitted(glm_out)
    }
  ),
  active = list(
    prob_obs = function(value) {private$rho},
    vExpec = function() {
      res <- sum(log(private$rho[private$N_obs])) + sum(log(1-private$rho[!private$N_obs]))
      res
    }
  )
)

doubleStandardSampling_fit <-
  R6::R6Class(classname = "doubleStandardSampling_fit",
  inherit = networkSamplingDyads_fit,
  private = list(
    So     = NULL, ## statistics only requiring the observed part of the network
    So.bar = NULL, ## can be computed once for all during the initialization
    Sm     = NULL, ## these ones will be updated during the optimization
    Sm.bar = NULL
  ),
  public = list(
    initialize = function(sampledNetwork, ...) {
      super$initialize(sampledNetwork, "double_standard")
      private$So      <- sum(    sampledNetwork$adjacencyMatrix[sampledNetwork$observedDyads])
      private$So.bar  <- sum(1 - sampledNetwork$adjacencyMatrix[sampledNetwork$observedDyads])
      ## can we do better than that?
      imputedNet      <- matrix(mean(sampledNetwork$adjacencyMatrix, na.rm = TRUE), sampledNetwork$nNodes, sampledNetwork$nNodes)
      self$update_parameters(imputedNet)
    },
    update_parameters = function(imputedNet, ...) {
      private$Sm     <- sum(    imputedNet[private$D_miss])
      private$Sm.bar <- sum(1 - imputedNet[private$D_miss])
      private$psi    <- c(private$So.bar / (private$So.bar + private$Sm.bar), private$So / (private$So + private$Sm))
    },
    update_imputation = function(Z, pi) {
      nu <- check_boundaries(logistic(log((1 - private$psi[2]) / (1 - private$psi[1])) + Z %*% log(pi/(1 - pi)) %*% t(Z)))
      nu
    }
  ),
  active = list(
    vExpec = function(value) {
      res <- log(private$psi[2]) * private$So + log(private$psi[1]) * private$So.bar +
        log(1 - private$psi[2]) * private$Sm + log(1 - private$psi[1]) * private$Sm.bar
      res
    }
  )
)

blockDyadSampling_fit <-
  R6::R6Class(classname = "blockDyadSampling_fit",
  inherit  = networkSamplingDyads_fit,
  private  = list(
    prob     = NULL,  ## for calculation of the log-likelihood
    NAs      = NULL,  ## localisation of NAs
    R        = NULL,  ## sampling matrix
    directed = NULL   ##
  ),
  public = list(
    initialize = function(sampledNetwork, blockInit) {
      super$initialize(sampledNetwork, "block_dyad")
      private$NAs      <- sampledNetwork$NAs
      private$R        <- sampledNetwork$samplingMatrix
      private$directed <- sampledNetwork$is_directed
      imputedNet       <- matrix(mean(sampledNetwork$adjacencyMatrix, na.rm = TRUE), sampledNetwork$nNodes, sampledNetwork$nNodes)
      self$update_parameters(imputedNet, blockInit)
    },
    update_parameters = function(imputedNet, Z) {
      private$psi    <- check_boundaries((t(Z) %*% private$R %*% Z) / (t(Z) %*% (1 - diag(nrow(imputedNet))) %*% Z))
      private$prob   <- check_boundaries(Z %*% private$psi %*% t(Z))
    }
  ),
  active = list(
    vExpec = function(value) {
      factor       <- ifelse(private$directed, 1, .5)
      sampMat      <- private$R ; diag(sampMat) <- 0
      sampMat_bar  <- 1 - private$R ; diag(sampMat_bar) <- 0
      res          <- factor * sum(sampMat * log(private$prob) + sampMat_bar *  log(1 - private$prob))
      res
    }
  )
)

blockSampling_fit <-
  R6::R6Class(classname = "blockSampling_fit",
  inherit = networkSamplingNodes_fit,
  private = list(
    So     = NULL, ## sum_(i in Nobs ) Z_iq
    Sm     = NULL  ## sum_(i in Nmiss) Z_iq
  ),
  public = list(
    initialize = function(sampledNetwork, blockInit) {
      super$initialize(sampledNetwork, "block")
      self$update_parameters(NA, blockInit)
    },
    update_parameters = function(imputedNet, Z) {
      private$So <- colSums(Z[ private$N_obs, , drop = FALSE])
      private$Sm <- colSums(Z[!private$N_obs, , drop = FALSE])
      private$psi <- check_boundaries(private$So / (private$So + private$Sm))
    }
  ),
  active = list(
    vExpec = function() {
      res <- c(crossprod(private$So, log(private$psi)) +  crossprod(private$Sm, log(1 - private$psi)))
      res
    },
    log_lambda = function(value) {
      res <- sapply(private$psi, function(psi) ifelse(private$N_obs, log(psi), log(1 - psi)))
      res
    }
  )
)

degreeSampling_fit <-
  R6::R6Class(classname = "degreeSampling_fit",
  inherit = networkSamplingNodes_fit,
  private = list(
    NAs = NULL,
    Dij = NULL, # stat required to perform imputation
    ksi = NULL, # additional variational parameter for approximation in the logistic function
    D   = NULL  # estimation of the degrees on the observed part of the network
  ),
  public = list(
    initialize = function(sampledNetwork, blockInit, connectInit) {
      super$initialize(sampledNetwork, "degree")

      private$NAs <- sampledNetwork$NAs
      ## will remain the same
      private$D <- rowSums(sampledNetwork$adjacencyMatrix, na.rm = TRUE)

      ## will fluctuate along the algorithm
      private$psi <- coefficients(glm(1*(private$N_obs) ~ private$D, family = binomial(link = "logit")))

      imputedNet <- blockInit %*% connectInit %*% t(blockInit)
      imputedNet[!private$NAs] <- sampledNetwork$adjacencyMatrix[!private$NAs]
      self$update_parameters(imputedNet)
    },
    update_parameters = function(imputedNet, ...) {
      private$D  <- rowSums(imputedNet)
      nu <- imputedNet
      nu[!private$NAs] <- NA
      D2 <- rowSums(nu * (1 - nu), na.rm = TRUE) + private$D^2
      private$ksi <- check_boundaries(sqrt( private$psi[1]^2 + private$psi[2]^2 * D2  + 2 * private$psi[1] * private$psi[2] * private$D))
      C <- .5 * nrow(imputedNet) - sum(!private$N_obs)
      s_hksi   <- sum(h(private$ksi))
      s_hksiD  <- sum(h(private$ksi) * private$D)
      s_hksiD2 <- sum(h(private$ksi) * D2)
      b <- (2 * C * s_hksiD  - (.5 * sum(private$D) - sum(private$D[!private$N_obs])) * s_hksi) / (2 * s_hksiD2 * s_hksi - (2 * s_hksiD)^2)
      a <- -(b * s_hksiD + C) / s_hksi
      private$psi <- c(a, b)

      ## update stat required to perform imputation
      nu[!private$NAs] <- 0
      private$Dij <- matrix(private$D, nrow(imputedNet), ncol(imputedNet)) - nu
    },
    update_imputation = function(Z, pi) {
      C <- 2 * h(private$ksi) * (private$psi[1] * private$psi[2] + private$psi[2]^2 * (1 + private$Dij))
      nu <- check_boundaries((Z %*% log(pi/(1 - pi)) %*% t(Z) - private$psi[2] + C + t(C) ))
      nu
    }
  ),
  active = list(
    vExpec = function() {
      prob <-  check_boundaries(logistic(private$psi[1] + private$psi[2] * private$D))
      res  <-  sum(private$N_obs * log(prob)) + sum( (!private$N_obs) * log(1 - prob))
      res
    }
  )
)
