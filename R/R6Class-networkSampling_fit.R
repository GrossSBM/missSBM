#' Virtual class used to define a family of networkSamplingDyads_fit
#' @include R6Class-networkSampling.R
#' @import R6
networkSamplingDyads_fit <-
  R6::R6Class(classname = "networkSamplingDyads_fit",
  inherit = networkSampling,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    card_D   = NULL, # number of possible dyads in the network
    card_D_o = NULL, # number of observed dyads in the network
    card_D_m = NULL, # number of missing dyads in the network
    D_miss   = NULL,  # where are the missing dyads
    D_obs    = NULL  # observed dyads
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description constructor for networkSampling_fit
    #' @param partlyObservedNetwork a object with class partlyObservedNetwork representing the observed data with possibly missing entries
    #' @param name a character for the name of sampling to fit on the partlyObservedNetwork
    initialize = function(partlyObservedNetwork, name) {
      private$name     <- name
      private$D_miss   <- partlyObservedNetwork$missingDyads
      private$D_obs    <- partlyObservedNetwork$observedDyads
      private$card_D   <- partlyObservedNetwork$nbDyads
      private$card_D_o <- sum(partlyObservedNetwork$samplingMatrix)
      private$card_D_m <- sum(partlyObservedNetwork$samplingMatrixBar)
    },
    #' @description show method
    show = function() {
      super$show()
      cat("  $penalty, $vExpec\n")
    },
    #' @description a method to update the estimation of the parameters. By default, nothing to do (corresponds to MAR sampling)
    #' @param ... use for compatibility
    update_parameters = function(...) {invisible(NULL)},
    #' @description a method to update the imputation of the missing entries.
    #' @param nu the matrix of (uncorrected) imputation for missing entries
    update_imputation = function(nu) {nu} # good for MCAR on node, dyads and NMAR with blocks
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## ACTIVE BINDING
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field penalty double, value of the penalty term in ICL
    penalty = function(value) {log(private$card_D) * self$df},
    #' @field log_lambda double, term for adjusting the imputation step which depends on the type of sampling
    log_lambda = function(value) {0}
  )
)

#' Virtual class used to define a family of networkSamplingNodes_fit
#' @include R6Class-networkSampling.R
#' @import R6
networkSamplingNodes_fit <-
  R6::R6Class(classname = "networkSamplingNodes_fit",
  inherit = networkSampling,
  private = list(
    card_N   = NULL, # number of nodes in the network
    card_N_m = NULL, # number of observed nodes
    card_N_o = NULL, # number of missing nodes
    N_obs    = NULL  # boolean for observed nodes
  ),
  public = list(
    #' @description constructor
    #' @param partlyObservedNetwork a object with class partlyObservedNetwork representing the observed data with possibly missing entries
    #' @param name a character for the name of sampling to fit on the partlyObservedNetwork
    initialize = function(partlyObservedNetwork, name) {
      private$name     <- name
      private$N_obs    <- partlyObservedNetwork$observedNodes
      private$card_N   <- partlyObservedNetwork$nbNodes
      private$card_N_o <- sum( partlyObservedNetwork$observedNodes)
      private$card_N_m <- sum(!partlyObservedNetwork$observedNodes)
    },
    #' @description show method
    show = function() {
      super$show()
      cat("  $penalty, $vExpec\n")
    },
    #' @description a method to update the estimation of the parameters. By default, nothing to do (corresponds to MAR sampling)
    #' @param ... use for compatibility
    update_parameters = function(...) {invisible(NULL)},
    #' @description a method to update the imputation of the missing entries.
    #' @param nu the matrix of (uncorrected) imputation for missing entries
    update_imputation = function(nu) {nu} ## good for MCAR on node, dyads and NMAR with blocks
  ),
  active = list(
    #' @field penalty double, value of the penalty term in ICL
    penalty = function(value) {log(private$card_N) * self$df},
    #' @field log_lambda double, term for adjusting the imputation step which depends on the type of sampling
    log_lambda = function(value) {0}
  )
)

#' Class for fitting a dyad sampling
dyadSampling_fit <-
  R6::R6Class(classname = "dyadSampling_fit",
  inherit = networkSamplingDyads_fit,
  public = list(
    #' @description constructor
    #' @param partlyObservedNetwork a object with class partlyObservedNetwork representing the observed data with possibly missing entries
    #' @param ... used for compatibility
    initialize = function(partlyObservedNetwork, ...) {
      super$initialize(partlyObservedNetwork, "dyad")
      private$psi <- check_boundaries(private$card_D_o / (private$card_D_m + private$card_D_o))
    }
  ),
  active = list(
    #' @field vExpec variational expectation of the sampling
    vExpec = function(value) {
      res <- private$card_D_o * log(private$psi + 1 * (private$psi == 0)) +
        private$card_D_m * log(1 - private$psi + 1 * (private$psi == 1))
      res
    }
  )
)

#' Class for fitting a dyad sampling with covariates
covarDyadSampling_fit <-
  R6::R6Class(classname = "covarDyadSampling_fit",
  inherit = networkSamplingDyads_fit,
  public = list(
    #' @description constructor
    #' @param partlyObservedNetwork a object with class partlyObservedNetwork representing the observed data with possibly missing entries
    #' @param ... used for compatibility
    initialize = function(partlyObservedNetwork, ...) {
      super$initialize(partlyObservedNetwork, "covar-dyad")
      dyads <- rbind(partlyObservedNetwork$observedDyads, partlyObservedNetwork$missingDyads)
      X <- cbind(1, apply(partlyObservedNetwork$covarArray, 3, function(x) x[dyads]))
      y <- partlyObservedNetwork$samplingMatrix[dyads]
      glm_out     <- glm.fit(X, y, family = binomial())
      private$psi <- coefficients(glm_out)
      y_hat <- fitted(glm_out)
      private$rho <- list(obs = y_hat[y == 1], miss = y_hat[y == 0])

    }
  ),
  active = list(
    #' @field vExpec variational expectation of the sampling
    vExpec = function(value) {
      res <- sum(log(1 - private$rho$miss)) + sum(log(private$rho$obs))
      res
    }
  )
)

#' Class for fitting a node sampling
nodeSampling_fit <-
  R6::R6Class(classname = "nodeSampling_fit",
  inherit = networkSamplingNodes_fit,
  public = list(
    #' @description constructor
    #' @param partlyObservedNetwork a object with class partlyObservedNetwork representing the observed data with possibly missing entries
    #' @param ... used for compatibility
    initialize = function(partlyObservedNetwork, ...) {
      super$initialize(partlyObservedNetwork, "node")
      private$psi <- private$card_N_o / (private$card_N_o + private$card_N_m)
    }
  ),
  active = list(
    #' @field vExpec variational expectation of the sampling
    vExpec = function() {
      res <- private$card_N_o * log(private$psi + 1 * (private$psi == 0)) +
        private$card_N_m * log(1 - private$psi + 1 * (private$psi == 1) )
      res
    }
  )
)

#' Class for fitting a node-centered sampling with covariate
covarNodeSampling_fit <-
  R6::R6Class(classname = "covarNodeSampling_fit",
  inherit = networkSamplingNodes_fit,
  private = list(
    rho = NULL # vector of predicted probabilities of observation
  ),
  #' @description constructor
  #' @param partlyObservedNetwork a object with class partlyObservedNetwork representing the observed data with possibly missing entries
  #' @param ... used for compatibility
  public = list(
    initialize = function(partlyObservedNetwork, ...) {
      super$initialize(partlyObservedNetwork, "covar-node")
      y <- 1 * (partlyObservedNetwork$observedNodes)
      glm_out     <- glm.fit(cbind(1, partlyObservedNetwork$covarMatrix), y, family = binomial())
      private$psi <- coefficients(glm_out)
      private$rho <- fitted(glm_out)
    }
  ),
  active = list(
    #' @field vExpec variational expectation of the sampling
    vExpec = function() {
      res <- sum(log(private$rho[private$N_obs])) + sum(log(1 - private$rho[!private$N_obs]))
      res
    }
  )
)

#' Class for fitting a double-standard sampling
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
    #' @description constructor
    #' @param partlyObservedNetwork a object with class partlyObservedNetwork representing the observed data with possibly missing entries
    #' @param ... used for compatibility
    initialize = function(partlyObservedNetwork, ...) {
      super$initialize(partlyObservedNetwork, "double-standard")
      private$So      <- sum(partlyObservedNetwork$networkData)
      private$So.bar  <- private$card_D_o - private$So
      self$update_parameters(partlyObservedNetwork$imputation("average") * partlyObservedNetwork$samplingMatrixBar)
    },
    #' @description a method to update the estimation of the parameters. By default, nothing to do (corresponds to MAR sampling)
    #' @param nu an adjacency matrix with imputed values (only)
    #' @param ... use for compatibility
    update_parameters = function(nu, ...) {
      private$Sm     <- sum(nu)
      private$Sm.bar <- private$card_D_m - private$Sm
      private$psi    <- c(private$So.bar / (private$So.bar + private$Sm.bar), private$So / (private$So + private$Sm))
    },
    #' @description a method to update the imputation of the missing entries.
    #' @param nu the matrix of (uncorrected) imputation for missing entries
    update_imputation = function(nu) {
      nu@x <- .logistic(log((1 - private$psi[2]) / (1 - private$psi[1])) + .logit(nu@x) )
      nu
    }
  ),
  active = list(
    #' @field vExpec variational expectation of the sampling
    vExpec = function(value) {
      res <- log(private$psi[2]) * private$So + log(private$psi[1]) * private$So.bar +
        log(1 - private$psi[2]) * private$Sm + log(1 - private$psi[1]) * private$Sm.bar
      res
    }
  )
)

#' Class for fitting a block-dyad sampling
blockDyadSampling_fit <-
  R6::R6Class(classname = "blockDyadSampling_fit",
  inherit  = networkSamplingDyads_fit,
  private  = list(
    R        = NULL,  ## sampling matrix
    S        = NULL,
    Z        = NULL,
    directed = NULL
  ),
  public = list(
    #' @description constructor
    #' @param partlyObservedNetwork a object with class partlyObservedNetwork representing the observed data with possibly missing entries
    #' @param blockInit n x Q matrix of initial block indicators
    initialize = function(partlyObservedNetwork, blockInit) {
      super$initialize(partlyObservedNetwork, "block-dyad")
      private$R <- partlyObservedNetwork$samplingMatrix
      private$S <- partlyObservedNetwork$samplingMatrixBar
      private$directed <- partlyObservedNetwork$is_directed
      self$update_parameters(NA, blockInit)
    },
    #' @description a method to update the estimation of the parameters. By default, nothing to do (corresponds to MAR sampling)
    #' @param nu the matrix of (uncorrected) imputation for missing entries
    #' @param Z probabilities of block memberships
    update_parameters = function(nu, Z) {
      private$Z <- Z
      ZtRZ <- as.matrix(t(Z) %*% private$R %*% Z)
      Zbar <- colSums(Z)
      if(private$directed) {
        private$psi <- check_boundaries(ZtRZ / ( Zbar %o% Zbar - Zbar ))
      } else {
        private$psi <- check_boundaries(( ZtRZ + t(ZtRZ) ) / ( Zbar %o% Zbar - Zbar ))
      }
    }
  ),
  active = list(
    #' @field vExpec variational expectation of the sampling
    vExpec = function(value) {
      rho <- check_boundaries(private$Z %*% private$psi %*% t(private$Z))
      sum(private$R * log(rho) + private$S *  log(1 - rho))
    },
    #' @field log_lambda matrix, term for adjusting the imputation step which depends on the type of sampling
    log_lambda = function(value) {
      res <- private$R %*% private$Z %*% t(log(private$psi)) + private$S %*% private$Z %*% t(log(1 - private$psi)) +
        t(private$R) %*% private$Z %*% log(private$psi) + t(private$S) %*% private$Z %*% log(1 - private$psi)
      res
    }
  )
)

#' Class for fitting a block-node sampling
blockNodeSampling_fit <-
  R6::R6Class(classname = "blockNodeSampling_fit",
  inherit = networkSamplingNodes_fit,
  private = list(
    So     = NULL, ## sum_(i in Nobs ) Z_iq
    Sm     = NULL  ## sum_(i in Nmiss) Z_iq
  ),
  public = list(
    #' @description constructor
    #' @param partlyObservedNetwork a object with class partlyObservedNetwork representing the observed data with possibly missing entries
    #' @param blockInit n x Q matrix of initial block indicators
    initialize = function(partlyObservedNetwork, blockInit) {
      super$initialize(partlyObservedNetwork, "block-node")
      self$update_parameters(NA, blockInit)
    },
    #' @description a method to update the estimation of the parameters. By default, nothing to do (corresponds to MAR sampling)
    #' @param imputedNet an adjacency matrix where missing values have been imputed
    #' @param Z indicator of blocks
    update_parameters = function(imputedNet, Z) {
      private$So  <- colSums(Z[ private$N_obs, , drop = FALSE])
      private$Sm  <- colSums(Z[!private$N_obs, , drop = FALSE])
      private$psi <- check_boundaries(private$So / (private$So + private$Sm))
    }
  ),
  active = list(
    #' @field vExpec variational expectation of the sampling
    vExpec = function() {
      res <- c(crossprod(private$So, log(private$psi)) +  crossprod(private$Sm, log(1 - private$psi)))
      res
    },
    #' @field log_lambda double, term for adjusting the imputation step which depends on the type of sampling
    log_lambda = function(value) {
      res <- sapply(private$psi, function(psi) ifelse(private$N_obs, log(psi), log(1 - psi)))
      res
    }
  )
)

#' Class for fitting a degree sampling
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
    #' @description constructor
    #' @param partlyObservedNetwork a object with class partlyObservedNetwork representing the observed data with possibly missing entries
    #' @param blockInit n x Q matrix of initial block indicators
    #' @param connectInit Q x Q matrix of initial block probabilities of connection
    initialize = function(partlyObservedNetwork, blockInit, connectInit) {
      super$initialize(partlyObservedNetwork, "degree")

      private$NAs <- as.matrix(partlyObservedNetwork$NAs)
      ## the node degrees - will remain the same
      private$D <- rowSums(partlyObservedNetwork$networkData, na.rm = TRUE)

      ## will fluctuate along the algorithm
      private$psi <- coefficients(glm(1*(private$N_obs) ~ private$D, family = binomial(link = "logit")))

      imputedNet <- blockInit %*% connectInit %*% t(blockInit)
      imputedNet[!private$NAs] <- partlyObservedNetwork$networkData[!private$NAs]
      self$update_parameters(imputedNet)
    },
    #' @description a method to update the estimation of the parameters. By default, nothing to do (corresponds to MAR sampling)
    #' @param imputedNet an adjacency matrix where missing values have been imputed
    #' @param ... used for compatibility
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
    #' @description a method to update the imputation of the missing entries.
    #' @param PI the matrix of inter/intra class probability of connection
    #' @param ... use for compatibility
    update_imputation = function(PI,...) {
      C <- 2 * h(private$ksi) * (private$psi[1] * private$psi[2] + private$psi[2]^2 * (1 + private$Dij))
      nu <- check_boundaries((.logit(PI) - private$psi[2] + C + t(C) ))
      nu
    }
  ),
  active = list(
    #' @field vExpec variational expectation of the sampling
    vExpec = function() {
      prob <-  check_boundaries(.logistic(private$psi[1] + private$psi[2] * private$D))
      res  <-  sum(private$N_obs * log(prob)) + sum( (!private$N_obs) * log(1 - prob))
      res
    }
  )
)
