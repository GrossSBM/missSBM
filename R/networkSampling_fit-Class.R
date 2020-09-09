#' Virtual class used to define a family of networkSamplingDyads_fit
#' @include networkSampling-Class.R
#' @import R6
networkSamplingDyads_fit <-
  R6::R6Class(classname = "networkSamplingDyads_fit",
  inherit = networkSampling,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    card_D = NULL, # number of possible dyads in the network
    D_miss = NULL  # where are the missing dyads
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description constructor for networkSampling_fit
    #' @param sampledNetwork a object with class sampledNetwork representing the observed data with possibly missing entries
    #' @param name a character for the name of sampling to fit on the sampledNetwork
    initialize = function(sampledNetwork, name) {
      private$name    <- name
      private$D_miss  <- sampledNetwork$missingDyads
      private$card_D  <- sampledNetwork$nDyads
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
    #' @param PI the matrix of inter/intra class probability of connection
    update_imputation = function(PI) {PI} ## good for MCAR on node, dyads and NMAR with blocks
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## ACTIVE BINDING
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field penalty double, value of the penalty term in vICL
    penalty = function(value) {log(private$card_D) * self$df},
    #' @field log_lambda double, term for adjusting the imputation stepn which depends on the type of sampling
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
    #' @description constructor
    #' @param sampledNetwork a object with class sampledNetwork representing the observed data with possibly missing entries
    #' @param name a character for the name of sampling to fit on the sampledNetwork
    initialize = function(sampledNetwork, name) {
      private$name   <- name
      private$N_obs  <- sampledNetwork$observedNodes
      private$card_N <- sampledNetwork$nNodes
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
    #' @param PI the matrix of inter/intra class probability of connection
    update_imputation = function(PI) {PI} ## good for MCAR on node, dyads and NMAR with blocks
  ),
  active = list(
    #' @field penalty double, value of the penalty term in vICL
    penalty = function(value) {log(private$card_N) * self$df},
    #' @field log_lambda double, term for adjusting the imputation stepn which depends on the type of sampling
    log_lambda = function(value) {0}
  )
)

#' Class for fitting a dyad sampling
dyadSampling_fit <-
  R6::R6Class(classname = "dyadSampling_fit",
  inherit = networkSamplingDyads_fit,
  private = list(
    card_D_o = NULL, # number of observed dyads
    card_D_m = NULL  # number of missing dyads
  ),
  public = list(
    #' @description constructor
    #' @param sampledNetwork a object with class sampledNetwork representing the observed data with possibly missing entries
    #' @param ... used for compatibility
    initialize = function(sampledNetwork, ...) {
      super$initialize(sampledNetwork, "dyad")
      private$card_D_o <- length(sampledNetwork$observedDyads)
      private$card_D_m <- length(sampledNetwork$missingDyads )
      private$psi      <- check_boundaries(private$card_D_o / (private$card_D_m + private$card_D_o))
      private$rho      <- rep(private$psi, private$card_D)
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
  private = list(
    D_obs = NULL, # observed dyads
    rho   = NULL  # matrix of predicted probabilities of observation
  ),
  public = list(
    #' @description constructor
    #' @param sampledNetwork a object with class sampledNetwork representing the observed data with possibly missing entries
    #' @param ... used for compatibility
    initialize = function(sampledNetwork, ...) {
      super$initialize(sampledNetwork, "covar-dyad")
      X <- cbind(1, apply(sampledNetwork$covarArray, 3, as.vector))
      y <- 1 * as.vector(!sampledNetwork$NAs)
      glm_out       <- glm.fit(X, y, family = binomial())
      private$psi   <- coefficients(glm_out)
      private$rho   <- fitted(glm_out)
      private$D_obs <- sampledNetwork$D_obs
    }
  ),
  active = list(
    #' @field prob_obs sampling rate
    prob_obs = function(value) {private$rho},
    #' @field vExpec variational expectation of the sampling
    vExpec = function(value) {
      res <- sum(log(1 - private$rho[private$D_miss])) + sum(log(private$rho[private$D_obs]))
      res
    }
  )
)


#' Class for fitting a node sampling
nodeSampling_fit <-
  R6::R6Class(classname = "nodeSampling_fit",
  inherit = networkSamplingNodes_fit,
  private = list(
    card_N_o = NULL, # number of observed nodes
    card_N_m = NULL  # number of missing nodes
  ),
  public = list(
    #' @description constructor
    #' @param sampledNetwork a object with class sampledNetwork representing the observed data with possibly missing entries
    #' @param ... used for compatibility
    initialize = function(sampledNetwork, ...) {
      super$initialize(sampledNetwork, "node")
      private$card_N_o <- sum( sampledNetwork$observedNodes)
      private$card_N_m <- sum(!sampledNetwork$observedNodes)
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
  #' @param sampledNetwork a object with class sampledNetwork representing the observed data with possibly missing entries
  #' @param ... used for compatibility
  public = list(
    initialize = function(sampledNetwork, ...) {
      super$initialize(sampledNetwork, "covar-node")
      y <- 1 * (sampledNetwork$observedNodes)
      glm_out     <- glm.fit(cbind(1, sampledNetwork$covarMatrix), y, family = binomial())
      private$psi <- coefficients(glm_out)
      private$rho <- fitted(glm_out)
    }
  ),
  active = list(
    #' @field prob_obs sampling rate
    prob_obs = function(value) {private$rho},
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
    #' @param sampledNetwork a object with class sampledNetwork representing the observed data with possibly missing entries
    #' @param ... used for compatibility
    initialize = function(sampledNetwork, ...) {
      super$initialize(sampledNetwork, "double-standard")
      private$So      <- sum(    sampledNetwork$adjacencyMatrix[sampledNetwork$observedDyads])
      private$So.bar  <- sum(1 - sampledNetwork$adjacencyMatrix[sampledNetwork$observedDyads])
      ## can we do better than that?
      imputedNet      <- matrix(mean(sampledNetwork$adjacencyMatrix, na.rm = TRUE), sampledNetwork$nNodes, sampledNetwork$nNodes)
      self$update_parameters(imputedNet)
    },
    #' @description a method to update the estimation of the parameters. By default, nothing to do (corresponds to MAR sampling)
    #' @param imputedNet an adjacency matrix where missing values have been imputed
    #' @param ... use for compatibility
    update_parameters = function(imputedNet, ...) {
      private$Sm     <- sum(    imputedNet[private$D_miss])
      private$Sm.bar <- sum(1 - imputedNet[private$D_miss])
      private$psi    <- c(private$So.bar / (private$So.bar + private$Sm.bar), private$So / (private$So + private$Sm))
    },
    #' @description a method to update the imputation of the missing entries.
    #' @param PI the matrix of inter/intra class probability of connection
    update_imputation = function(PI) {
      nu <- check_boundaries(logistic(log((1 - private$psi[2]) / (1 - private$psi[1])) + logit(PI) ))
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
    prob     = NULL,  ## for calculation of the log-likelihood
    NAs      = NULL,  ## localisation of NAs
    R        = NULL,  ## sampling matrix
    directed = NULL   ##
  ),
  public = list(
    #' @description constructor
    #' @param sampledNetwork a object with class sampledNetwork representing the observed data with possibly missing entries
    #' @param blockInit n x Q matrix of initial block indicators
    initialize = function(sampledNetwork, blockInit) {
      super$initialize(sampledNetwork, "block-dyad")
      private$NAs      <- sampledNetwork$NAs
      private$R        <- sampledNetwork$samplingMatrix
      private$directed <- sampledNetwork$is_directed
      imputedNet       <- matrix(mean(sampledNetwork$adjacencyMatrix, na.rm = TRUE), sampledNetwork$nNodes, sampledNetwork$nNodes)
      self$update_parameters(imputedNet, blockInit)
    },
    #' @description a method to update the estimation of the parameters. By default, nothing to do (corresponds to MAR sampling)
    #' @param imputedNet an adjacency matrix where missing values have been imputed
    #' @param Z indicator of blocks
    update_parameters = function(imputedNet, Z) {
      private$psi    <- check_boundaries((t(Z) %*% private$R %*% Z) / (t(Z) %*% (1 - diag(nrow(imputedNet))) %*% Z))
      private$prob   <- check_boundaries(Z %*% private$psi %*% t(Z))
    }
  ),
  active = list(
    #' @field vExpec variational expectation of the sampling
    vExpec = function(value) {
      factor       <- ifelse(private$directed, 1, .5)
      sampMat      <- private$R ; diag(sampMat) <- 0
      sampMat_bar  <- 1 - private$R ; diag(sampMat_bar) <- 0
      res          <- factor * sum(sampMat * log(private$prob) + sampMat_bar *  log(1 - private$prob))
      res
    }
  )
)

#' Class for fitting a block-node sampling
blockSampling_fit <-
  R6::R6Class(classname = "blockSampling_fit",
  inherit = networkSamplingNodes_fit,
  private = list(
    So     = NULL, ## sum_(i in Nobs ) Z_iq
    Sm     = NULL  ## sum_(i in Nmiss) Z_iq
  ),
  public = list(
    #' @description constructor
    #' @param sampledNetwork a object with class sampledNetwork representing the observed data with possibly missing entries
    #' @param blockInit n x Q matrix of initial block indicators
    initialize = function(sampledNetwork, blockInit) {
      super$initialize(sampledNetwork, "block-node")
      self$update_parameters(NA, blockInit)
    },
    #' @description a method to update the estimation of the parameters. By default, nothing to do (corresponds to MAR sampling)
    #' @param imputedNet an adjacency matrix where missing values have been imputed
    #' @param Z indicator of blocks
    update_parameters = function(imputedNet, Z) {
      private$So <- colSums(Z[ private$N_obs, , drop = FALSE])
      private$Sm <- colSums(Z[!private$N_obs, , drop = FALSE])
      private$psi <- check_boundaries(private$So / (private$So + private$Sm))
    }
  ),
  active = list(
    #' @field vExpec variational expectation of the sampling
    vExpec = function() {
      res <- c(crossprod(private$So, log(private$psi)) +  crossprod(private$Sm, log(1 - private$psi)))
      res
    },
    #' @field log_lambda double, term for adjusting the imputation stepn which depends on the type of sampling
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
    #' @param sampledNetwork a object with class sampledNetwork representing the observed data with possibly missing entries
    #' @param blockInit n x Q matrix of initial block indicators
    #' @param connectInit Q x Q matrix of initial block probabilities of connection
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
    update_imputation = function(PI) {
      C <- 2 * h(private$ksi) * (private$psi[1] * private$psi[2] + private$psi[2]^2 * (1 + private$Dij))
      nu <- check_boundaries((logit(PI) - private$psi[2] + C + t(C) ))
      nu
    }
  ),
  active = list(
    #' @field vExpec variational expectation of the sampling
    vExpec = function() {
      prob <-  check_boundaries(logistic(private$psi[1] + private$psi[2] * private$D))
      res  <-  sum(private$N_obs * log(prob)) + sum( (!private$N_obs) * log(1 - prob))
      res
    }
  )
)
