#' This internal class is designed to adjust a binary Stochastic Block Model in the context of missSBM.
#'
#' It is not designed not be call by the user
#'
#' @import R6
SimpleSBM_fit <-
R6::R6Class(classname = "SimpleSBM_fit",
  inherit = sbm::SimpleSBM,
  private = list(
    variant = NULL, # model variant
    R       = NULL, # the sampling matrix (sparse encoding)
    S       = NULL, # the "anti" sampling matrix (sparse encoding)
    M_step  = NULL, # pointing to the appropriate M_step function
    E_step  = NULL, # pointing to the appropriate E_step function
    vLL_complete = NULL # pointing to the appropriate Expected complete LL function
  ),
  public = list(
    #' @description constructor for simpleSBM_fit for missSBM purpose
    #' @param networkData a structure to store network under missing data condition: either a matrix possibly with NA, or a missSBM:::partlyObservedNetwork
    #' @param clusterInit Initial clustering: a vector with size \code{ncol(adjacencyMatrix)}, providing a user-defined clustering with \code{nbBlocks} levels.
    #' @param covarList An optional list with M entries (the M covariates).
    initialize = function(networkData, clusterInit, covarList = list()) {

      ## networkData
      if (inherits(networkData, "matrix")) networkData <- partlyObservedNetwork$new(networkData, covarList)
      stopifnot(inherits(networkData, "partlyObservedNetwork"))

      ## Initial Clustering
      private$Z <- clustering_indicator(clusterInit)

      # Basic fields initialization by call to super constructor
### TODO: determine the model according to the adjacency matrix
      model <- "bernoulli"
      super$initialize(model        = model,
                       directed     = networkData$is_directed,
                       nbNodes      = networkData$nbNodes,
                       blockProp    = vector("numeric", ncol(private$Z)) + .Machine$double.eps,
### FIXME: this should be model-dependent
                       connectParam = list(mean = matrix(0, ncol(private$Z), ncol(private$Z))),
                       covarList    = covarList)

      # Storing data
      private$R <- networkData$samplingMatrix
      private$S <- networkData$samplingMatrixBar
      private$Y <- networkData$networkData

      ## point to the functions that performs E/M steps and compute the likelihood
      private$variant <-
        paste(model, ifelse(self$directed, "directed", "undirected"),
          ifelse(self$nbCovariates>0, "covariates", "nocovariate"), sep="_")
      private$M_step       <- get(paste("M_step_sparse"      , model, ifelse(self$nbCovariates>0, "covariates", "nocovariate"), sep = "_"))
      private$E_step       <- get(paste("E_step_sparse"      , model, ifelse(self$nbCovariates>0, "covariates", "nocovariate"), sep = "_"))
      private$vLL_complete <- get(paste("vLL_complete_sparse", model, ifelse(self$nbCovariates>0, "covariates", "nocovariate"), sep = "_"))

### TODO:
###  - check if parameters are not already initialize
###  - specialize the initialization to each model (this should be done in sbm::SimpleSBM...)

      ## Initialize estimation of the model parameters
      private$theta$mean <- matrix(0.5, ncol(private$Z), ncol(private$Z))
      private$beta       <- numeric(self$nbCovariates)
      self$update_parameters()

      invisible(self)
    },
    #' @description method to perform estimation via variational EM
    #' @param threshold stop when an optimization step changes the objective function by less than threshold. Default is 1e-4.
    #' @param maxIter V-EM algorithm stops when the number of iteration exceeds maxIter. Default is 10
    #' @param fixPointIter number of fix-point iterations in the Variational E step. Default is 5.
    #' @param trace logical for verbosity. Default is \code{FALSE}.
    doVEM = function(threshold = 1e-4, maxIter = 10, fixPointIter = 3, trace = FALSE) {

      ## Initialization of quantities that monitor convergence
      delta_par <- vector("numeric", maxIter)
      delta_obj <- vector("numeric", maxIter)
      objective <- vector("numeric", maxIter)
      objective[1] <- self$loglik
      iterate <- 1; stop <- ifelse(self$nbBlocks > 1, FALSE, TRUE)

      ## Starting the variational EM algorithm
      if (trace) cat("\n Adjusting Variational EM for Stochastic Block Model\n")
      while (!stop) {
        iterate <- iterate + 1
        if (trace) cat(" iteration #:", iterate, "\r")

        theta_old <- private$theta # save old value of parameters to assess convergence

        # Variational E-Step
        for (i in seq.int(fixPointIter)) self$update_blocks()

        # M-step
        self$update_parameters()

        # Assess convergence
        objective[iterate] <- self$loglik
        delta_par[iterate] <- sqrt(sum((private$theta$mean - theta_old$mean)^2)) / sqrt(sum((theta_old$mean)^2))
        delta_obj[iterate] <- (objective[iterate] - objective[iterate-1]) / abs(objective[iterate])
        stop <- (iterate > maxIter) |  ((delta_par[iterate] < threshold) & (delta_obj[iterate] < threshold))
      }
      self$reorder()
      if (trace) cat("\n")
      res <- data.frame(delta_pararameters = delta_par[1:iterate], delta_objective = delta_obj[1:iterate],  elbo = objective[1:iterate])
      res
    },
    #' @description permute group labels by order of decreasing probability
    reorder = function(){
      o <- order(private$theta$mean %*% private$pi, decreasing = TRUE)
      private$pi <- private$pi[o]
      private$theta$mean <- private$theta$mean[o, o, drop = FALSE]
      private$Z <- private$Z[, o, drop = FALSE]
    }
  ),
  active = list(
    #' @field type the type of SBM (distribution of edges values, network type, presence of covariates)
    type = function(value) {private$variant},
    #' @field penalty double, value of the penalty term in ICL
    penalty  = function(value) {unname((self$nbConnectParam + self$nbCovariates) * log(self$nbDyads) + (self$nbBlocks-1) * log(self$nbNodes))},
    #' @field entropy double, value of the entropy due to the clustering distribution
    entropy  = function(value) {-sum(.xlogx(private$Z))},
    #' @field loglik double: approximation of the log-likelihood (variational lower bound) reached
    loglik = function(value) {self$vExpec + self$entropy},
    #' @field ICL double: value of the integrated classification log-likelihood
    ICL    = function(value) {-2 * self$vExpec + self$penalty}
  )
)

#' This internal class is designed to adjust a binary Stochastic Block Model in the context of missSBM.
#'
#' It is not designed not be call by the user
#'
#' @import R6
SimpleSBM_fit_noCov <-
R6::R6Class(classname = "SimpleSBM_fit_noCov",
  inherit = SimpleSBM_fit,
  public = list(
    #' @description update parameters estimation (M-step)
    #' @param ... additional arguments, only required for NMAR cases
    update_parameters = function(...) {
      res <- private$M_step(private$Y, private$R, private$Z, !self$directed)
      private$theta <- res$theta
      private$pi    <- as.numeric(res$pi)
      invisible(res)
    },
    #' @description update variational estimation of blocks (VE-step)
    #' @param ... additional arguments, only required for NMAR cases
    update_blocks =   function(...) {
      private$Z <- private$E_step(private$Y, private$R, private$Z, private$theta$mean, private$pi)
    }
  ),
  active = list(
    #' @field imputation the matrix of imputed values
    imputation = function(value) {
      as(.logistic(private$Z %*% log(private$theta$mean/(1-private$theta$mean)) %*% t(private$Z)) * private$S, "dgCMatrix")
    },
    #' @field vExpec double: variational approximation of the expectation complete log-likelihood
    vExpec = function(value) {
      private$vLL_complete(private$Y, private$R, private$Z, private$theta$mean, private$pi)
    }
  )
)

#' This internal class is designed to adjust a binary Stochastic Block Model in the context of missSBM.
#'
#' It is not designed not be call by the user
#'
#' @import R6
SimpleSBM_fit_withCov <-
R6::R6Class(classname = "SimpleSBM_fit_withCov",
  inherit = SimpleSBM_fit,
  public = list(
    #' @description update parameters estimation (M-step)
    #' @param control a list to tune nlopt for optimization, see documentation of nloptr
    #' @param ... use for compatibility
    update_parameters = function(...) {
      control <- list(maxeval = 50, xtol_rel = 1e-4, algorithm = "CCSAQ")
      res <- private$M_step(
        init_param = list(Gamma = .logit(private$theta$mean), beta = private$beta),
        Y = private$Y,
        R = private$R,
        X = self$covarArray,
        Z = private$Z,
        configuration = control
      )
      private$beta  <- as.numeric(res$beta)
      private$theta <- res$theta
      private$pi    <- as.numeric(res$pi)
      invisible(res)
    },
    #' @description update variational estimation of blocks (VE-step)
    #' @param ... use for compatibility
    update_blocks =   function(...) {
       private$Z <- private$E_step(private$Y, private$R, self$covarEffect, private$Z, .logit(private$theta$mean), private$pi, !self$directed, TRUE)
    }
  ),
  active = list(
    #' @field imputation the matrix of imputed values
    imputation = function(value) {
      as(.logistic(private$Z %*% log(private$theta$mean/(1-private$theta$mean)) %*% t(private$Z)) * private$S, "dgCMatrix")
    },
    #' @field vExpec double: variational approximation of the expectation complete log-likelihood
    vExpec = function(value) {
      private$vLL_complete(private$Y, private$R, self$covarEffect, private$Z, .logit(private$theta$mean), private$pi)
    }
  )
)

# vExpec = function(value) {
#   vLL_MAR <- private$vLL_complete(private$Y, private$R, private$Z, private$theta$mean, private$pi)
#   vLL_IMP <- private$vLL_complete(self$imputation, private$S, private$Z, private$theta$mean, private$pi)
#   vLL <- vLL_MAR + vLL_IMP - sum(private$Z %*% log(private$pi)) # counted twice
#   vLL
# }


#' This internal class is designed to adjust a binary Stochastic Block Model in the context of missSBM.
#'
#' It is not designed not be call by the user
#'
#' @import R6
SimpleSBM_fit_NMAR_noCov <-
R6::R6Class(classname = "SimpleSBM_NMAR_noCov",
  inherit = SimpleSBM_fit_noCov,
  private = list(
    V = NULL # matrix of imputed values
  ),
  public = list(
    #' @description constructor for simpleSBM_fit for missSBM purpose
    #' @param networkData a structure to store network under missing data condition: either a matrix possibly with NA, or a missSBM:::partlyObservedNetwork
    #' @param clusterInit Initial clustering: a vector with size \code{ncol(adjacencyMatrix)}, providing a user-defined clustering with \code{nbBlocks} levels.
    initialize = function(networkData, clusterInit) {
      super$initialize(networkData, clusterInit)
      private$V <- self$imputation
    },
    #' @description update parameters estimation (M-step)
    #' @param nu currently imputed values
    update_parameters = function(nu = NULL) {
      if (is.null(nu)) { # fall back to the MAR case
        super$update_parameters()
      } else {
        private$V <- nu
        ## only Bernoulli for the moment !!!
        Zbar <- colSums(private$Z)
        tZYZ <- t(private$Z) %*% private$Y %*% private$Z
        tZVZ <- t(private$Z) %*% private$V %*% private$Z
        if (self$directed) {
          private$theta$mean <- missSBM:::check_boundaries(as.matrix ( (tZYZ + tZVZ) / ( Zbar %o% Zbar - Zbar ) ))
        } else {
          private$theta$mean <-  missSBM:::check_boundaries(as.matrix ( (tZYZ + t(tZYZ) + tZVZ + t(tZVZ)) / ( Zbar %o% Zbar - Zbar ) ) )
        }
        private$pi <- colMeans(private$Z)
      }
    },
    #' @description update variational estimation of blocks (VE-step)
    #' @param log_lambda additional term sampling dependent used to de-bias estimation of tau
    update_blocks =   function(log_lambda = 0) {
      if (self$nbBlocks > 1) {
        log_tau_obs  <- private$E_step(private$Y, private$R, private$Z, private$theta$mean, private$pi, rescale = FALSE)
        log_tau_miss <- private$E_step(private$V, private$S, private$Z, private$theta$mean, private$pi, rescale = FALSE)
        private$Z    <- t(apply(log_tau_obs  + log_tau_miss + log_lambda, 1, .softmax))
      }
    }
  ),
  active = list(
    #' @field imputation the matrix of imputed values
    imputation = function(value) {
      as(.logistic(private$Z %*% log(private$theta$mean/(1-private$theta$mean)) %*% t(private$Z)) * private$S, "dgCMatrix")
    },
    #' @field vExpec double: variational approximation of the expectation complete log-likelihood
    vExpec = function(value) {
      vLL_MAR <- private$vLL_complete(private$Y, private$R, private$Z, private$theta$mean, private$pi)
      vLL_IMP <- private$vLL_complete(private$V, private$S, private$Z, private$theta$mean, private$pi)
      vLL <- vLL_MAR + vLL_IMP - sum(private$Z %*% log(private$pi)) # counted twice
      vLL
    }
  )
)

