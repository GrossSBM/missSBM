#' This internal class is designed to adjust a binary Stochastic Block Model in the context of missSBM.
#'
#' It is not designed not be call by the user
#'
#' @import R6
SimpleSBM_fit <-
R6::R6Class(classname = "SimpleSBM_fit",
  inherit = sbm::SimpleSBM,
  private = list(
    variant= NULL, # model variant
    R      = NULL, # the sampling matrix (sparse encoding)
    M_step = NULL, # pointing to the appropriate M_step function
    E_step = NULL, # pointing to the appropriate E_step function
    vLL_complete = NULL # pointing to the appropriate Expected complete LL function
  ),
  public = list(
    #' @description constructor for simpleSBM_fit for missSBM purpose
    #' @param adjacencyMatrix a matrix encoding the graph
    #' @param clusterInit Initial clustering: a vector with size \code{ncol(adjacencyMatrix)}, providing a user-defined clustering with \code{nbBlocks} levels.
    #' @param covarList An optional list with M entries (the M covariates).
    initialize = function(adjacencyMatrix, clusterInit, covarList = list()) {

### TODO: determine the model according to the adjacency matrix
### Choose in Bernoulli, Poisson, Gaussian, ZIGaussian
      model <- "bernoulli"

      ## SANITY CHECKS (on data)
      stopifnot(is.matrix(adjacencyMatrix) | inherits(adjacencyMatrix, "Matrix")) # must be a matrix or sparseMatrix
      stopifnot(all.equal(nrow(adjacencyMatrix), ncol(adjacencyMatrix)))          # matrix must be square
      ## Covariates are tested elsewhere

      ## Initial Clustering
      private$Z <- clustering_indicator(clusterInit)

      # Basic fields initialization and call to super constructor
      super$initialize(model        = model,
                       directed     = ifelse(isSymmetric(adjacencyMatrix), FALSE, TRUE),
                       nbNodes      = nrow(adjacencyMatrix),
                       blockProp    = vector("numeric", ncol(private$Z)) + .Machine$double.eps,
### FIXME: this should be model-dependent
                       connectParam = list(mean = matrix(0, ncol(private$Z), ncol(private$Z))),
                       covarList    = covarList)

      # Storing data
### TODO: handle the case where Y is a sparse Matrix

      ## where are my observations?
      diag(adjacencyMatrix) <- NA
      obs <- which(!is.na(adjacencyMatrix), arr.ind = TRUE)
      private$R <- Matrix::sparseMatrix(obs[,1], obs[,2],x = 1, dims = dim(adjacencyMatrix))

      ## where are my non-zero entries?
      nzero <- which(!is.na(adjacencyMatrix) & adjacencyMatrix != 0, arr.ind = TRUE)
      private$Y   <- Matrix::sparseMatrix(nzero[,1], nzero[,2], x = 1, dims = dim(adjacencyMatrix))

      ## point to the functions that performs E/M steps and compute the likelihood
      private$variant <-
        paste(model, ifelse(self$directed, "directed", "undirected"),
          ifelse(self$nbCovariates>0, "covariates", "nocovariate"), sep="_")
      private$M_step       <- get(paste("M_step_sparse"      , private$variant, sep = "_"))
      private$E_step       <- get(paste("E_step_sparse"      , private$variant, sep = "_"))
      private$vLL_complete <- get(paste("vLL_complete_sparse", private$variant, sep = "_"))

### TODO:
###  - check if parameters are not already intialize
###  - specialize the initialization to each model (this should be done in sbm::SImpleSBM...)

      ## Initialize estimation of the model parameters
      private$theta$mean <- matrix(0.5, ncol(private$Z), ncol(private$Z))
      private$beta       <- numeric(self$nbCovariates)
      self$update_parameters()

      invisible(self)
    },
    #' @description method to perform estimation via variational EM
    #' @param threshold stop when an optimization step changes the objective function by less than threshold. Default is 1e-4.
    #' @param maxIter V-EM algorithm stops when the number of iteration exceeds maxIter. Default is 10
    #' @param fixPointIter number of fix-point iterations in the Variational E step. Default is 3.
    #' @param trace logical for verbosity. Default is \code{FALSE}.
    doVEM = function(threshold = 1e-4, maxIter = 10, fixPointIter = 3, trace = FALSE) {

      ## Initialization of quantities that monitor convergence
      delta     <- vector("numeric", maxIter)
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
        delta[iterate] <- sqrt(sum((private$theta$mean - theta_old$mean)^2)) / sqrt(sum((theta_old$mean)^2))
        stop <- (iterate > maxIter) |  (delta[iterate] < threshold)
        objective[iterate] <- self$loglik
      }
      if (trace) cat("\n")
      res <- data.frame(delta = delta[1:iterate], objective = objective[1:iterate])
      res
    },
    #' @description permute group labels by order of decreasing probability
    reorder = function(){
      o <- order(private$theta$mean %*% private$pi, decreasing = TRUE)
      private$pi <- private$pi[o]
      private$theta$mean <- private$theta$mean[o,o]
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
    update_parameters = function() {
      res <- private$M_step(private$Y, private$R, private$Z)
      private$theta <- res$theta
      private$pi    <- as.numeric(res$pi)
      invisible(res)
    },
    #' @description update variational estimation of blocks (VE-step)
    update_blocks =   function() {
      private$Z <- private$E_step(private$Y, private$R, private$Z, private$theta$mean, private$pi)
    }
  ),
  active = list(
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
SimpleSBM_fit_NMAR_noCov <-
R6::R6Class(classname = "SimpleSBM_fit_noCov",
  inherit = SimpleSBM_fit,
  private = list(
    V = NULL, # matrix of imputed values  (sparse encoding)
    S = NULL  # the "anti" sampling matrix (sparse encoding)
  ),
  public = list(
    #' @description constructor for simpleSBM_fit for missSBM purpose
    #' @param adjacencyMatrix a matrix encoding the graph
    #' @param clusterInit Initial clustering: a vector with size \code{ncol(adjacencyMatrix)}, providing a user-defined clustering with \code{nbBlocks} levels.
    #' @param imputedValues a matrix encoding the imputed part of the network
    initialize = function(adjacencyMatrix, clusterInit) {
      super$initialize(adjacencyMatrix, clusterInit)
      diag(adjacencyMatrix) <- 0 ## not auto-loop: neither observed nor imputed
      miss <- which(is.na(adjacencyMatrix), arr.ind = TRUE)
      private$S <- Matrix::sparseMatrix(miss[,1], miss[,2], x = 1, dims = dim(adjacencyMatrix))
      self$update_parameters()
      private$V <- Matrix::sparseMatrix(miss[,1], miss[,2], x = 1, dims = dim(adjacencyMatrix))
    },
    #' @description update parameters estimation (M-step)
    update_parameters = function(imputed_values = NULL) {
      if (is.null(imputed_values)) {
        res <- private$M_step(private$Y, private$R, private$Z)
        private$theta <- res$theta
        private$pi    <- as.numeric(res$pi)
      } else {
        ## only Bernoulli for the moment !!!
        private$V@x <- imputed_values
        Zbar <- colSums(private$Z)
        private$theta$mean <- ( t(private$Z) %*% private$Y %*% private$Z + t(private$Z) %*% private$V %*% private$Z  ) / ( Zbar %o% Zbar - Zbar )
        private$pi         <- colMeans(private$Z)
      }
      invisible(res)
    },
    #' @description update variational estimation of blocks (VE-step)
    #' @param log_lambda additional term sampling dependent used to de-bias estimation of tau
    update_blocks =   function(log_lambda = 0) {
      log_tau_obs  <- private$E_step(private$Y, private$R   , private$Z, private$theta$mean, private$pi)
      log_tau_miss <- private$E_step(private$V, private$S, private$Z, private$theta$mean, private$pi)
      private$Z    <- t(apply(log_tau_obs  + log_tau_miss + log_lambda, 1, .softmax))
      private$Z    <- t(apply(log_tau_obs  + log_tau_miss + log_lambda, 1, .softmax))
    }
  ),
  active = list(
    #' @field vExpec double: variational approximation of the expectation complete log-likelihood
    vExpec = function(value) {
      vLL_MAR <- private$vLL_complete(private$Y, private$R, private$Z, private$theta$mean, private$pi)
      vLL_IMP <- private$vLL_complete(private$V, private$S, private$Z, private$theta$mean, private$pi)
      vLL <- vLL_MAR + vLL_IMP - sum(private$Z %*% log(private$pi)) # counted twice
      vLL
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
    update_parameters = function(control = list(maxeval = 50, xtol_rel = 1e-4, algorithm = "CCSAQ")) {
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
    #' @param log_lambda double use to adjust the parameter estimation according to the sampling design
    update_blocks =   function() {
       private$Z <- private$E_step(private$Y, private$R, roundProduct(private$X, private$beta), private$Z, .logit(private$theta$mean), private$pi)
    }
  ),
  active = list(
    #' @field vExpec double: variational approximation of the expectation complete log-likelihood
    vExpec = function(value) {
      private$vLL_complete(private$Y, private$R, roundProduct(private$X, private$beta), private$Z, .logit(private$theta$mean), private$pi)
    }
  )
)
