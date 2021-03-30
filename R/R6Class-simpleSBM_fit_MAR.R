#' This internal class is designed to adjust a binary Stochastic Block Model in the context of missSBM.
#'
#' It is not designed not be call by the user
#'
#' @import R6
#' @import nloptr
#' @import Matrix
SimpleSBM_fit_MAR <-
R6::R6Class(classname = "SimpleSBM_fit_MAR",
  inherit = SimpleSBM_fit,
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
    #' @param clusterInit Initial clustering: either a character in "hierarchical", "spectral" or "kmeans", or a vector with size \code{ncol(adjacencyMatrix)}, providing a user-defined clustering with \code{nbBlocks} levels. Default is "hierarchical".
    #' @param covarList An option list with M entries (the M covariates).
    initialize = function(adjacencyMatrix, clusterInit, covarList = list(), model = "bernoulli") {

      super$initialize(adjacencyMatrix, clusterInit, covarList, model)

      # Storing data
      ## where are my observations?
      diag(adjacencyMatrix) <- NA
      obs <- which(!is.na(adjacencyMatrix), arr.ind = TRUE)
      private$R <- Matrix::sparseMatrix(obs[,1], obs[,2],x = 1, dims = dim(adjacencyMatrix))

      ## where are my non-zero entries?
      nzero <- which(!is.na(adjacencyMatrix) & adjacencyMatrix != 0, arr.ind = TRUE)
      private$Y   <- Matrix::sparseMatrix(nzero[,1], nzero[,2], x = 1, dims = dim(adjacencyMatrix))

      private$variant <-
        paste(model,
          ifelse(self$directed, "directed", "undirected"),
          ifelse(self$nbCovariates>0, "covariates", "nocovariate"), sep="_")
      private$M_step       <- get(paste("M_step_sparse"      , private$variant, sep = "_"))
      private$E_step       <- get(paste("E_step_sparse"      , private$variant, sep = "_"))
      private$vLL_complete <- get(paste("vLL_complete_sparse", private$variant, sep = "_"))

      ## Initialize estimation of the model parameters
      private$theta$mean <- matrix(0.5, ncol(private$Z), ncol(private$Z))
      private$beta       <- numeric(self$nbCovariates)
      self$update_parameters()

      invisible(self)
    },
    #' @description update parameters estimation (M-step)
    #' @param control a list to tune nlopt for optimization, see documentation of nloptr. Only relevant when covariates are present
    update_parameters = function(control = list(maxeval = 50, xtol_rel = 1e-4, algorithm = "CCSAQ")) {
      if (self$nbCovariates > 0) {
        res <- private$M_step(
          init_param = list(Gamma = .logit(private$theta$mean), beta = private$beta),
                   Y = private$Y,
                   R = private$R,
                   X = self$covarArray,
                   Z = private$Z,
          configuration = control
         )
         private$beta  <- as.numeric(res$beta)
      } else {
        res <- private$M_step(private$Y, private$R, private$Z)
      }
      private$theta <- res$theta
      private$pi    <- as.numeric(res$pi)
      invisible(res)
    },
    #' @description update variational estimation of blocks (VE-step)
    #' @param log_lambda double use to adjust the parameter estimation according to the sampling design
    update_blocks =   function(log_lambda = 0) {
      if (self$nbBlocks > 1) {
        if (self$nbCovariates > 0) {
          private$Z <- private$E_step(private$Y, private$R, roundProduct(private$X, private$beta), private$Z, .logit(private$theta$mean), private$pi, log_lambda)
        } else {
          private$Z <- private$E_step(private$Y, private$R, private$Z, private$theta$mean, private$pi, log_lambda)
        }
      }
    }
  ),
  active = list(
    #' @field vExpec double: variational approximation of the expectation complete log-likelihood
    vExpec = function(value) {
      if (self$nbCovariates > 0) {
        res <- private$vLL_complete(private$Y, private$R, roundProduct(private$X, private$beta), private$Z, .logit(private$theta$mean), private$pi)
      } else {
        res <- private$vLL_complete(private$Y, private$R, private$Z, private$theta$mean, private$pi)
      }
      res
    }
  )
)
