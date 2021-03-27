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
    R      = NULL, # the sampling matrix (sparse encoding)
    M_step = NULL, # pointing to the appropriate M_step function
    E_step = NULL, # pointing to the appropriate E_step function
    vbound = NULL  # pointing to the appropriate E_step function
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

      private$M_step <- M_step_sparse_bernoulli_undirected_nocovariate
      private$E_step <- E_step_sparse_bernoulli_undirected_nocovariate

      ## Initialize estimation of the model parameters
      self$update_parameters()

      invisible(self)
    },
    #' @description update parameters estimation (M-step)
    update_parameters = function() {
        private$theta$mean <- private$M_step(private$Y, private$R, private$Z)
        private$pi         <- check_boundaries(colMeans(private$Z))
    },
    #' @description update variational estimation of blocks (VE-step)
    #' @param log_lambda double use to adjust the parameter estimation according to the sampling design
    update_blocks =   function(log_lambda = 0) {
      if (self$nbBlocks > 1) {
        private$Z <- private$E_step(private$Y, private$R, private$Z, private$theta$mean, private$pi, log_lambda)
      }
    }
  ),
  active = list(
    #' @field vExpec double: variational approximation of the expectation complete log-likelihood
    vExpec = function(value) {
      log_piql_2 <- log(1 - private$theta$mean)
      log_piql_1 <- log(private$theta$mean) - log_piql_2
      res <- sum(t(private$Z) %*% private$Y %*% private$Z * log_piql_1) +
             sum( t(private$Z) %*% private$R %*% private$Z * log_piql_2)
      res <- ifelse(private$directed_, 1, .5) * res + sum(private$Z %*% log(private$pi))
      res
    }
  )
)
