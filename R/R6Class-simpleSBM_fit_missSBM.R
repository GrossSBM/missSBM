#' This internal class is designed to adjust a binary Stochastic Block Model in the context of missSBM.
#'
#' It is not designed not be call by the user
#'
#' @import R6
#' @import nloptr
SimpleSBM_fit_missSBM <-
R6::R6Class(classname = "SimpleSBM_fit_missSBM",
  inherit = SimpleSBM_fit,
  public = list(
    #' @description constructor for simpleSBM_fit for missSBM purpose
    #' @param adjacencyMatrix a matrix encoding the graph
    #' @param clusterInit Initial clustering: either a character in "hierarchical", "spectral" or "kmeans", or a vector with size \code{ncol(adjacencyMatrix)}, providing a user-defined clustering with \code{nbBlocks} levels. Default is "hierarchical".
    #' @param covarList An option list with M entries (the M covariates).
    initialize = function(adjacencyMatrix, clusterInit, covarList = list()) {

      super$initialize(adjacencyMatrix, clusterInit, covarList, "bernoulli")

      # Storing data
      private$Y <- adjacencyMatrix

      ## Initialize estimation of the model parameters
      private$theta <- list(mean = check_boundaries(quad_form(adjacencyMatrix, private$Z) / quad_form(1 - diag(self$nbNodes), private$Z)))
      private$pi    <- check_boundaries(colMeans(private$Z))

      invisible(self)
    },
    #' @description update parameters estimation (M-step)
    update_parameters = function() { # NA not allowed in adjMatrix (should be imputed)
      if (self$nbCovariates > 0) {
        optim_out  <-
          nloptr::nloptr(
            # starting parameters
            c(as.vector(.logit(private$theta$mean)), private$beta),
            # objective function + gradient
            ifelse(private$directed_, Mstep_covariates_directed, Mstep_covariates_undirected),
            # optimizer parameters
            opts = list("algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-4),
            # additional argument for objective/gradient function
            Y = private$Y, cov = self$covarArray, Tau = private$Z,
          )
        private$beta  <- optim_out$solution[-(1:(self$nbBlocks^2))]
        private$theta <- list(mean = matrix(.logistic(optim_out$solution[1:(self$nbBlocks^2)]), self$nbBlocks, self$nbBlocks))
      } else {
        private$theta <- list(mean = check_boundaries(quad_form(private$Y, private$Z) / quad_form(1 - diag(self$nbNodes), private$Z)))
      }
      private$pi    <- check_boundaries(colMeans(private$Z))
    },
    #' @description update variational estimation of blocks (VE-step)
    #' @param log_lambda double use to adjust the parameter estimation according to the sampling design
    update_blocks =   function(log_lambda = 0) {
      if (self$nbCovariates > 0) {
        private$Z <-
          E_step_covariates(
            private$Y,
            self$covarEffect,
            matrix(.logit(self$connectParam$mean),self$nbBlocks, self$nbBlocks),
            self$probMemberships,
            self$blockProp
          )
      } else {
        if (self$nbBlocks > 1) {
          adjMatrix_bar <- bar(private$Y)
          ## Bernoulli undirected
          tau <- private$Y %*% private$Z %*% t(log(private$theta$mean)) + adjMatrix_bar %*% private$Z %*% t(log(1 - private$theta$mean)) + log_lambda
          if (private$directed_) {
            ## Bernoulli directed
            tau <- tau + t(private$Y) %*% private$Z %*% t(log(t(private$theta$mean))) + t(adjMatrix_bar) %*% private$Z %*% t(log(1 - t(private$theta$mean)))
          }
          private$Z <- t(apply(sweep(tau, 2, log(private$pi), "+"), 1, .softmax))
        }
      }
    }
  ),
  active = list(
    #' @field vExpec double: variational approximation of the expectation complete log-likelihood
    vExpec = function(value) {
      if (self$nbCovariates > 0) {
        res <- vExpec_covariates(
          private$Y,
          self$covarEffect,
          matrix(.logit(self$connectParam$mean),self$nbBlocks, self$nbBlocks),
          self$probMemberships,
          self$blockProp
        )
      } else {
        factor <- ifelse(private$directed_, 1, .5)
        adjMat <- private$Y ; diag(adjMat) <- 0
        tmp <- factor * sum( adjMat * private$Z %*% log(private$theta$mean) %*% t(private$Z) +
                               bar(private$Y)  *  private$Z %*% log(1 - private$theta$mean) %*% t(private$Z))
        res <- sum(private$Z %*% log(private$pi)) +  tmp
      }
      res
    }
  )
)
