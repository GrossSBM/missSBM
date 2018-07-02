#' @import R6
#' @include SBM_fit-Class.R
#' @export
SBM_fit_covariates <-
R6::R6Class(classname = "SBM_fit_covariates",
  inherit = SBM_fit,
  public = list(
    init_parameters = function(adjMatrix) { ## NA allowed in adjMatrix
      NAs           <- is.na(adjMatrix); adjMatrix[NAs] <- 0
      pi_           <- check_boundaries((t(private$tau) %*% (adjMatrix * !NAs) %*% private$tau) / (t(private$tau) %*% ((1 - diag(self$nNodes)) * !NAs) %*% private$tau))
      private$pi    <- log(pi_/(1 - pi_))
      private$alpha <- check_boundaries(colMeans(private$tau))
      private$beta  <- numeric(private$M)
    },
    update_parameters = function(adjMatrix) { # NA not allowed in adjMatrix (should be imputed)
        param <- c(as.vector(private$pi),private$beta)
        optim_out <- optim(param, objective_Mstep_covariates, gradient_Mstep_covariates,
          Y = adjMatrix, cov = private$X, Tau = private$tau, directed = private$directed,
          method = "BFGS", control = list(fnscale = -1)
        )
        private$beta  <- optim_out$par[-(1:(Q^2))]
        private$pi    <- matrix(optim_out$par[1:(private$Q^2)], private$Q, private$Q)
        private$alpha <- check_boundaries(colMeans(private$tau))
    },
    vExpec = function(adjMatrix) {
      # prob   <- private$tau %*% private$pi %*% t(private$tau)
      # factor <- ifelse(private$directed, 1, .5)
      # adjMatrix_zeroDiag     <- adjMatrix ; diag(adjMatrix_zeroDiag) <- 0           ### Changement ici ###
      # adjMatrix_zeroDiag_bar <- 1 - adjMatrix ; diag(adjMatrix_zeroDiag_bar) <- 0   ### Changement ici ###
      # sum(private$tau %*% log(private$alpha)) +  factor * sum( adjMatrix_zeroDiag * log(prob) + adjMatrix_zeroDiag_bar *  log(1 - prob))
    }
  )
)

### !!!TODO!!! Write this in C++ and update each individual coordinate-wise
SBM_fit_covariates$set("public", "update_blocks",
  function(adjMatrix, fixPointIter, log_lambda = 0) {

    private$tau <- E_step(adjMatrix, private$X, private$pi, private$beta, private$tau, private$alpha)

    # adjMatrix_bar <- bar(adjMatrix)
    # for (i in 1:fixPointIter) {
    #   ## Bernoulli undirected
    #   tau <- adjMatrix %*% private$tau %*% t(log(private$pi)) + adjMatrix_bar %*% private$tau %*% t(log(1 - private$pi)) + log_lambda
    #   if (private$directed) {
    #     ## Bernoulli directed
    #     tau <- tau + t(adjMatrix) %*% private$tau %*% t(log(t(private$pi))) + t(adjMatrix_bar) %*% private$tau %*% t(log(1 - t(private$pi)))
    #   }
    #   tau <- exp(sweep(tau, 2, log(private$alpha),"+"))
    #   tau <- tau/rowSums(tau)
    # }
    # private$tau <- tau
  }
)
