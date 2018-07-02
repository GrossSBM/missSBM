#' @import R6
#' @include SBM_fit-Class.R
#' @export
SBM_fit_covariates <-
R6::R6Class(classname = "SBM_fit_covariates",
  inherit = SBM_fit,
  public = list(
    init_parameters = function(adjMatrix) { ## NA allowed in adjMatrix
      # NAs <- is.na(adjMatrix); adjMatrix[NAs] <- 0
      # private$pi    <- check_boundaries((t(private$tau) %*% (adjMatrix * !NAs) %*% private$tau) / (t(private$tau) %*% ((1 - diag(self$nNodes)) * !NAs) %*% private$tau))
      # private$alpha <- check_boundaries(colMeans(private$tau))
    },
    update_parameters = function(adjMatrix) { # NA not allowed in adjMatrix (should be imputed)
      # private$pi    <- check_boundaries((t(private$tau) %*% adjMatrix %*% private$tau) / (t(private$tau) %*% (1 - diag(self$nNodes)) %*% private$tau))
      # private$alpha <- check_boundaries(colMeans(private$tau))
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


