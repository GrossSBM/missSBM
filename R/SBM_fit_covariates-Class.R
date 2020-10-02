#' @import R6
#' @import nloptr
#' @include SBM_fit-Class.R
SBM_fit_covariates <-
R6::R6Class(classname = "SBM_fit_covariates",
  inherit = SBM_fit,
  public = list(
    initialize = function(adjacencyMatrix, clusterInit, covarList) {

      # Basic fields intialization and call to super constructor
      super$initialize(
        adjacencyMatrix = adjacencyMatrix,
        model           = "bernoulli",
        directed        = ifelse(isSymmetric(adjacencyMatrix), FALSE, TRUE),
        covarList       = covarList
      )

      ## Initial Clustering
      Z <- clustering_indicator(clusterInit)

      ## Initialize parameters
      # private$theta <- list(mean = .logit(check_boundaries(quad_form(adjacencyMatrix, Z) / quad_form(1 - diag(self$nbNodes), Z))))
      private$theta <- list(mean = check_boundaries(quad_form(adjacencyMatrix, Z) / quad_form(1 - diag(self$nbNodes), Z)))
      private$pi    <- check_boundaries(colMeans(Z))
      private$beta  <- numeric(self$nbCovariates)
      private$tau   <- Z

      invisible(self)
    },
    update_parameters = function() { # NA not allowed in adjMatrix (should be imputed)
      optim_out  <-
        nloptr::nloptr(
          # starting parameters
          c(as.vector(.logit(private$theta$mean)),private$beta),
          # objective function + gradient
          ifelse(private$directed_, Mstep_covariates_directed, Mstep_covariates_undirected),
          # optimizer parameters
          opts = list("algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-4),
          # additional argument for objective/gradient function
          Y = private$Y, cov = self$covarArray, Tau = private$tau,
        )
      private$beta  <- optim_out$solution[-(1:(self$nbBlocks^2))]
      private$theta <- list(mean = matrix(.logistic(optim_out$solution[1:(self$nbBlocks^2)]), self$nbBlocks, self$nbBlocks))
      private$pi    <- check_boundaries(colMeans(private$tau))
    },
    update_blocks =   function(log_lambda = NULL) {
      private$tau <-
        E_step_covariates(
          private$Y,
          self$covarEffect,
          .logit(self$connectParam$mean),
          self$probMemberships,
          self$blockProp
        )
    }
  ),
  active = list(
    vExpec = function(value) {
      vExpec_covariates(
        private$Y,
        self$covarEffect,
        .logit(self$connectParam$mean),
        self$probMemberships,
        self$blockProp
      )
    }
  )
)
