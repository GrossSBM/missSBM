#' @import R6
#' @import nloptr
#' @include SBM_fit-Class.R
SBM_fit_covariates <-
R6::R6Class(classname = "SBM_fit_covariates",
  inherit = SBM_fit,
  public = list(
    initialize = function(adjacencyMatrix, clusterInit, covarList) {

      # Basic fields intialization and call to super constructor
      nbBlocks <- length(unique(clusterInit))
      super$initialize(
        directed     = ifelse(isSymmetric(adjacencyMatrix), FALSE, TRUE),
        nbNodes       = nrow(adjacencyMatrix),
        blockProp = rep(NA, nbBlocks),
        connectParam = matrix(NA, nbBlocks, nbBlocks),
        covarList   = covarList
      )
      private$Y <- adjacencyMatrix

      ## Initial Clustering
      Z <- clustering_indicator(clusterInit)

      ## Initialize parameters
      private$theta <- .logit(check_boundaries(quad_form(adjacencyMatrix, Z) / quad_form(1 - diag(self$nbNodes), Z)))
      private$alpha <- check_boundaries(colMeans(Z))
      private$beta  <- numeric(self$nbCovariates)
      private$tau   <- Z

      invisible(self)
    },
    update_parameters = function() { # NA not allowed in adjMatrix (should be imputed)
      optim_out  <-
        nloptr::nloptr(
          # starting parameters
          c(as.vector(private$theta),private$beta),
          # objective function + gradient
          ifelse(private$directed, Mstep_covariates_directed, Mstep_covariates_undirected),
          # optimizer parameters
          opts = list("algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-4),
          # additional argument for objective/gradient function
          Y = private$Y, cov = self$covarArray, Tau = private$tau,
        )
      private$beta  <- optim_out$solution[-(1:(self$nbBlocks^2))]
      private$theta <- matrix(optim_out$solution[1:(self$nbBlocks^2)], self$nbBlocks, self$nbBlocks)
      private$alpha <- check_boundaries(colMeans(private$tau))
    },
    update_blocks =   function(log_lambda = NULL) {
      private$tau <-
        E_step_covariates(
          private$Y,
          self$covarEffect,
          self$connectParam,
          private$tau,
          private$alpha
        )
    }
  ),
  active = list(
    vExpec = function(value) {
      vExpec_covariates(
        private$Y,
        self$covarEffect,
        self$connectParam,
        private$tau,
        private$alpha
      )
    }
  )
)
