#' @import R6
#' @import nloptr
#' @include SBM_fit-Class.R
SBM_fit_covariates <-
R6::R6Class(classname = "SBM_fit_covariates",
  inherit = SBM_fit,
  public = list(
    initialize = function(adjacencyMatrix, clusterInit, covarArray) {

      # Basic fields intialization and call to super constructor
      nbBlocks <- length(unique(clusterInit))
      super$initialize(
        directed     = ifelse(isSymmetric(adjacencyMatrix), FALSE, TRUE),
        nbNodes       = nrow(adjacencyMatrix),
        blockProp = rep(NA, nbBlocks),
        connectParam = matrix(NA, nbBlocks, nbBlocks),
        covarParam   = numeric(dim(covarArray)[3]),
        covarArray   = covarArray
      )
      private$Y <- adjacencyMatrix

      ## Initial Clustering
      Z <- clustering_indicator(clusterInit)

      ## Initialize parameters
      private$pi    <- logit(check_boundaries(quad_form(adjacencyMatrix, Z) / quad_form(1 - diag(self$nbNodes), Z)))
      private$alpha <- check_boundaries(colMeans(Z))
      private$beta  <- numeric(private$M)
      private$tau   <- Z

      invisible(self)
    },
    update_parameters = function() { # NA not allowed in adjMatrix (should be imputed)
      optim_out  <-
        nloptr::nloptr(
          # starting parameters
          c(as.vector(private$pi),private$beta),
          # objective function + gradient
          ifelse(private$directed, Mstep_covariates_directed, Mstep_covariates_undirected),
          # optimizer parameters
          opts = list("algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-4),
          # additional argument for objective/gradient function
          Y = private$Y, cov = private$X, Tau = private$tau,
        )
      private$beta  <- optim_out$solution[-(1:(private$Q^2))]
      private$pi    <- matrix(optim_out$solution[1:(private$Q^2)], private$Q, private$Q)
      private$alpha <- check_boundaries(colMeans(private$tau))
    },
    update_blocks =   function(log_lambda = NULL) {
      private$tau <-
        E_step_covariates(
          private$Y,
          roundProduct(private$X, private$beta),
          private$pi,
          private$tau,
          private$alpha
        )
    }
  ),
  active = list(
    vExpec = function(value) {
      vExpec_covariates(
        private$Y,
        roundProduct(private$X, private$beta),
        private$pi,
        private$tau,
        private$alpha
      )
    }
  )
)
