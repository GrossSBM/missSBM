#' @import R6
#' @import nloptr
#' @include SBM_fit-Class.R
#' @export
SBM_fit_covariates <-
R6::R6Class(classname = "SBM_fit_covariates",
  inherit = SBM_fit,
  public = list(
    initialize = function(adjacencyMatrix, covariates, nBlocks, clusterInit = "spectral") {

      # Basic fields intialization and call to super constructor
      super$initialize(
        directed        = ifelse(isSymmetric(adjacencyMatrix), FALSE, TRUE),
        nNodes          = nrow(adjacencyMatrix),
        mixtureParam    = rep(NA,nBlocks),
        connectParam    = matrix(NA,nBlocks,nBlocks),
        covariates      = covariates,
        covarParam      = numeric(dim(covariates)[3])
      )

      ## Initial Clustering
      if (self$nBlocks > 1) {
        if (is.character(clusterInit)) {
          y <- as.vector(adjacencyMatrix)
          X <- apply(covariates, 3, as.vector)
          out_logistic <- glm.fit(X, y, family = binomial())
          adjacencyResiduals <- matrix(logistic(residuals(out_logistic)), self$nNodes, self$nNodes)
          clusterInit <-
            switch(clusterInit,
                   "hierarchical" = init_hierarchical(adjacencyResiduals, self$nBlocks),
                   "kmeans"       = init_kmeans(      adjacencyResiduals, self$nBlocks),
                                    init_spectral(    adjacencyResiduals, self$nBlocks)
            )
          Z <- matrix(0,self$nNodes,self$nBlocks)
          Z[cbind(1:self$nNodes, clusterInit)] <- 1
        } else if (is.numeric(clusterInit) | is.factor(clusterInit)) {
          Z <- matrix(0,self$nNodes,self$nBlocks)
          Z[cbind(1:self$nNodes, as.numeric(clusterInit))] <- 1
        } else {
          stop("unknown type for initial clustering")
        }
      } else {
        Z <- matrix(1, self$nNodes, self$nBlocks)
      }
      private$tau <- Z

      ## Initialize parameters
      self$init_parameters(adjacencyMatrix)

      invisible(self)
    },
    init_parameters = function(adjMatrix) { ## NA allowed in adjMatrix
      NAs           <- is.na(adjMatrix); adjMatrix[NAs] <- 0
      private$pi    <- logit(check_boundaries((t(private$tau) %*% (adjMatrix * !NAs) %*% private$tau) / (t(private$tau) %*% ((1 - diag(self$nNodes)) * !NAs) %*% private$tau)))
      private$alpha <- check_boundaries(colMeans(private$tau))
      private$beta  <- numeric(private$M)
    },
    update_parameters = function(adjMatrix) { # NA not allowed in adjMatrix (should be imputed)

      optim_out  <-
        nloptr::nloptr(
          # starting parameters
          c(as.vector(private$pi),private$beta),
          # objective function + gradient
          ifelse(private$directed, Mstep_covariates_directed, Mstep_covariates_undirected),
          # optimizer parameters
          opts = list("algorithm" = "NLOPT_LD_LBFGS", "xtol_rel" = 1.0e-4),
          # additional argument for objective/gradient function
          Y = adjMatrix, cov = private$phi, Tau = private$tau
        )
      private$beta  <- optim_out$solution[-(1:(private$Q^2))]
      private$pi    <- matrix(optim_out$solution[1:(private$Q^2)], private$Q, private$Q)
      private$alpha <- check_boundaries(colMeans(private$tau))

    },
    vExpec = function(adjMatrix) {
      vExpec_covariates(
        adjMatrix,
        roundProduct(private$phi, private$beta),
        private$pi,
        private$tau,
        private$alpha
      )
    },
    update_blocks =   function(adjMatrix, fixPointIter, log_lambda = NULL) {
      private$tau <-
        E_step_covariates(
          adjMatrix,
          roundProduct(private$phi, private$beta),
          private$pi,
          private$tau,
          private$alpha,
          fixPointIter
        )
    }
  )
)
