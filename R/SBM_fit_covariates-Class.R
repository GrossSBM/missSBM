#' @import R6
#' @include SBM_fit-Class.R
#' @export
SBM_fit_covariates <-
R6::R6Class(classname = "SBM_fit_covariates",
  inherit = SBM_fit,
  public = list(
    initialize = function(adjacencyMatrix, covariates, nBlocks, clusterInit = "spectral", covarSimilarity = sim_abs) {

      # Basic fields intialization and call to super constructor
      super$initialize(
        directed        = ifelse(isSymmetric(adjacencyMatrix), FALSE, TRUE),
        nNodes          = nrow(adjacencyMatrix),
        covariates      = covariates,
        mixtureParam    = rep(NA,nBlocks),
        connectParam    = matrix(NA,nBlocks,nBlocks),
        covarParam      = numeric(ncol(covariates)),
        covarSimilarity = covarSimilarity
      )

      ## Initial Clustering
      ## TODO: have an initialization specific to covariates
      if (self$nBlocks > 1) {
        if (is.character(clusterInit)) {
          clusterInit <-
            switch(clusterInit,
                   "hierarchical" = init_hierarchical(adjacencyMatrix, self$nBlocks),
                   "kmeans"       = init_kmeans(      adjacencyMatrix, self$nBlocks),
                   init_spectral(    adjacencyMatrix, self$nBlocks)
            )
          Z <- matrix(0,self$nNodes,self$nBlocks)
          Z[cbind(1:self$nNodes, clusterInit)] <- 1
        } else if (is.numeric(clusterInit)) {
          Z <- matrix(0,self$nNodes,self$nBlocks)
          Z[cbind(1:self$nNodes, clusterInit)] <- 1
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
      pi_           <- check_boundaries((t(private$tau) %*% (adjMatrix * !NAs) %*% private$tau) / (t(private$tau) %*% ((1 - diag(self$nNodes)) * !NAs) %*% private$tau))
      private$pi    <- log(pi_/(1 - pi_))
      private$alpha <- check_boundaries(colMeans(private$tau))
      private$beta  <- numeric(private$M)
    },
    update_parameters = function(adjMatrix) { # NA not allowed in adjMatrix (should be imputed)
        param <- c(as.vector(private$pi),private$beta)
        optim_out <- optim(param, objective_Mstep_covariates, gradient_Mstep_covariates,
          Y = adjMatrix, cov = private$cov, Tau = private$tau, directed = private$directed,
          method = "BFGS", control = list(fnscale = -1)
        )
        private$beta  <- optim_out$par[-(1:(Q^2))]
        private$pi    <- matrix(optim_out$par[1:(private$Q^2)], private$Q, private$Q)
        private$alpha <- check_boundaries(colMeans(private$tau))
    },
    vExpec = function(adjMatrix) {
      J <- vBound_covariates(adjMatrix, private$pi, private$beta, private$cov, private$tau, priavte$alpha)
      J
    },
    update_blocks =   function(adjMatrix, fixPointIter, log_lambda = 0) {
      ## TODO: check how log_lambda hsould be handle...
      ## TODO: include the loop in C++
      for (i in 1:fixPointIter) {
        private$tau <- E_step_covariates(adjMatrix, private$cov, private$pi, private$beta, private$tau, private$alpha)
      }
    }
  )
)
