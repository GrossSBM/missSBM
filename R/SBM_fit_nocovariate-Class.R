#' @import R6
#' @include SBM_fit-Class.R
#' @export
SBM_fit_nocovariate <-
R6::R6Class(classname = "SBM_fit_nocovariate",
  inherit = SBM_fit,
  public = list(
    initialize =   function(adjacencyMatrix, nBlocks, clusterInit = "spectral") {

      super$initialize( # call to super constructor
        directed     = ifelse(isSymmetric(adjacencyMatrix), FALSE, TRUE),
        nNodes       = nrow(adjacencyMatrix),
        mixtureParam = rep(NA,nBlocks),
        connectParam = matrix(NA,nBlocks,nBlocks)
      )

      ## Initial Clustering
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

      ## Initialize estimation of the model parameters
      self$init_parameters(adjacencyMatrix)

      invisible(self)
    },
    init_parameters = function(adjMatrix) { ## NA allowed in adjMatrix
      NAs <- is.na(adjMatrix); adjMatrix[NAs] <- 0
      private$pi    <- check_boundaries((t(private$tau) %*% (adjMatrix * !NAs) %*% private$tau) / (t(private$tau) %*% ((1 - diag(self$nNodes)) * !NAs) %*% private$tau))
      private$alpha <- check_boundaries(colMeans(private$tau))
    },
    update_parameters = function(adjMatrix) { # NA not allowed in adjMatrix (should be imputed)
      private$pi    <- check_boundaries((t(private$tau) %*% adjMatrix %*% private$tau) / (t(private$tau) %*% (1 - diag(self$nNodes)) %*% private$tau))
      private$alpha <- check_boundaries(colMeans(private$tau))
    },
    vExpec = function(adjMatrix) {
      prob   <- private$tau %*% private$pi %*% t(private$tau)
      factor <- ifelse(private$directed, 1, .5)
      adjMatrix_zeroDiag     <- adjMatrix ; diag(adjMatrix_zeroDiag) <- 0           ### Changement ici ###
      adjMatrix_zeroDiag_bar <- 1 - adjMatrix ; diag(adjMatrix_zeroDiag_bar) <- 0   ### Changement ici ###
      sum(private$tau %*% log(private$alpha)) +  factor * sum( adjMatrix_zeroDiag * log(prob) + adjMatrix_zeroDiag_bar *  log(1 - prob))
    },
    update_blocks = function(adjMatrix, fixPointIter, log_lambda = 0) {
      ## TODO: code this function in Rcpp/C++ ...
      adjMatrix_bar <- bar(adjMatrix)
      tau_old <- private$tau
      for (i in 1:fixPointIter) {
        ## Bernoulli undirected
        tau <- adjMatrix %*% tau_old %*% t(log(private$pi)) + adjMatrix_bar %*% tau_old %*% t(log(1 - private$pi)) + log_lambda
        if (private$directed) {
          ## Bernoulli directed
          tau <- tau + t(adjMatrix) %*% tau_old %*% t(log(t(private$pi))) + t(adjMatrix_bar) %*% tau_old %*% t(log(1 - t(private$pi)))
        }
        tau <- exp(sweep(tau, 2, log(private$alpha),"+"))
        tau <- tau/rowSums(tau)
        tau_old <- tau
      }
      private$tau <- tau
    }
  )

)

