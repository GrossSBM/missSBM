#' @import R6
#' @include SBM_fit-Class.R
#' @export
SBM_fit_nocovariate <-
R6::R6Class(classname = "SBM_fit_nocovariate",
  inherit = SBM_fit,
  public = list(
    initialize =   function(adjacencyMatrix, clusterInit) {

      stopifnot(all.equal(nrow(adjacencyMatrix), length(clusterInit)))
      nBlocks <- length(unique(clusterInit))
      super$initialize( # call to super constructor
        directed     = ifelse(isSymmetric(adjacencyMatrix), FALSE, TRUE),
        nNodes       = nrow(adjacencyMatrix),
        mixtureParam = rep(NA,nBlocks),
        connectParam = matrix(NA,nBlocks,nBlocks)
      )

      ## Initial Clustering
      private$tau <- clustering_indicator(clusterInit)

      ## Initialize estimation of the model parameters
      self$update_parameters(adjacencyMatrix)

      invisible(self)
    },
    update_parameters = function(adjMatrix) { # NA not allowed in adjMatrix (should be imputed)
      private$pi    <- check_boundaries((t(private$tau) %*% adjMatrix %*% private$tau) / (t(private$tau) %*% (1 - diag(self$nNodes)) %*% private$tau))
      private$alpha <- check_boundaries(colMeans(private$tau))
    },
    vExpec = function(adjMatrix) {
      prob   <- private$tau %*% private$pi %*% t(private$tau)
      factor <- ifelse(private$directed, 1, .5)
      adjMatrix_zeroDiag     <- adjMatrix ; diag(adjMatrix_zeroDiag) <- 0
      adjMatrix_zeroDiag_bar <- 1 - adjMatrix ; diag(adjMatrix_zeroDiag_bar) <- 0
      sum(private$tau %*% log(private$alpha)) +  factor * sum( adjMatrix_zeroDiag * log(prob) + adjMatrix_zeroDiag_bar *  log(1 - prob))
    },
    update_blocks = function(adjMatrix, log_lambda = 0) {
      if (private$Q > 1) {
        adjMatrix_bar <- bar(adjMatrix)
        ## Bernoulli undirected
        tau <- adjMatrix %*% private$tau %*% t(log(private$pi)) + adjMatrix_bar %*% private$tau %*% t(log(1 - private$pi)) + log_lambda
        if (private$directed) {
          ## Bernoulli directed
          tau <- tau + t(adjMatrix) %*% private$tau %*% t(log(t(private$pi))) + t(adjMatrix_bar) %*% private$tau %*% t(log(1 - t(private$pi)))
        }
        private$tau <- t(apply(sweep(tau, 2, log(private$alpha), "+"), 1, .softmax))
      }
    }
  )
)

