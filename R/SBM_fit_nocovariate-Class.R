#' @import R6
#' @include SBM_fit-Class.R
SBM_fit_nocovariate <-
R6::R6Class(classname = "SBM_fit_nocovariate",
  inherit = SBM_fit,
  public = list(
    initialize = function(adjacencyMatrix, clusterInit) {

      # Basic fields intialization and call to super constructor
      nBlocks <- length(unique(clusterInit))
      super$initialize(
        directed     = ifelse(isSymmetric(adjacencyMatrix), FALSE, TRUE),
        nNodes       = nrow(adjacencyMatrix),
        mixtureParam = rep(NA, nBlocks),
        connectParam = matrix(NA, nBlocks, nBlocks)
      )
      private$Y <- adjacencyMatrix

      ## Initial Clustering
      private$tau <- clustering_indicator(clusterInit)

      ## Initialize estimation of the model parameters
      self$update_parameters()

      invisible(self)
    },
    update_parameters = function() { # NA not allowed in adjMatrix (should be imputed)
      private$pi    <- check_boundaries(quad_form(private$Y, private$tau) / quad_form(1 - diag(self$nNodes), private$tau))
      private$alpha <- check_boundaries(colMeans(private$tau))
    },
    update_blocks = function(log_lambda = 0) {
      if (private$Q > 1) {
        adjMatrix_bar <- bar(private$Y)
        ## Bernoulli undirected
        tau <- private$Y %*% private$tau %*% t(log(private$pi)) + adjMatrix_bar %*% private$tau %*% t(log(1 - private$pi)) + log_lambda
        if (private$directed) {
          ## Bernoulli directed
          tau <- tau + t(private$Y) %*% private$tau %*% t(log(t(private$pi))) + t(adjMatrix_bar) %*% private$tau %*% t(log(1 - t(private$pi)))
        }
        private$tau <- t(apply(sweep(tau, 2, log(private$alpha), "+"), 1, .softmax))
      }
    }
  ),
  active = list(
    vExpec = function(value) {
      prob   <- quad_form(private$pi, t(private$tau))
      factor <- ifelse(private$directed, 1, .5)
      adjMatrix_zeroDiag     <- private$Y ; diag(adjMatrix_zeroDiag) <- 0
      adjMatrix_zeroDiag_bar <- 1 - private$Y ; diag(adjMatrix_zeroDiag_bar) <- 0
      sum(private$tau %*% log(private$alpha)) +  factor * sum( adjMatrix_zeroDiag * log(prob) + adjMatrix_zeroDiag_bar *  log(1 - prob))
    }
  )
)

