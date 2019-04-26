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
      factor <- ifelse(private$directed, 1, .5)
      adjMat <- private$Y ; diag(adjMat) <- 0
      tmp <- factor * sum( adjMat * private$tau %*% log(private$pi) %*% t(private$tau) +
                             bar(private$Y)  *  private$tau %*% log(1 - private$pi) %*% t(private$tau))
      sum(private$tau %*% log(private$alpha)) +  tmp
    }
  )
)

