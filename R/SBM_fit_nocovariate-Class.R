#' @import R6
#' @include SBM_fit-Class.R
SBM_fit_nocovariate <-
R6::R6Class(classname = "SBM_fit_nocovariate",
  inherit = SBM_fit,
  public = list(
    initialize = function(adjacencyMatrix, clusterInit) {

      # Basic fields intialization and call to super constructor
      nbBlocks <- length(unique(clusterInit))
      super$initialize(
        directed     = ifelse(isSymmetric(adjacencyMatrix), FALSE, TRUE),
        nbNodes       = nrow(adjacencyMatrix),
        blockProp = rep(NA, nbBlocks),
        connectParam = matrix(NA, nbBlocks, nbBlocks)
      )
      private$Y <- adjacencyMatrix

      ## Initial Clustering
      private$tau <- clustering_indicator(clusterInit)

      ## Initialize estimation of the model parameters
      self$update_parameters()

      invisible(self)
    },
    update_parameters = function() { # NA not allowed in adjMatrix (should be imputed)
      private$theta <- check_boundaries(quad_form(private$Y, private$tau) / quad_form(1 - diag(self$nbNodes), private$tau))
      private$pi <- check_boundaries(colMeans(private$tau))
    },
    update_blocks = function(log_lambda = 0) {
      if (self$nbBlocks > 1) {
        adjMatrix_bar <- bar(private$Y)
        ## Bernoulli undirected
        tau <- private$Y %*% private$tau %*% t(log(private$theta)) + adjMatrix_bar %*% private$tau %*% t(log(1 - private$theta)) + log_lambda
        if (private$directed_) {
          ## Bernoulli directed
          tau <- tau + t(private$Y) %*% private$tau %*% t(log(t(private$theta))) + t(adjMatrix_bar) %*% private$tau %*% t(log(1 - t(private$theta)))
        }
        # tau <- check_boundaries(t(apply(sweep(tau, 2, log(private$pi), "+"), 1, .softmax)), zero = 1e-4)
        # private$tau <- tau / matrix(rowSums(tau), self$nbNodes, self$nbBlocks, byrow = FALSE)
        private$tau <- t(apply(sweep(tau, 2, log(private$pi), "+"), 1, .softmax))
      }
    }
  ),
  active = list(
    vExpec = function(value) {
      factor <- ifelse(private$directed_, 1, .5)
      adjMat <- private$Y ; diag(adjMat) <- 0
      tmp <- factor * sum( adjMat * private$tau %*% log(private$theta) %*% t(private$tau) +
                             bar(private$Y)  *  private$tau %*% log(1 - private$theta) %*% t(private$tau))
      sum(private$tau %*% log(private$pi)) +  tmp
    }
  )
)

