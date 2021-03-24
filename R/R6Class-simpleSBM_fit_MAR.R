#' This internal class is designed to adjust a binary Stochastic Block Model in the context of missSBM.
#'
#' It is not designed not be call by the user
#'
#' @import R6
#' @import nloptr
SimpleSBM_fit_MAR <-
R6::R6Class(classname = "SimpleSBM_fit_MAR",
  inherit = sbm::SimpleSBM,
  private = list(
    R   = NULL, # the sampling matrix (sparse encoding)
    one = NULL, # indexes of ones in Y
    obs = NULL  # indexes of observed dyads in Y
  ),
  public = list(
    #' @description constructor for simpleSBM_fit for missSBM purpose
    #' @param adjacencyMatrix a matrix encoding the graph
    #' @param clusterInit Initial clustering: either a character in "hierarchical", "spectral" or "kmeans", or a vector with size \code{ncol(adjacencyMatrix)}, providing a user-defined clustering with \code{nbBlocks} levels. Default is "hierarchical".
    #' @param covarList An option list with M entries (the M covariates).
    initialize = function(adjacencyMatrix, clusterInit, covarList = list()) {

      ## SANITY CHECKS (on data)
      stopifnot(is.matrix(adjacencyMatrix))                   # must be a matrix
      stopifnot(all.equal(nrow(adjacencyMatrix),
                          ncol(adjacencyMatrix)))             # matrix must be square
      stopifnot(all(sapply(covarList, nrow) == nrow(adjacencyMatrix))) # consistency of the covariates
      stopifnot(all(sapply(covarList, ncol) == ncol(adjacencyMatrix))) # with the network data

      # Basic fields initialization and call to super constructor
      super$initialize(model        = "bernoulli",
                       directed     = ifelse(isSymmetric(adjacencyMatrix), FALSE, TRUE),
                       nbNodes      = nrow(adjacencyMatrix),
                       blockProp    = vector("numeric", 0),
                       connectParam = list(mean = matrix(0, 0, 0)),
                       covarList    = covarList)

      # Storing data
      ## where are my observations?
      diag(adjacencyMatrix) <- NA
      obs <- which(!is.na(adjacencyMatrix), arr.ind = TRUE)
      private$R <- Matrix::sparseMatrix(obs[,1], obs[,2],x = 1, dims = dim(adjacencyMatrix))

      ## where are my non-zero entries?
      one <- which(!is.na(adjacencyMatrix) & adjacencyMatrix != 0, arr.ind = TRUE)
      private$Y   <- Matrix::sparseMatrix(one[,1], one[,2], x = 1, dims = dim(adjacencyMatrix))

      ## Initial Clustering
      private$Z <- clustering_indicator(clusterInit)

      ## Initialize estimation of the model parameters
      self$update_parameters()

      invisible(self)
    },
    #' @description method to perform estimation via variational EM
    #' @param threshold stop when an optimization step changes the objective function by less than threshold. Default is 1e-4.
    #' @param maxIter V-EM algorithm stops when the number of iteration exceeds maxIter. Default is 10
    #' @param fixPointIter number of fix-point iterations in the Variational E step. Default is 3.
    #' @param trace logical for verbosity. Default is \code{FALSE}.
    doVEM = function(threshold = 1e-4, maxIter = 10, fixPointIter = 3, trace = FALSE) {

      ## Initialization of quantities that monitor convergence
      delta     <- vector("numeric", maxIter)
      objective <- vector("numeric", maxIter)
      iterate <- 0; stop <- FALSE

      ## Starting the variational EM algorithm
      if (trace) cat("\n Adjusting Variational EM for Stochastic Block Model\n")
      while (!stop) {
        iterate <- iterate + 1
        if (trace) cat(" iteration #:", iterate, "\r")

        theta_old <- private$theta # save old value of parameters to assess convergence

        # Variational E-Step
        for (i in seq.int(fixPointIter)) self$update_blocks()
        # M-step
        self$update_parameters()

        # Assess convergence
        delta[iterate] <- sqrt(sum((private$theta$mean - theta_old$mean)^2)) / sqrt(sum((theta_old$mean)^2))
        stop <- (iterate > maxIter) |  (delta[iterate] < threshold)
        objective[iterate] <- self$loglik
      }
      if (trace) cat("\n")
      res <- data.frame(delta = delta[1:iterate], objective = objective[1:iterate])
      res
    },
    #' @description update parameters estimation (M-step)
    update_parameters = function() {
        private$theta$mean <- as.matrix(t(private$Z) %*% private$Y %*% private$Z / t(private$Z) %*% private$R %*% private$Z)
        private$pi         <- check_boundaries(colMeans(private$Z))
    },
    #' @description update variational estimation of blocks (VE-step)
    #' @param log_lambda double use to adjust the parameter estimation according to the sampling design
    update_blocks =   function(log_lambda = 0) {
      if (self$nbBlocks > 1) {
        Z_log_piql_2 <- tcrossprod(private$Z, log(1 - private$theta$mean))
        Z_log_piql_1 <- tcrossprod(private$Z, log(private$theta$mean/(1 - private$theta$mean)))
        log_tau <- matrix(log_lambda + log(private$pi), self$nbNodes, self$nbBlocks, byrow = TRUE)
        ## Bernoulli undirected
        log_tau <- log_tau + private$Y %*% Z_log_piql_1 + private$R %*% Z_log_piql_2
        if (private$directed_) {
          ## Bernoulli directed
          ## log_tau <- log_tau + t(private$Y) %*% private$Z %*% t(log_piql_1) + t(private$R) %*% private$Z %*% t(log_piql_2)
          log_tau <- log_tau + Matrix::crossprod(private$Y, Z_log_piql_1) +  Matrix::crossprod(private$R, Z_log_piql_2)
        }
        private$Z <- t(apply(log_tau, 1, .softmax))
      }
    },
    #' @description permute group labels by order of decreasing probability
    reorder = function(){
      o <- order(private$theta$mean %*% private$pi, decreasing = TRUE)
      private$pi <- private$pi[o]
      private$theta$mean <- private$theta$mean[o,o]
      private$Z <- private$Z[, o, drop = FALSE]
    }
  ),
  active = list(
    #' @field vExpec double: variational approximation of the expectation complete log-likelihood
    vExpec = function(value) {
      log_piql_2 <- log(1 - private$theta$mean)
      log_piql_1 <- log(private$theta$mean) - log_piql_2
      res <- sum(t(private$Z) %*% private$Y %*% private$Z * log_piql_1) +
             sum(t(private$Z) %*% private$R %*% private$Z * log_piql_2)
      res <- ifelse(private$directed_, 1, .5) * res + sum(private$Z %*% log(private$pi))
      res
    },
    #' @field penalty double, value of the penalty term in ICL
    penalty  = function(value) {unname((self$nbConnectParam + self$nbCovariates) * log(self$nbDyads) + (self$nbBlocks-1) * log(self$nbNodes))},
    #' @field entropy double, value of the entropy due to the clustering distribution
    entropy  = function(value) {-sum(.xlogx(private$Z))},
    #' @field loglik double: approximation of the log-likelihood (variational lower bound) reached
    loglik = function(value) {self$vExpec + self$entropy},
    #' @field ICL double: value of the integrated classification log-likelihood
    ICL    = function(value) {-2 * self$vExpec + self$penalty}
  )
)
