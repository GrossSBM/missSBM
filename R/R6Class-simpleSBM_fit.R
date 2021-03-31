#' This internal class is designed to adjust a binary Stochastic Block Model in the context of missSBM.
#'
#' It is not designed not be call by the user
#'
#' @import R6
SimpleSBM_fit <-
R6::R6Class(classname = "SimpleSBM_fit",
  inherit = sbm::SimpleSBM,
  public = list(
    #' @description constructor for simpleSBM_fit for missSBM purpose
    #' @param adjacencyMatrix a matrix encoding the graph
    #' @param clusterInit Initial clustering: either a character in "hierarchical", "spectral" or "kmeans", or a vector with size \code{ncol(adjacencyMatrix)}, providing a user-defined clustering with \code{nbBlocks} levels. Default is "hierarchical".
    #' @param covarList An option list with M entries (the M covariates).
    initialize = function(adjacencyMatrix, clusterInit, covarList = list(), model = "bernoulli") {

      ## SANITY CHECKS (on data)
      stopifnot(is.matrix(adjacencyMatrix) | inherits(adjacencyMatrix, "Matrix")) # must be a matrix
      stopifnot(all.equal(nrow(adjacencyMatrix),
                          ncol(adjacencyMatrix)))             # matrix must be square
      stopifnot(all(sapply(covarList, nrow) == nrow(adjacencyMatrix))) # consistency of the covariates
      stopifnot(all(sapply(covarList, ncol) == ncol(adjacencyMatrix))) # with the network data

      # Basic fields initialization and call to super constructor
      super$initialize(model        = model,
                       directed     = ifelse(isSymmetric(adjacencyMatrix), FALSE, TRUE),
                       nbNodes      = nrow(adjacencyMatrix),
                       blockProp    = vector("numeric", 0),
                       connectParam = list(mean = matrix(0, 0, 0)),
                       covarList    = covarList)

      ## Initial Clustering
      private$Z <- clustering_indicator(clusterInit)
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
      objective[1] <- self$loglik
      iterate <- 1; stop <- FALSE

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
    #' @description permute group labels by order of decreasing probability
    reorder = function(){
      o <- order(private$theta$mean %*% private$pi, decreasing = TRUE)
      private$pi <- private$pi[o]
      private$theta$mean <- private$theta$mean[o,o]
      private$Z <- private$Z[, o, drop = FALSE]
    }
  ),
  active = list(
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
