#' This internal class is designed to adjust a binary Stochastic Block Model in the context of missSBM.
#'
#' It is not designed not be call by the user
#'
#' @import R6
#' @import nloptr
SimpleSBM_fit_missSBM <-
R6::R6Class(classname = "SimpleSBM_fit_missSBM",
  inherit = sbm::SimpleSBM,
  public = list(
    #' @description constructor for simpleSBM_fit for missSBM purpose
    #' @param adjacencyMatrix a matrix encoding the graph
    #' @param clusterInit Initial clustering: either a character in "hierarchical", "spectral" or "kmeans", or a vector with size \code{ncol(adjacencyMatrix)}, providing a user-defined clustering with \code{nbBlocks} levels. Default is "hierarchical".
    #' @param covarList An option list with M entries (the M covariates).
    initialize = function(adjacencyMatrix, clusterInit, covarList = list()) {

      ## SANITY CHECKS (on data)
      stopifnot(is.matrix(adjacencyMatrix) | inherits(adjacencyMatrix, "Matrix")) # must be a matrix
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

      ## Initial Clustering
      private$Z <- clustering_indicator(clusterInit)

      # Storing data
      private$Y <- adjacencyMatrix

      ## Initialize estimation of the model parameters
      private$theta <- list(mean = check_boundaries(quad_form(adjacencyMatrix, private$Z) / quad_form(1 - diag(self$nbNodes), private$Z)))
      private$pi    <- check_boundaries(colMeans(private$Z))

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
    },
    #' @description update parameters estimation (M-step)
    update_parameters = function() { # NA not allowed in adjMatrix (should be imputed)
      if (self$nbCovariates > 0) {
        optim_out  <-
          nloptr::nloptr(
            # starting parameters
            c(as.vector(.logit(private$theta$mean)), private$beta),
            # objective function + gradient
            ifelse(private$directed_, Mstep_covariates_directed, Mstep_covariates_undirected),
            # optimizer parameters
            opts = list("algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-4),
            # additional argument for objective/gradient function
            Y = private$Y, cov = self$covarArray, Tau = private$Z,
          )
        private$beta  <- optim_out$solution[-(1:(self$nbBlocks^2))]
        private$theta <- list(mean = matrix(.logistic(optim_out$solution[1:(self$nbBlocks^2)]), self$nbBlocks, self$nbBlocks))
      } else {
        private$theta <- list(mean = check_boundaries(quad_form(private$Y, private$Z) / quad_form(1 - diag(self$nbNodes), private$Z)))
      }
      private$pi    <- check_boundaries(colMeans(private$Z))
    },
    #' @description update variational estimation of blocks (VE-step)
    #' @param log_lambda double use to adjust the parameter estimation according to the sampling design
    update_blocks =   function(log_lambda = 0) {
      if (self$nbCovariates > 0) {
        private$Z <-
          E_step_covariates(
            private$Y,
            self$covarEffect,
            matrix(.logit(self$connectParam$mean),self$nbBlocks, self$nbBlocks),
            self$probMemberships,
            self$blockProp
          )
      } else {
        if (self$nbBlocks > 1) {
          adjMatrix_bar <- bar(private$Y)
          ## Bernoulli undirected
          tau <- private$Y %*% private$Z %*% t(log(private$theta$mean)) + adjMatrix_bar %*% private$Z %*% t(log(1 - private$theta$mean)) + log_lambda
          if (private$directed_) {
            ## Bernoulli directed
            tau <- tau + t(private$Y) %*% private$Z %*% t(log(t(private$theta$mean))) + t(adjMatrix_bar) %*% private$Z %*% t(log(1 - t(private$theta$mean)))
          }
          private$Z <- t(apply(sweep(tau, 2, log(private$pi), "+"), 1, .softmax))
        }
      }
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
    ICL    = function(value) {-2 * self$vExpec + self$penalty},
    #' @field vExpec double: variational approximation of the expectation complete log-likelihood
    vExpec = function(value) {
      if (self$nbCovariates > 0) {
        res <- vExpec_covariates(
          private$Y,
          self$covarEffect,
          matrix(.logit(self$connectParam$mean),self$nbBlocks, self$nbBlocks),
          self$probMemberships,
          self$blockProp
        )
      } else {
        factor <- ifelse(private$directed_, 1, .5)
        adjMat <- private$Y ; diag(adjMat) <- 0
        tmp <- factor * sum( adjMat * private$Z %*% log(private$theta$mean) %*% t(private$Z) +
                               bar(private$Y)  *  private$Z %*% log(1 - private$theta$mean) %*% t(private$Z))
        res <- sum(private$Z %*% log(private$pi)) +  tmp
      }
      res
    }
  )
)
