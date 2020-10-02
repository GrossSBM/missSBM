#' R6 Class definition of an SBM-fit
#'
#' This class is designed to adjust a Stochastic Block Model on a fully observed network.
#' The doVEM method performs inference via Variational EM.
#'
#' This class is virtual: inference is effective only for instance of one of the two child
#' classes 'SBM_fit_nocovariate' and 'SBM_fit_covariates'
#'
#' @inherit SBM details
#'
#' @import R6
SBM_fit <-
R6::R6Class(classname = "SBM_fit",
  inherit = SBM,
  private = list(
    tau  = NULL  # variational parameters for posterior probablility of class belonging
  ),
  public = list(
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
        delta[iterate] <- sqrt(sum((private$theta - theta_old)^2)) / sqrt(sum((theta_old)^2))
        stop <- (iterate > maxIter) |  (delta[iterate] < threshold)
        objective[iterate] <- self$loglik
      }
      if (trace) cat("\n")
      res <- data.frame(delta = delta[1:iterate], objective = objective[1:iterate])
      res
    },
    #' @description show method
    #' @param type character used to specify the type of SBM
    show = function(type = "Fit of a Stochastic Block Model\n") {
      super$show(type)
      cat("* Additional fields \n")
      cat("  $probMemberships, $memberships, $adjMatrix, $expectation, $loglik, $ICL, $penalty\n")
      cat("* S3 methods \n")
      cat("  plot, print, summary, coef, fitted \n")
    }
  ),
  active = list(
    #' @field loglik double: approximation of the log-likelihood (variational lower bound) reached
    loglik      = function(value) {self$vExpec + self$entropy},
    #' @field ICL double: value of the integrated classification log-likelihood
    ICL        = function(value) {-2 * self$vExpec + self$penalty},
    #' @field probMemberships matrix for clustering memberships
    probMemberships      = function(value) {if (missing(value)) return(private$tau) else  private$tau <- value},
    #' @field memberships vector of clustering
    memberships = function(value) {apply(private$tau, 1, which.max)},
    #' @field penalty double, value of the penalty term in ICL
    penalty     = function(value) {(self$nbConnectParam + self$nbCovariates) * log(self$nbDyads) + (self$nbBlocks-1) * log(self$nbNodes)},
    #' @field entropy double, value of the entropy due to the clustering distribution
    entropy     = function(value) {-sum(xlogx(private$tau))},
    #' @field expectation expected values of connection under the current model
    expectation = function(value) {
      if (self$nbCovariates > 0) {## theta is gamma in covariate SBM
        Prob <- check_boundaries(.logistic(private$tau %*% private$theta %*% t(private$tau) + self$covarEffect))
      } else {
        Prob <- check_boundaries(.logistic(private$tau %*% .logit(private$theta) %*% t(private$tau)))
      }
      Prob
    }
  )
)

