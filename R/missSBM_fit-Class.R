#' An R6 class to represent an SBM fit with missing data
#'
#' @description The function [estimate()] fits a collection of SBM for varying number of block.
#' Each fitted SBM is an instance of an R6 object with class ['missSBM_fit'], described here.
#'
#' Fields are accessed via active binding and cannot be changed by the user.
#'
#' This class comes with a set of R6 methods, some of them being useful for the user and exported
#' as S3 methods. See the documentation for  [show()], [print()], [fitted()], [predict()], [plot()].
#'
#' @examples
#' ## Sample 75% of dyads in  French political Blogosphere's network data
#' adjMatrix <- missSBM::frenchblog2007 %>%
#'   igraph::as_adj (sparse = FALSE) %>%
#'   missSBM::sample(sampling = "dyad", parameters = 0.25)
#' collection <- missSBM::estimate(adjMatrix, 3:5, sampling = "dyad")
#' my_missSBM_fit <- collection$bestModel
#' class(my_missSBM_fit)
#' plot(my_missSBM_fit, "network")
#'
#' @include SBM_fit-Class.R
#' @include networkSampling_fit-Class.R
#' @export
missSBM_fit <-
  R6::R6Class(classname = "missSBM_fit",
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## fields for internal use (referring to mathematical notations)
  private = list(
    sampledNet = NULL, # network data with convenient encoding (object of class 'sampledNetwork')
    imputedNet = NULL, # imputed network data (a matrix possibly with NA when MAR sampling is used)
    sampling   = NULL, # fit of the current sampling model (object of class 'networkSampling_fit')
    SBM        = NULL, # fit of the current stochastic block model (object of class 'SBM_fit')
    optStatus  = NULL  # status of the optimization process
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description constructor for networkSampling
    #' @param sampledNet An object with class [`sampledNetwork`], typically obtained with the function [`prepare_data`] (real-word data) or [`sample`] (simulation).
    #' @param nBlocks integer, the number of blocks in the SBM
    #' @param netSampling The sampling design for the modelling of missing data: MAR designs ("dyad", "node") and NMAR designs ("double-standard", "block-dyad", "block-node" ,"degree")
    #' @param clusterInit Initial clustering: either a character in "hierarchical", "spectral" or "kmeans", or a vector with size \code{ncol(adjacencyMatrix)}, providing a user-defined clustering with \code{nBlocks} levels. Default is "hierarchical".
    #' @param useCov logical. If covariates are present in sampledNet, should they be used for the inference or of the network sampling design, or just for the SBM inference? default is TRUE.
    initialize = function(sampledNet, nBlocks, netSampling, clusterInit, useCov) {

      ## Basic sanity checks
      stopifnot(netSampling %in% available_samplings)
      stopifnot(inherits(sampledNet, "sampledNetwork"))
      stopifnot(length(nBlocks) == 1 & nBlocks >= 1 & is.numeric(nBlocks))

      ## Initial Clustering
      if (is.numeric(clusterInit) | is.factor(clusterInit))
        clusterInit <- as.integer(clusterInit)
      else
        clusterInit <- sampledNet$clustering(nBlocks, clusterInit)

      ## network data with basic imputation at start-up
      private$imputedNet <- sampledNet$imputation(clusterInit)

      ## Save the sampledNetwork object in the current environment
      private$sampledNet <- sampledNet

      ## Initialize the SBM fit
      if (is.null(sampledNet$covarArray) | !useCov) {
        private$SBM <- SBM_fit_nocovariate$new(private$imputedNet, clusterInit)
      } else {
        private$SBM <- SBM_fit_covariates$new(private$imputedNet, clusterInit, sampledNet$covarArray)
      }

      ## Initialize the sampling fit
      private$sampling <- switch(netSampling,
        "dyad"            = dyadSampling_fit$new(private$sampledNet),
        "node"            = nodeSampling_fit$new(private$sampledNet),
        "covar-dyad"      = covarDyadSampling_fit$new(private$sampledNet),
        "covar-node"      = covarNodeSampling_fit$new(private$sampledNet),
        "block-node"      = blockSampling_fit$new(private$sampledNet, Z),
        "double-standard" = doubleStandardSampling_fit$new(private$sampledNet),
        "block-dyad"      = blockDyadSampling_fit$new(private$sampledNet, Z),
        "degree"          = degreeSampling_fit$new(private$sampledNet, Z, private$SBM$connectParam),
        "snowball"        = nodeSampling_fit$new(private$sampledNet) # estimated sampling parameter not relevant
      )
    },
    #' @description a method to perform inference of the current missSBM fit with variational EM
    #' @param control a list of parameters controlling the variational EM algorithm. See details of function [`estimate`]
    doVEM = function(control = list(threshold = 1e-3, maxIter = 100, fixPointIter = 5, trace = 1)) {

      ## Initialization of quantities that monitor convergence
      delta     <- vector("numeric", control$maxIter)
      objective <- vector("numeric", control$maxIter)
      i <- 0; cond <- FALSE
      ## Starting the Variational EM algorithm
      if (control$trace) cat("\n Adjusting Variational EM for Stochastic Block Model\n")
      if (control$trace) cat("\n\tDyads are distributed according to a '", private$SBM$direction,"' SBM.\n", sep = "")
      if (control$trace) cat("\n\tImputation assumes a '", private$sampling$type,"' network-sampling process\n", sep = "")

      while (!cond) {
        i <- i + 1
        if (control$trace) cat(" iteration #:", i, "\r")
        pi_old <- private$SBM$connectParam # save current value of the parameters to assess convergence

        ## ______________________________________________________
        ## Variational E-Step
        #
        for (k in seq.int(control$fixPointIter)) {
          # update the variational parameters for missing entries (a.k.a nu)
          nu <- private$sampling$update_imputation(private$SBM$connectProb)
          private$imputedNet[private$sampledNet$NAs] <- nu[private$sampledNet$NAs]
          # update the variational parameters for block memberships (a.k.a tau)
          private$SBM$adjacencyMatrix <- private$imputedNet
          private$SBM$update_blocks(log_lambda = private$sampling$log_lambda)
        }

        ## ______________________________________________________
        ## M-step
        #
        # update the parameters of the SBM (a.k.a alpha and pi)
        private$SBM$update_parameters()
        # update the parameters of network sampling process (a.k.a psi)
        private$sampling$update_parameters(private$imputedNet, private$SBM$blocks)

        ## Check convergence
        delta[i] <- sqrt(sum((private$SBM$connectParam - pi_old)^2)) / sqrt(sum((pi_old)^2))
        cond     <- (i > control$maxIter) |  (delta[i] < control$threshold)
        objective[i] <- self$vBound

      }
      if (control$trace) cat("\n")
      private$optStatus <- data.frame(delta = delta[1:i], objective = c(NA, objective[2:i]), iteration = 1:i)
      invisible(private$optStatus)
    },
    #' @description show method for missSBM_fit
    show = function() {
      cat("missSBM-fit\n")
      cat("==================================================================\n")
      cat("Structure for storing a SBM fitted under missing data condition   \n")
      cat("==================================================================\n")
      cat("* Useful fields (most are special object themselves with methods) \n")
      cat("  $fittedSBM (the adjusted stochastic block model)                \n")
      cat("  $fittedSampling (the estimated sampling process)                \n")
      cat("  $imputedNetwork (the adjacency matrix with imputed values)      \n")
      cat("  $monitoring, $vICL, $vBound, $vExpec, $penalty                  \n")
      cat("* S3 methods                                                      \n")
      cat("  plot, coef, fitted, predict, print                              \n")
    },
    #' @description User friendly print method
    print = function() { self$show() }
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## ACTIVE BINDING
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field fittedSBM the fitted SBM, inheriting from class [`SBM_fit`]
    fittedSBM = function(value) {private$SBM},
    #' @field fittedSampling  the fitted sampling, inheriting from class [`networkSampling_fit`]
    fittedSampling = function(value) {private$sampling}  ,
    #' @field sampledNetwork The original network data used for the fit, with class [`sampledNetwork`]
    sampledNetwork = function(value) {private$sampledNet},
    #' @field imputedNetwork The network data with NAs values imputed with the model.
    imputedNetwork = function(value) {private$imputedNet},
    #' @field monitoring a list carrying information about the optimization process
    monitoring     = function(value) {private$optStatus},
    #' @field entropyImputed the entropy of the distribution of the imputed dyads
    entropyImputed = function(value) {
      nu <- private$imputedNet[private$sampledNet$missingDyads]
      res <- -sum(xlogx(nu) + xlogx(1 - nu))
      res
    },
    #' @field vBound double: approximation of the log-likelihood (variational lower bound) reached
    vBound  = function(value) {private$SBM$vBound + self$entropyImputed + private$sampling$vExpec},
    #' @field vExpec double: variational expectation of the complete log-likelihood
    vExpec  = function(value) {private$SBM$vExpec + private$sampling$vExpec},
    #' @field penalty double, value of the penalty term in vICL
    penalty = function(value) {private$SBM$penalty + private$sampling$penalty},
    #' @field vICL double: value of the integrated classification log-likelihood
    vICL    = function(value) {-2 * self$vExpec + self$penalty}
  )
)

## PUBLIC S3 METHODS FOR missSBMfit
## =========================================================================================

## Auxiliary functions to check the given class of an objet
is_missSBMfit <- function(Robject) {inherits(Robject, "missSBM_fit")}

#' Extract model fitted values from object  ['missSBM_fit'], return by [estimate()]
#'
#' @name fitted.missSBM_fit
#'
#' @param object an R6 object with class [`missSBM_fit`]
#' @param ... additional parameters for S3 compatibility.
#' @return a matrix of estimated probability of connection
#'
#' @export
fitted.missSBM_fit <- function(object, ...) {
  stopifnot(is_missSBMfit(object))
  fitted(object$fittedSBM)
}

#' Prediction of a ['missSBM_fit'] (i.e. network with imputed missing dyads)
#'
#' @name predicted.missSBM_fit
#'
#' @param object an R6 object with class [`missSBM_fit`]
#' @param ... additional parameters for S3 compatibility.
#'
#' @return an adjacency matrix between pairs of nodes. Missing dyads are imputed with
#' their expected values, i.e. by there estimated probabilities of connection under the missing SBM.
#'
#' @export
predict.missSBM_fit <- function(object, ...) {
  stopifnot(is_missSBMfit(object))
  object$imputedNetwork
}

#' Summary method for a ['missSBM_fit']
#'
#' @name summary.missSBM_fit
#'
#' @param object an R6 object with class [`missSBM_fit`]
#' @param ... additional parameters for S3 compatibility.
#'
#' @export
summary.missSBM_fit <- function(object, ...) {
  stopifnot(is_missSBMfit(object))
  object$show()
}

#' Visualization for an object [`missSBM_fit`]
#'
#' @description Plot function for the various fields of a [`missSBM_fit`]: the fitted SBM (network or connectivity),
#' the original sampled network data, and a plot monitoring the optimization.
#'
#' @name plot.missSBM_fit
#'
#' @param x an object with class [`missSBM_fit`]
#' @param type the type specifies the field to plot, either "network", "connectivity", "sampled" of "monitoring"
#' @param ... additional parameters for S3 compatibility. Not used
#' @export
#' @import ggplot2
#' @importFrom rlang .data
plot.missSBM_fit <- function(x, type = c("network", "connectivity", "sampled", "monitoring"), ...) {
  stopifnot(is_missSBMfit(x))
  type <- match.arg(type)
  if (type == "network")
    x$fittedSBM$plot("network")
  if (type == "connectivity")
    x$fittedSBM$plot("connectivity")
  if (type == "sampled")
    x$sampledNetwork$plot(x$fittedSBM$memberships)
  if (type == "monitoring")
    ggplot(x$monitoring, aes(x = .data$iteration, y = .data$objective)) + geom_line() + theme_bw()
}

#' Extract model coefficients
#'
#' @description Extracts model coefficients from objects [`missSBM_fit`] returned by [estimate()]
#'
#' @name coef.missSBM_fit
#'
#' @param object an R6 object with class [`missSBM_fit`]
#' @param type type of parameter that should be extracted. Either "mixture" (default), "connectivity", "covariates" or "sampling"
#' @param ... additional parameters for S3 compatibility. Not used
#' @return A vector or matrix of coefficients extracted from the missSBM_fit model.
#'
#' @export
coef.missSBM_fit <- function(object, type = c("mixture", "connectivity", "covariates", "sampling"), ...) {
  stopifnot(is_missSBMfit(object))
  switch(match.arg(type),
         mixture      = object$fittedSBM$mixtureParam,
         connectivity = object$fittedSBM$connectParam,
         covariates   = object$fittedSBM$covarParam,
         sampling     = object$fittedSampling$parameters)
}
