#' An R6 class to represent an SBM fit with missing data
#'
#' @description The function [estimateMissSBM()] fits a collection of SBM for varying number of block.
#' Each fitted SBM is an instance of an R6 object with class [`missSBM_fit`], described here.
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
#'   missSBM::observeNetwork(sampling = "dyad", parameters = 0.75)
#' collection <- estimateMissSBM(adjMatrix, 3:5, sampling = "dyad")
#' my_missSBM_fit <- collection$bestModel
#' class(my_missSBM_fit)
#' plot(my_missSBM_fit, "imputed")
#'
#' @include R6Class-simpleSBM_fit.R
#' @include R6Class-networkSampling_fit.R
#' @export
missSBM_fit <-
  R6::R6Class(classname = "missSBM_fit",
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## fields for internal use (referring to mathematical notations)
  private = list(
    nu         = NULL, # imputed values (a sparse Matrix with entries only for imputed values, "dgCMatrix)
    sampling   = NULL, # fit of the current sampling model (object of class 'networkSampling_fit')
    SBM        = NULL, # fit of the current stochastic block model (object of class 'SBM_fit')
    optStatus  = NULL  # status of the optimization process
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description constructor for networkSampling
    #' @param partlyObservedNet An object with class [`partlyObservedNetwork`].
    #' @param netSampling The sampling design for the modelling of missing data: MAR designs ("dyad", "node") and NMAR designs ("double-standard", "block-dyad", "block-node" ,"degree")
    #' @param clusterInit Initial clustering: a vector with size \code{ncol(adjacencyMatrix)}, providing a user-defined clustering. The number of blocks is deduced from the number of levels in with \code{clusterInit}.
    #' @param useCov logical. If covariates are present in partlyObservedNet, should they be used for the inference or of the network sampling design, or just for the SBM inference? default is TRUE.
    initialize = function(partlyObservedNet, netSampling, clusterInit, useCov = TRUE) {

      ## Basic sanity checks
      stopifnot(netSampling %in% available_samplings)
      stopifnot(inherits(partlyObservedNet, "partlyObservedNetwork"))
      stopifnot(is.numeric(clusterInit))

      ## Initialize the SBM fit
      covariates <- array2list(partlyObservedNet$covarArray)
      if (!useCov) covariates <- list()

      if (length(covariates) == 0) {
        if (netSampling %in% c("dyad", "node", "snowball")) {
          private$SBM <- SimpleSBM_fit_noCov$new(partlyObservedNet, clusterInit)
        }
        if (netSampling %in% c("double-standard", "block-node", "block-dyad")) {
          private$SBM <- SimpleSBM_fit_NMAR$new(partlyObservedNet, clusterInit)
        }
      } else {
        private$SBM <- SimpleSBM_fit_withCov$new(partlyObservedNet, clusterInit, covariates)
      }

      ## Initialize the sampling fit
      private$sampling <- switch(netSampling,
        "dyad"            = dyadSampling_fit$new(partlyObservedNet),
        "node"            = nodeSampling_fit$new(partlyObservedNet),
        "covar-dyad"      = covarDyadSampling_fit$new(partlyObservedNet),
        "covar-node"      = covarNodeSampling_fit$new(partlyObservedNet),
        "double-standard" = doubleStandardSampling_fit$new(partlyObservedNet),
        "block-node"      = blockNodeSampling_fit$new(partlyObservedNet, clustering_indicator(clusterInit)),
        "block-dyad"      = blockDyadSampling_fit$new(partlyObservedNet, clustering_indicator(clusterInit)),
        "degree"          = degreeSampling_fit$new(partlyObservedNet, clustering_indicator(clusterInit), private$SBM$connectParam$mean),
        "snowball"        = nodeSampling_fit$new(partlyObservedNet) # estimated sampling parameter not relevant
      )
    },
    #' @description a method to perform inference of the current missSBM fit with variational EM
    #' @param control a list of parameters controlling the variational EM algorithm. See details of function [estimateMissSBM()]
    doVEM = function(control = list(threshold = 1e-2, maxIter = 100, fixPointIter = 3, trace = 1)) {

      ## Initialization of quantities that monitor convergence
      delta_par <- vector("numeric", control$maxIter); delta_par[1] <- NA
      delta_obj <- vector("numeric", control$maxIter); delta_obj[1] <- NA
      objective <- vector("numeric", control$maxIter); objective[1] <- self$loglik
      iterate <- 1; stop <- FALSE

      ## Starting the Variational EM algorithm
      if (control$trace) cat("\n Adjusting Variational EM for Stochastic Block Model\n")
      if (control$trace) cat("\n\tDyads are distributed according to a '",
                             ifelse(private$SBM$directed, "directed", "undirected"),"' SBM.\n", sep = "")
      if (control$trace) cat("\n\tImputation assumes a '", private$sampling$type,"' network-sampling process\n", sep = "")

      while (!stop) {
        iterate <- iterate + 1
        if (control$trace) cat(" iteration #:", iterate, "\r")

        theta_old <- private$SBM$connectParam # save current value of the parameters to assess convergence
        ## ______________________________________________________
        ## Variational E-Step
        #
        for (i in seq.int(control$fixPointIter)) {
          # update the variational parameters for missing entries (a.k.a nu)
          private$nu <- private$sampling$update_imputation(private$SBM$imputation)
          # update the variational parameters for block memberships (a.k.a tau)
          private$SBM$update_blocks(log_lambda = private$sampling$log_lambda)
        }

        ## ______________________________________________________
        ## M-step
        #
        # update the parameters of the SBM (a.k.a pi and theta)
        private$SBM$update_parameters(private$nu)
        # update the parameters of network sampling process (a.k.a psi)
        private$sampling$update_parameters(private$nu, private$SBM$probMemberships)

        ## Check convergence
        objective[iterate] <- self$loglik
        delta_par[iterate] <- sqrt(sum((private$SBM$connectParam$mean - theta_old$mean)^2)) / sqrt(sum((theta_old$mean)^2))
        delta_obj[iterate] <- abs(objective[iterate] - objective[iterate-1]) / abs(objective[iterate])
        stop <- (iterate > control$maxIter) |  ((delta_par[iterate] < control$threshold) & (delta_obj[iterate] < control$threshold))
      }
      private$SBM$reorder()

      if (control$trace) cat("\n")
      private$optStatus <- data.frame(iteration = 1:iterate, delta_pararameters = delta_par[1:iterate], delta_objective = delta_obj[1:iterate], elbo = objective[1:iterate])
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
      cat("  $monitoring, $ICL, $loglik, $vExpec, $penalty                  \n")
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
    #' @field fittedSBM the fitted SBM with class [`SimpleSBM_fit_noCov`], [`SimpleSBM_fit_withCov`] or
    #' [`SimpleSBM_fit_NMAR`] inheriting from class [`sbm::SimpleSBM_fit`]
    fittedSBM = function(value) {private$SBM},
    #' @field fittedSampling  the fitted sampling, inheriting from class [`networkSampling`] and corresponding fits
    fittedSampling = function(value) {private$sampling},
    #' @field imputedNetwork The network data as a matrix with NAs values imputed with the current model
    imputedNetwork = function(value) {
      res <- private$SBM$networkData + private$nu
      if (!private$SBM$directed) res <- res + t(res)
      res
    },
    #' @field monitoring a list carrying information about the optimization process
    monitoring     = function(value) {private$optStatus},
    #' @field entropyImputed the entropy of the distribution of the imputed dyads
    entropyImputed = function(value) {
      res <- -sum(xlogx(private$nu) + xlogx(1 - private$nu))
      res
    },
    #' @field entropy the entropy due to the distribution of the imputed dyads and of the clustering
    entropy = function(value) {private$SBM$entropy + self$entropyImputed},
    #' @field vExpec double: variational expectation of the complete log-likelihood
    vExpec  = function(value) {private$SBM$vExpec + private$sampling$vExpec},
    #' @field penalty double, value of the penalty term in ICL
    penalty = function(value) {private$SBM$penalty + private$sampling$penalty},
    #' @field loglik double: approximation of the log-likelihood (variational lower bound) reached
    loglik  = function(value) {self$vExpec + self$entropy},
    #' @field ICL double: value of the integrated classification log-likelihood
    ICL    = function(value) {-2 * self$vExpec + self$penalty}
  )
)

## PUBLIC S3 METHODS FOR missSBMfit
## =========================================================================================

## Auxiliary functions to check the given class of an objet
is_missSBMfit <- function(Robject) {inherits(Robject, "missSBM_fit")}

#' Extract model fitted values from object  [`missSBM_fit`], return by [estimateMissSBM()]
#'
#' @name fitted.missSBM_fit
#'
#' @param object an R6 object with class [`missSBM_fit`]
#' @param ... additional parameters for S3 compatibility.
#'
#' @return A matrix of estimated probabilities of connection
#'
#' @importFrom stats fitted
#' @export
fitted.missSBM_fit <- function(object, ...) {
  stopifnot(is_missSBMfit(object))
  fitted(object$fittedSBM)
}

#' Prediction of a [`missSBM_fit`] (i.e. network with imputed missing dyads)
#'
#' @name predicted.missSBM_fit
#'
#' @param object an R6 object with class [`missSBM_fit`]
#' @param ... additional parameters for S3 compatibility.
#'
#' @return an adjacency matrix between pairs of nodes. Missing dyads are imputed with
#' their expected values, i.e. by there estimated probabilities of connection under the missing SBM.
#'
#' @importFrom stats predict
#' @export
predict.missSBM_fit <- function(object, ...) {
  stopifnot(is_missSBMfit(object))
  object$imputedNetwork
}

#' Summary method for a [`missSBM_fit`]
#'
#' @name summary.missSBM_fit
#'
#' @param object an R6 object with class [`missSBM_fit`]
#' @param ... additional parameters for S3 compatibility.
#'
#' @return a basic printing output
#'
#' @method summary missSBM_fit
#' @export
summary.missSBM_fit <- function(object, ...) {
  stopifnot(is_missSBMfit(object))
  object$show()
}

#' Visualization for an object [`missSBM_fit`]
#'
#' @description Plot function for the various fields of a [`missSBM_fit`]: the fitted
#' SBM (network or connectivity), and a plot monitoring the optimization.
#'
#' @return a ggplot object
#'
#' @name plot.missSBM_fit
#'
#' @param x an object with class [`missSBM_fit`]
#' @param type the type specifies the field to plot, either "expected", "imputed", "meso",  or "monitoring"
#' @param ... additional parameters for S3 compatibility. Not used
#' @export
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom sbm plotMyMatrix
plot.missSBM_fit <- function(x, type = c("expected", "imputed", "meso", "monitoring"), ...) {
  stopifnot(is_missSBMfit(x))
  gg_obj <- switch(match.arg(type),
    "expected"   = x$fittedSBM$plot("expected"),
    "meso"       = x$fittedSBM$plot("meso"),
    "imputed"    = plotMyMatrix(as.matrix(predict(x)),  clustering = list(row = x$fittedSBM$memberships)),
    "monitoring" = ggplot(x$monitoring, aes(x = .data$iteration, y = .data$elbo)) + geom_line() + theme_bw()
  )
  gg_obj
}

#' Extract model coefficients
#'
#' @description Extracts model coefficients from objects [`missSBM_fit`] returned by [estimateMissSBM()]
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
         mixture      = object$fittedSBM$blockProp,
         connectivity = object$fittedSBM$connectParam,
         covariates   = object$fittedSBM$covarParam,
         sampling     = object$fittedSampling$parameters)
}
