#' An object to represent a collection of missSBM_fit
#'
#' This R6 class stores a collection of missSBM_fit. Comes with basic printing methods and field access.
#'
#' @seealso The function [`estimate`], which produces an instance of this class.
#' The function [`smooth`] can be used to smooth the ICL on a collection of model, as post-treatment.
#'
#' @include missSBM_fit-Class.R
#' @export
missSBM_collection <-
  R6::R6Class(classname = "missSBM_collection",
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## fields for internal use (referring to mathematical notations)
  private = list(
    missSBM_fit = NULL, # a list of models
    # method for performing forward smoothing of the ICL
    # a list of parameters controlling the variational EM algorithm. See details of function [`estimate`]
    smoothing_forward = function(control) {
      trace <- control$trace > 0; control$trace <- FALSE
      sampledNet  <- private$missSBM_fit[[1]]$sampledNetwork
      sampling    <- private$missSBM_fit[[1]]$fittedSampling$type
      useCov      <- private$missSBM_fit[[1]]$useCovariates
      adjacencyMatrix <- sampledNet$adjacencyMatrix
      if (!is.null(sampledNet$covarArray)) {
        y <- as.vector(adjacencyMatrix)
        X <- apply(sampledNet$covarArray, 3, as.vector)
        adjacencyMatrix <- matrix(NA, sampledNet$nNodes, sampledNet$nNodes)
        NAs <- is.na(y)
        adjacencyMatrix[!NAs] <- logistic(residuals(glm.fit(X[!NAs, ], y[!NAs], family = binomial())))
      }

      if (trace) cat("   Going forward ")
      for (i in self$vBlocks[-length(self$vBlocks)]) {
        if (trace) cat("+")
        cl0 <- private$missSBM_fit[[i]]$fittedSBM$memberships
        if (length(unique(cl0)) == i) { # when would this not happens ?
          candidates <- mclapply(1:i, function(j) {
            cl <- cl0
            J  <- which(cl == j)
            if (length(J) > 1) {
              J1 <- base::sample(J, floor(length(J)/2))
              J2 <- setdiff(J, J1)
              cl[J1] <- j; cl[J2] <- i + 1
              model <- missSBM_fit$new(sampledNet, i + 1, sampling, cl, useCov)
              model$doVEM(control)
            } else {
              model <- private$missSBM_fit[[i + 1]]$clone()
            }
            model
          }, mc.cores = control$cores)
          best_one <- candidates[[which.min(sapply(candidates, function(candidate) candidate$vICL))]]
          if (is.na(private$missSBM_fit[[i + 1]]$vICL)) {
            private$missSBM_fit[[i + 1]] <- best_one
          } else if (best_one$vICL < private$missSBM_fit[[i + 1]]$vICL) {
            private$missSBM_fit[[i + 1]] <- best_one
          }

          # best_one <- candidates[[which.max(sapply(candidates, function(candidate) candidate$vBound))]]
          # if (is.na(private$missSBM_fit[[i + 1]]$vBound)) {
          #   private$missSBM_fit[[i + 1]] <- best_one
          # } else if (best_one$vBound > private$missSBM_fit[[i + 1]]$vBound) {
          #   private$missSBM_fit[[i + 1]] <- best_one
          # }
        }
      }
      if (trace) cat("\r                                                                                                    \r")
    },
    # method for performing backward smoothing of the ICL
    # control a list of parameters controlling the variational EM algorithm. See details of function [`estimate`]
    smoothing_backward = function(control) {

      trace <- control$trace > 0; control$trace <- FALSE
      sampledNet  <- private$missSBM_fit[[1]]$sampledNetwork
      sampling    <- private$missSBM_fit[[1]]$fittedSampling$type
      useCov      <- private$missSBM_fit[[1]]$useCovariates

      if (trace) cat("   Going backward ")
      for (i in rev(self$vBlocks[-1])) {
        if (trace) cat('+')
        cl0 <- factor(private$missSBM_fit[[i]]$fittedSBM$memberships)
        if (nlevels(cl0) == i) {
          candidates <- mclapply(combn(i, 2, simplify = FALSE), function(couple) {
            cl_fusion <- cl0
            levels(cl_fusion)[which(levels(cl_fusion) == paste(couple[1]))] <- paste(couple[2])
            levels(cl_fusion) <- as.character(1:(i - 1))
            model <- missSBM_fit$new(sampledNet, i - 1, sampling, cl_fusion, useCov)
            model$doVEM(control)
            model
          }, mc.cores = control$cores)

          best_one <- candidates[[which.min(sapply(candidates, function(candidate) candidate$vICL))]]
          if (is.na(private$missSBM_fit[[i - 1]]$vICL)) {
            private$missSBM_fit[[i - 1]] <- best_one
          } else if (best_one$vICL < private$missSBM_fit[[i - 1]]$vICL) {
            private$missSBM_fit[[i - 1]] <- best_one
          }

          # best_one <- candidates[[which.max(sapply(candidates, function(candidate) candidate$vBound))]]
          # if (is.na(private$missSBM_fit[[i - 1]]$vBound)) {
          #   private$missSBM_fit[[i - 1]] <- best_one
          # } else if (best_one$vBound > private$missSBM_fit[[i - 1]]$vBound) {
          #   private$missSBM_fit[[i - 1]] <- best_one
          # }
        }
      }
      if (trace) cat("\r                                                                                                    \r")
    }
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description constructor for networkSampling
    #' @param sampledNet An object with class [`sampledNetwork`], typically obtained with the function [`prepare_data`] (real-word data) or [`sample`] (simulation).
    #' @param vBlocks vector of integer with the number of blocks in the successively fitted models
    #' @param sampling The sampling design for the modelling of missing data: MAR designs ("dyad", "node") and NMAR designs ("double-standard", "block-dyad", "block-node" ,"degree")
    #' @param clusterInit Initial method for clustering: either a character in "hierarchical", "spectral" or "kmeans", or a list with \code{length(vBlocks)} vectors, each with size \code{ncol(adjacencyMatrix)}, providing a user-defined clustering. Default is "hierarchical".
    #' @param cores integer for number of cores used. Default is 1.
    #' @param trace integer for verbosity (0, 1, 2). Default is 1. Useless when \code{cores} > 1
    #' @param useCov logical. If covariates are present in sampledNet, should they be used for the inference or of the network sampling design, or just for the SBM inference? default is TRUE.
    initialize = function(sampledNet, vBlocks, sampling, clusterInit, cores, trace, useCov) {

      if (trace) cat("\n")
      if (trace) cat("\n Adjusting Variational EM for Stochastic Block Model\n")
      if (trace) cat("\n\tImputation assumes a '", sampling,"' network-sampling process\n", sep = "")
      if (trace) cat("\n")
      if (!is.list(clusterInit)) clusterInit <- rep(list(clusterInit), length(vBlocks))

      private$missSBM_fit <- mcmapply(
        function(nBlock, cl0) {
          if (trace) cat(" Initialization of model with", nBlock,"blocks.", "\r")
          missSBM_fit$new(sampledNet, nBlock, sampling, cl0, useCov)
        }, nBlock = vBlocks, cl0 = clusterInit, mc.cores = cores
      )
    },
    #' @description method to launch the estimation of the collection of models
    #' @param control a list of parameters controlling the variational EM algorithm. See details of function [`estimate`]
    estimate = function(control) {
      trace_main <- control$trace > 0
      control$trace <- ifelse (control$trace > 1, TRUE, FALSE)
      if (trace_main) cat("\n")
      invisible(
        mclapply(private$missSBM_fit, function(model) {
          if (trace_main) cat(" Performing VEM inference for model with", model$fittedSBM$nBlocks,"blocks.\r")
            model$doVEM(control)
          }, mc.cores = control$cores
        )
      )
      invisible(self)
    },
    #' @description method for performing smoothing of the ICL
    #' @param type character, the type of smoothing: forward, backward, both
    #' @param control a list of parameters controlling the smoothing. See details of regular function [`smooth`]
    smooth = function(type, control) {
      if (control$trace > 0) control$trace <- TRUE else control$trace <- FALSE
      if (control$trace) cat("\n Smoothing ICL\n")
      for (i in 1:control$iterates) {
        if (type %in% c('forward' , 'both')) private$smoothing_forward(control)
        if (type %in% c('backward', 'both')) private$smoothing_backward(control)
      }
    },
    #' @description show method for missSBM_collection
    show = function() {
      cat("--------------------------------------------------------\n")
      cat("COLLECTION OF", length(self$vBlocks), "SBM fits          \n")
      cat("========================================================\n")
      cat(" - Number of blocks considers: from ", min(self$vBlocks), " to ", max(self$vBlocks),"\n", sep = "")
      cat(" - Best model (smaller ICL): ", self$bestModel$fittedSBM$nBlocks, "\n", sep = "")
      cat(" - Fields: $models, $ICL, $vBlocks, $bestModel, $optimizationStatus\n")
    }
  ),
  active = list(
    #' @field models a list of models
    models = function(value) (private$missSBM_fit),
    #' @field ICL the vector of Integrated Classification Criterion (ICL) associated to the models in the collection (the smaller, the better)
    ICL = function(value) {setNames(sapply(self$models, function(model) model$vICL), names(self$vBlocks))},
    #' @field bestModel the best model according to the ICL
    bestModel = function(value) {self$models[[which.min(self$ICL)]]},
    #' @field vBlocks a vector with the number of blocks
    vBlocks = function(value) {sapply(self$models, function(model) model$fittedSBM$nBlocks)},
    #' @field optimizationStatus a data.frame summarizing the optimization process for all models
    optimizationStatus = function(value) {
      Reduce("rbind",
             lapply(self$models, function(model) {
               res <- model$monitoring
               res$nBlock <- model$fittedSBM$nBlocks
               res
             })
      )
    }
  )
)

#' Smooth the path ICL in a collection of missSBM_fit models
#'
#' Apply a split and/or merge strategy of the clustering in a path of models in a collection
#' of SBM ordered by number of block. The goal is to find better initialization. This results
#' in a "smoothing" of the ICL, that should be close to concave.
#'
#' @param Robject an object with class missSBM_collection, i.e. an output from \code{\link{estimate}}
#' @param type character indicating what kind of ICL smoothing should be use among "forward", "backward" or "both". Default is "forward".
#' @param control a list controlling the variational EM algorithm. See details.
#'
#' @details The list of parameters \code{control} controls the optimization process and the variational EM algorithm, with the following entries
#'  \itemize{
#'  \item{"iterates"}{integer for the number of iteration of smoothing. Default is 1.}
#'  \item{"threshold"}{stop when an optimization step changes the objective function by less than threshold. Default is 1e-4.}
#'  \item{"maxIter"}{V-EM algorithm stops when the number of iteration exceeds maxIter. Default is 200}
#'  \item{"fixPointIter"}{number of fix-point iteration for the Variational E step. Default is 5.}
#'  \item{"cores"}{integer for number of cores used. Default is 1.}
#'  \item{"trace"}{integer for verbosity. Useless when \code{cores} > 1}
#' }
#' @return an invisible missSBM_collection, in which the ICL has been smoothed
#' @export
smooth <- function(Robject, type = c("forward", "backward", "both"), control = list()) {

  stopifnot(inherits(Robject, "missSBM_collection"))

  ## defaut control parameter for VEM, overwritten by user specification
  ctrl <- list(threshold = 1e-3, maxIter = 50, fixPointIter = 2, cores = 1, trace = 1, iterates = 1)
  ctrl[names(control)] <- control

  ## Run the smoothing
  Robject$smooth(match.arg(type), ctrl)

  invisible(Robject)
}
