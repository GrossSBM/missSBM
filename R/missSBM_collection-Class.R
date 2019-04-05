#' An R6 Class to represent a collection of missSBM_fit
#'
#' The function \code{\link{estimate}} produces an instance of this class.
#' The function \code{\link{smooth}} (also available as ann R6 method of this class) can be used
#' to smooth the ICL on a collection of model, as post-treatment.
#'
#' @field models a list of models
#' @field ICL the vector of Integrated Classfication Criterion (ICL) associated to the models of the collection: the smaller, the better
#' @field bestModel the best model according to the ICL
#' @field optimizationStatus a data.frame summarizing the optimization process for all models
#'
#' @seealso \code{\link{estimate}}, \code{\link{smooth}}
#' @include missSBM_fit-Class.R
#' @export
missSBM_collection <-
  R6::R6Class(classname = "missSBM_collection",
  private = list(missSBM_fit = NULL), # a list of models
  active = list(
    models = function(value) (private$missSBM_fit),
    ICL = function(value) {setNames(sapply(self$models, function(model) model$vICL), names(self$vBlocks))},
    bestModel = function(value) {self$models[[which.min(self$ICL)]]},
    vBlocks = function(value) {sapply(self$models, function(model) model$fittedSBM$nBlocks)},
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

missSBM_collection$set("public", "initialize",
function(adjMatrix, vBlocks, sampling, clusterInit, covarMatrix, covarSimilarity, cores, trace) {

  if (trace) cat("\n")
  if (trace) cat("\n Adjusting Variational EM for Stochastic Block Model\n")
  if (trace) cat("\n\tImputation assumes a '", sampling,"' network-sampling process\n", sep = "")
  if (trace) cat("\n")
  if (!is.list(clusterInit)) clusterInit <- rep(list(clusterInit), length(vBlocks))

  covarArray <- getCovarArray(covarMatrix, covarSimilarity)
  sampledNet <- sampledNetwork$new(adjMatrix)
  private$missSBM_fit <- mcmapply(
    function(nBlock, cl0) {
      if (trace) cat(" Initialization of model with", nBlock,"blocks.", "\r")
      missSBM_fit$new(sampledNet, nBlock, sampling, cl0, covarMatrix, covarArray)
    }, nBlock = vBlocks, cl0 = clusterInit, mc.cores = cores
  )
})

missSBM_collection$set("public", "estimate",
function(control_VEM, mc.cores, trace) {
  if (trace) cat("\n")
  invisible(
    mclapply(private$missSBM_fit,
       function(model) {
         if (trace) cat(" Performing VEM inference for model with", model$fittedSBM$nBlocks,"blocks.\r")
           model$doVEM(control_VEM)
       }, mc.cores = mc.cores
    )
  )
  invisible(self)
})

missSBM_collection$set("public", "smooth",
function(type, control) {
  if (control$trace) cat("\n Smoothing ICL\n")
  if (type == "forward")
    private$smoothing_forward(control)
  if (type == "backward")
    private$smoothing_backward(control)
  if (type == "both")
    for (i in 1:control$iterates) {
      private$smoothing_forward(control)
      private$smoothing_backward(control)
    }
})

#' Smooth path of models in a collection
#'
#' Apply a split and/or merge strategy to the path of model in a collection of SBM, in order to find better initialization. This should result in
#' a "smoothing" of the ICL, that should be close to concave.
#'
#' @param Robject an object with class missSBM_collection, i.e. an output from \code{\link{estimate}}
#' @param type character indicating what kind of ICL smoothing should be use among "forward", "backward" or "both". Default is "foward".
#' @param control_VEM a list controlling the variational EM algorithm. See details in \code{\link{estimate}}.
#' @param cores integer, the number of cores to use when multiply model are fitted. Default is 1.
#' @param iterates integer for the number of iteration in case of foward-backward (aka both) smoothing. Default is 1.
#' @param trace logical, control the verbosity. Default to \code{TRUE}.
#'
#' @return an invisible missSBM_collection, in which the ICL has been smoothed
#' @export
smooth <- function(Robject, type = c("forward", "backward", "both"), control_VEM = list(), cores = 1, iterates = 1, trace = TRUE) {

  stopifnot(inherits(Robject, "missSBM_collection"))

  ## defaut control parameter for VEM, overwritten by user specification
  control <- list(threshold = 1e-4, maxIter = 100, fixPointIter = 5, trace = trace)
  control[names(control_VEM)] <- control_VEM
  ## add some additional control param to pass to smoothing
  control$iterates <- iterates
  control$mc.cores <- cores
  ## Run the smoothing
  Robject$smooth(match.arg(type), control)

  invisible(Robject)
}

missSBM_collection$set("private", "smoothing_forward",
function(control) {
  trace <- control$trace; control$trace <- FALSE
  sampledNet  <- private$missSBM_fit[[1]]$sampledNetwork
  sampling    <- private$missSBM_fit[[1]]$fittedSampling$type
  covarMatrix <- sampledNet$covarMatrix
  covarArray  <- private$missSBM_fit[[1]]$fittedSBM$covarArray
  adjacencyMatrix <- sampledNet$adjMatrix
  if (!is.null(covarArray)) {
    y <- as.vector(adjacencyMatrix)
    X <- apply(covarArray, 3, as.vector)
    adjacencyMatrix <- matrix(NA, N, N)
    NAs <- is.na(y)
    adjacencyMatrix[!NAs] <- logistic(residuals(glm.fit(X[!NAs, ], y[!NAs], family = binomial())))
  }

  if (trace) cat("   Going forward ")
  for (i in self$vBlocks[-length(self$vBlocks)]) {
    if (trace) cat("+")
    cl_split <- factor(private$missSBM_fit[[i]]$fittedSBM$memberships)
    levels(cl_split) <- c(levels(cl_split), as.character(i + 1))
    if (nlevels(cl_split) - 1 == i) { # when would this happens ?
      candidates <- mclapply(1:i, function(j) {
        cl <- as.numeric(cl_split); J <- which(cl == j)
        cut <- as.numeric(init_hierarchical(adjacencyMatrix[J, ], 2))
        if (length(unique(cut)) > 1) {
          cut <- as.numeric(init_hierarchical(adjacencyMatrix[J, ], 2))
          cl[J][which(cut == 1)] <- j; cl[J][which(cut == 2)] <- i + 1
          model <- missSBM_fit$new(sampledNet, i + 1, sampling, cl, covarMatrix, covarArray)
          model$doVEM(control)
          model
        } else {
          private$missSBM_fit[[i + 1]]$clone()
        }
      }, mc.cores = control$mc.cores)

      best_one <- candidates[[which.max(sapply(candidates, function(candidate) candidate$vBound))]]
      if (best_one$vBound > private$missSBM_fit[[i + 1]]$vBound) {
        private$missSBM_fit[[i + 1]] <- best_one
      }
    }
  }
  if (trace) cat("\r                                                                                                    \r")
})

missSBM_collection$set("private", "smoothing_backward",
function(control) {
  trace <- control$trace; control$trace <- FALSE
  sampledNet  <- private$missSBM_fit[[1]]$sampledNetwork
  sampling    <- private$missSBM_fit[[1]]$fittedSampling$type
  covarMatrix <- sampledNet$covarMatrix
  covarArray  <- private$missSBM_fit[[1]]$fittedSBM$covarArray

  if (trace) cat("   Going backward ")
  for (i in rev(self$vBlocks[-1])) {
    if (trace) cat('+')
    cl0 <- factor(private$missSBM_fit[[i]]$fittedSBM$memberships)
    if (nlevels(cl0) == i) {
      candidates <- mclapply(combn(i, 2, simplify = FALSE), function(couple) {
        cl_fusion <- cl0
        levels(cl_fusion)[which(levels(cl_fusion) == paste(couple[1]))] <- paste(couple[2])
        levels(cl_fusion) <- as.character(1:(i - 1))
        model <- missSBM_fit$new(sampledNet, i - 1, sampling, cl_fusion, covarMatrix, covarArray)
        model$doVEM(control)
        model
      }, mc.cores = control$mc.cores)

      best_one <- candidates[[which.max(sapply(candidates, function(candidate) candidate$vBound))]]
      if (best_one$vBound > private$missSBM_fit[[i - 1]]$vBound) {
        private$missSBM_fit[[i - 1]] <- best_one
      }
    }
  }
  if (trace) cat("\r                                                                                                    \r")
})

missSBM_collection$set("public", "show",
function() {
  cat("--------------------------------------------------------\n")
  cat("COLLECTION OF", length(self$vBlocks), "SBM fits          \n")
  cat("========================================================\n")
  cat(" - Number of blocks considers: from ", min(self$vBlocks), " to ", max(self$vBlocks),"\n", sep = "")
  cat(" - Best model (smaller ICL): ", self$bestModel$fittedSBM$nBlocks, "\n", sep = "")
  cat(" - Fields: $models, $ICL, $vBlocks, $bestModel, $optimizationStatus\n")
})
missSBM_collection$set("public", "print", function() self$show())

