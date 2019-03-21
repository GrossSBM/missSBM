#' An R6 Class to represent a collection of missSBM_fit
#'
#' The function \code{\link{estimate}} produces an instance of this class.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for \code{\link[=smoothing_ICL]{smoothing_ICL}},
#' \code{\link[=plot.missSBM_collection]{plot}}.
#'
#' @field ICL the vector of Integrated Classfication Criterion (ICL) associated to the models of the collection: the smaller, the better
#' @field bestModel the best model according to the ICL
#' @field optimizationStatus a data.frame summarizing the optimization process for all models
#'
#' @seealso \code{\link{estimate}}
#' @include utils_smoothing.R
#' @include missingSBM_fit-Class.R
#' @export
missSBM_collection <-
  R6::R6Class(classname = "missSBM_collection",
  public = list(
    models = NULL  # a list of models
  ),
  active = list(
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
  self$models <- mcmapply(
    function(nBlock, cl0) {
      if (trace) cat(" Initialization of model with", nBlock,"blocks.", "\r")
      missingSBM_fit$new(sampledNet, nBlock, sampling, cl0, covarMatrix, covarArray)
    }, nBlock = vBlocks, cl0 = clusterInit, mc.cores = cores
  )
})

missSBM_collection$set("public", "estimate",
function(control_VEM, mc.cores, trace) {
  if (trace) cat("\n")
  invisible(
    mclapply(self$models,
       function(model) {
         if (trace) cat(" Performing VEM inference for model with", model$fittedSBM$nBlocks,"blocks.\r")
           model$doVEM(control_VEM)
       }, mc.cores = mc.cores
    )
  )
  invisible(self)
})

missSBM_collection$set("public", "smooth_ICL",
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
#' @param split character indicating the function use to split nodes during the forward algorithm. Either "hierarchical", "spectral" or "kmeans". Default is "hierarchical".
#' @param control_VEM a list controlling the variational EM algorithm. See details in \code{\link{estimate}}.
#' @param cores integer, the number of cores to use when multiply model are fitted. Default is 2.
#' @param iterates integer for the number of iteration in case of foward-backward (aka both) smoothing. Default is 1.
#' @param trace logical, control the verbosity. Default to \code{TRUE}.
#'
#' @return an invisible missSBM_collection, in which the ICL has been smoothed
#' @export
smooth <- function(Robject, type = c("forward", "backward", "both"), split = c("hierarchical", "spectral", "kmeans"), control_VEM = list(), cores = 2, iterates = 1, trace = TRUE) {

  stopifnot(inherits(Robject, "missSBM_collection"))

  ## defaut control parameter for VEM, overwritten by user specification
  control <- list(threshold = 1e-4, maxIter = 100, fixPointIter = 5, trace = trace)
  control[names(control_VEM)] <- control_VEM
  ## add some additional control param to pass to smoothing
  control$iterates <- iterates
  control$mc.cores <- cores
  ## select which clustering method will be used for splitting a group
  control$split_fn <- switch(match.arg(split),
    "spectral"     = init_spectral,
    "hierarchical" = init_hierarchical,
    "kmeans"       = init_kmeans
  )

  Robject$smooth_ICL(match.arg(type), control)

  invisible(Robject)
}

missSBM_collection$set("private", "smoothing_forward",
function(control) {
  trace <- control$trace; control$trace <- FALSE
  sampledNet  <- self$models[[1]]$sampledNetwork
  sampling    <- self$models[[1]]$fittedSampling$type
  covarMatrix <- sampledNet$covarMatrix
  covarArray  <- self$models[[1]]$fittedSBM$covarArray

  if (trace) cat("   Going forward ")
  for (i in self$vBlocks[-length(self$vBlocks)]) {
    if (trace) cat("+")
    cl_split <- factor(self$models[[i]]$fittedSBM$memberships)
    levels(cl_split) <- c(levels(cl_split), as.character(i + 1))
    if (nlevels(cl_split) - 1 == i) {
      candidates <- mclapply(1:i, function(j) {
        cl <- as.numeric(cl_split); indices <- which(cl == j)
        if (length(cl[indices]) > 10) { ## ???? What is this "10"
          cut <- as.numeric(control$split_fn(sampledNet$adjMatrix[indices, indices],2))
          cl[which(cl == j)][which(cut == 1)] <- j; cl[which(cl == j)][which(cut == 2)] <- i + 1
          model <- missingSBM_fit$new(sampledNet, i + 1, sampling, cl, covarMatrix, covarArray)
          model$doVEM(control)
          model
        } else {
          self$models[[i + 1]]
        }
      }, mc.cores = control$mc.cores)
      vICLs <- sapply(candidates, function(candidate) candidate$vICL)
      best_one <- candidates[[which.min(vICLs)]]
      if (is.na(self$models[[i + 1]]$vICL)) {
        self$models[[i + 1]] <- best_one
      } else if (best_one$vICL < self$models[[i + 1]]$vICL) {
        self$models[[i + 1]] <- best_one
      }
    }
  }
  if (trace) cat("\r                                                                                                    \r")
})

missSBM_collection$set("private", "smoothing_backward",
function(control) {
  trace <- control$trace; control$trace <- FALSE
  sampledNet  <- self$models[[1]]$sampledNetwork
  sampling    <- self$models[[1]]$fittedSampling$type
  covarMatrix <- sampledNet$covarMatrix
  covarArray  <- self$models[[1]]$fittedSBM$covarArray

  if (trace) cat("   Going backward ")
  for (i in rev(self$vBlocks[-1])) {
    if (trace) cat('+')
    cl0 <- factor(self$models[[i]]$fittedSBM$memberships)
    if (nlevels(cl0) == i) {
      candidates <- mclapply(combn(i, 2, simplify = FALSE), function(couple) {
        cl_fusion <- cl0
        levels(cl_fusion)[which(levels(cl_fusion) == paste(couple[1]))] <- paste(couple[2])
        levels(cl_fusion) <- as.character(1:(i - 1))
        model <- missingSBM_fit$new(sampledNet, i - 1, sampling, cl_fusion, covarMatrix, covarArray)
        model$doVEM(control)
        model
      }, mc.cores = control$mc.cores)
      vICLs <- sapply(candidates, function(candidate) candidate$vICL)
      best_one <- candidates[[which.min(vICLs)]]
      if (is.na(self$models[[i - 1]]$vICL)) {
        self$models[[i - 1]] <- best_one
      } else if (best_one$vICL < self$models[[i - 1]]$vICL) {
        self$models[[i - 1]] <- best_one
      }
    }
  }
})

missSBM_collection$set("public", "show",
function() {
  cat("--------------------------------------------------------\n")
  cat("COLLECTION OF", length(self$models), "SBM fits          \n")
  cat("========================================================\n")
  cat(" - Number of blocks considers: from ", min(self$vBlocks), " to ", max(self$vBlocks),"\n", sep = "")
  cat(" - Best model (smaller ICL): ", self$bestModel$fittedSBM$nBlocks, "\n", sep = "")
})
missSBM_collection$set("public", "print", function() self$show())

