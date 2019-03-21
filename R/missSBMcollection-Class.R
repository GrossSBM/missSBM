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
    sampledNet   = NULL, # network data with convenient encoding (object of class 'sampledNetwork')
    vBlocks      = NULL, # vectors of block number
    sampling     = NULL, # the sampling design
    clusterInit  = NULL, # rule for initial clustering
    covarMatrix  = NULL, # the matrix of covariates
    covarArray   = NULL, # the array of covariates
    splitting_fn = NULL, # function use to split the clustering during the smoothing of the ICL
    models       = NULL  # a list of models
  ),
  active = list(
    ICL = function(value) {setNames(sapply(self$models, function(model) model$vICL), names(self$vBlocks))},
    bestModel = function(value) {self$models[[which.min(self$ICL)]]},
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
function(adjMatrix, vBlocks, sampling, clusterInit, covarMatrix, covarSimilarity, smoothing, mc.cores, trace) {

  ## Storage...
  self$sampledNet   <- sampledNetwork$new(adjMatrix)
  self$vBlocks      <- vBlocks
  self$sampling     <- sampling
  self$covarMatrix  <- covarMatrix
  self$covarArray   <- getCovarArray(covarMatrix, covarSimilarity)

  ## Function use to split the clusters durint the smoothing of the ICL
  if (!is.character(clusterInit)) {
    self$splitting_fn <- init_hierarchical
  } else {
    self$splitting_fn <- switch(clusterInit,
                                "spectral" = init_spectral,
                                "hierarchical" = init_hierarchical,
                                init_hierarchical)
  }
  if (!is.list(clusterInit)) clusterInit <- rep(list(clusterInit), length(vBlocks))
  self$clusterInit <- clusterInit

  if (trace) cat("\n")
  if (trace) cat("\n Adjusting Variational EM for Stochastic Block Model\n")
  if (trace) cat("\n\tImputation assumes a '", sampling,"' network-sampling process\n", sep = "")
  if (trace) cat("\n")
  self$models <- mcmapply(
    function(nBlock, clInit) {
      if (trace) cat(" Initialization of model with", nBlock,"blocks.", "\r")
      missingSBM_fit$new(self$sampledNet, nBlock, self$sampling, clInit, self$covarMatrix, self$covarArray)
    }, nBlock = vBlocks, clInit = clusterInit, mc.cores = mc.cores
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
function(type, control, mc.cores, iter, trace) {
  if (trace & type != "none") cat("\n Smoothing ICL\n")
  control$trace <- FALSE
  if (type == "forward")
    private$smoothing_forward(control, mc.cores, trace)
  if (type == "backward")
    private$smoothing_backward(control, mc.cores, trace)
  if (type == "both")
    for (i in 1:iter) {
      private$smoothing_forward(control, mc.cores, trace)
      private$smoothing_backward(control, mc.cores, trace)
    }
})

missSBM_collection$set("private", "smoothing_forward",
function(control, mc.cores, trace) {

  if (trace) cat("   Going forward ")
  for (i in self$vBlocks[-length(self$vBlocks)]) {
    if (trace) cat("+")
    cl_split <- factor(self$models[[i]]$fittedSBM$memberships)
    levels(cl_split) <- c(levels(cl_split), as.character(i + 1))
    if (nlevels(cl_split) - 1 == i) {
      candidates <- mclapply(1:i, function(j) {
        cl <- as.numeric(cl_split); indices <- which(cl == j)
        if (length(cl[indices]) > 10) { ## ???? What is this "10"
          cut <- as.numeric(self$splitting_fn(self$sampledNet$adjMatrix[indices, indices],2))
          cl[which(cl == j)][which(cut == 1)] <- j; cl[which(cl == j)][which(cut == 2)] <- i + 1
          model <- missingSBM_fit$new(self$sampledNet, i + 1, self$sampling, cl, self$covarMatrix, self$covarArray)
          model$doVEM(control)
          model
        } else {
          self$models[[i + 1]]
        }
      }, mc.cores = 1)
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
function(control, mc.cores, trace) {

  if (trace) cat("   Going backward ")
  for (i in rev(self$vBlocks[-1])) {
    if (trace) cat('+')
    cl0 <- factor(self$models[[i]]$fittedSBM$memberships)
    if (nlevels(cl0) == i) {
      candidates <- mclapply(combn(i, 2, simplify = FALSE), function(couple) {
        cl_fusion <- cl0
        levels(cl_fusion)[which(levels(cl_fusion) == paste(couple[1]))] <- paste(couple[2])
        levels(cl_fusion) <- as.character(1:(i - 1))
        model <- missingSBM_fit$new(self$sampledNet, i - 1, self$sampling, cl_fusion, self$covarMatrix, self$covarArray)
        model$doVEM(control)
        model
      }, mc.cores = mc.cores)
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

