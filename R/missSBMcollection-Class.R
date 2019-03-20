#' R6 Class definition of an missSBMcollection
#'
#' This class is to store a collection fo missSBM-fit
#'
#' @include utils_smoothing.R
#' @include missingSBM_fit-Class.R
#' @export
missSBM_collection <-
  R6::R6Class(classname = "missSBM_collection",
  public = list(
    sampledNet   = NULL, # network data with convenient encoding (object of class 'sampledNetwork')
    vBlocks      = NULL, # vectors of block number
    sampling     = NULL, # the sampling design
    covarMatrix  = NULL, # the matrix of covariates
    covarArray   = NULL, # the array of covariates
    splitting_fn = NULL, # function use to split the clustering during the smoothing of the ICL
    models       = NULL, # a list of models
    initialize   = function(adjMatrix, vBlocks, sampling, clusterInit, covarMatrix, covarSimilarity, smoothing, trace = FALSE) {

      ## Storage...
      self$sampledNet   <- sampledNetwork$new(adjMatrix)
      self$vBlocks      <- vBlocks
      self$sampling     <- sampling
      self$covarMatrix  <- covarMatrix
      self$covarArray   <- getCovarArray(covarMatrix, covarSimilarity)

      ## Function use to split the clusters durint the smoothing of the ICL
      if (!is.character(clusterInit)) {
        split_fn <- init_hierarchical
      } else {
        split_fn <- switch(clusterInit,
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

    },
    estimate = function(control_VEM, mc.cores, trace) {
      if (trace) cat("\n")
      mclapply(self$models,
        function(model) {
           if (trace) cat(" Performing VEM inference for model with", model$fittedSBM$nBlocks,"blocks.\r")
            model$doVEM(control_VEM)
        }, mc.cores = mc.cores
      )
    },
    smooth_ICL = function(type, control, mc.cores, iter, trace) {
      if (trace) cat("\n Smoothing ICL\n")
      control$trace <- FALSE
      if (type == "forward")
        self$smoothingForward(control, mc.cores, trace)
      if (type == "backward")
        self$smoothingBackward(control, mc.cores, trace)
      if (type == "booth")
        for (i in 1:iter) {
          self$smoothingForward(control, mc.cores, trace)
          self$smoothingBackward(control, mc.cores, trace)
        }
    },
    smoothingBackward = function(control, mc.cores, trace) {
      if (trace) cat("   Going backward ")
      for (i in rev(self$vBlocks[-1])) {
        if (trace) cat('+')
        cl0 <- factor(self$models[[i]]$fittedSBM$memberships)
        if (nlevels(cl0) == i) {
          candidates <- mclapply(combn(i, 2, simplify = FALSE), function(couple) {
            cl_fusion <- cl0
            levels(cl_fusion)[which(levels(cl_fusion) == paste(couple[1]))] <- paste(couple[2])
            levels(cl_fusion) <- as.character(1:(i - 1))
            model <- missingSBM_fit$new(sampledNet, i - 1, self$sampling, cl_fusion, self$covarMatrix, self$covarArray)
            model$doVEM(control)
            model
          }, mc.cores = mc.cores)
          vICLs <- sapply(candidates, function(candidate) candidate$vICL)
          best_one <- candidates[[which.min(vICLs)]]
          if (is.na(self$models[[i - 1]]$vICL)) {
            self$models[[i - 1]] <- best_one
          } else if (best_one$vICL < models[[i - 1]]$vICL) {
            self$models[[i - 1]] <- best_one
          }
        }
      }
    },
    smoothingForward = function(control, mc.cores, trace) {
      if (trace) cat("   Going forward ")
      for (i in self$vBlocks[-length(self$vBlocks)]) {
        if (trace) cat("+")
        cl_split <- factor(self$models[[i]]$fittedSBM$memberships)
        levels(cl_split) <- c(levels(cl_split), as.character(i + 1))
        if (nlevels(cl_split) - 1 == i) {
          candidates <- mclapply(1:i, function(j) {
            cl <- as.numeric(cl_split); indices <- which(cl == j)
            if (length(cl[indices]) > 10) { ## ????
              cut <- as.numeric(self$split_fn(sampledNet$adjMatrix[indices, indices],2))
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
          if (is.na(self$models[[i + 1]]$vICL)){
            self$models[[i + 1]] <- best_one
          } else if(best_one$vICL < models[[i + 1]]$vICL) {
            self$models[[i + 1]] <- best_one
          }
        }
      }
      if (trace) cat("\r                                                                                                        \r")
    }
  )
)
