#' An R6 class to represent a collection of SBM fits with missing data
#'
#' @description The function [estimateMissSBM()] fits a collection of SBM with missing data for
#' a varying number of block. These models with class [`missSBM_fit`]  are stored in an instance
#' of an object with class [`missSBM_collection`], described here.
#'
#' Fields are accessed via active binding and cannot be changed by the user.
#'
#' This class comes with a set of R6 methods, some of them being useful for the user and exported
#' as S3 methods. See the documentation for [show()] and [print()]
#'
#' @examples
#' ## Sample 75% of dyads in  French political Blogosphere's network data
#' adjacencyMatrix <- missSBM::frenchblog2007 %>%
#'   igraph::as_adj (sparse = FALSE) %>%
#'   missSBM::observeNetwork(sampling = "dyad", parameters = 0.75)
#' collection <- estimateMissSBM(adjacencyMatrix, 3:5, sampling = "dyad")
#' class(collection)
#'
#'
#' @rdname missSBM_collection
#' @importFrom parallel mclapply
#' @include R6Class-missSBM_fit.R
#' @export
missSBM_collection <-
  R6::R6Class(classname = "missSBM_collection",
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## fields for internal use (referring to mathematical notations)
  private = list(
    partlyObservedNet = NULL, # network data with convenient encoding (object of class 'partlyObservedNetwork')
    missSBM_fit       = NULL, # a list of models

    # method for performing forward smoothing of the ICL
    # a list of parameters controlling the variational EM algorithm. See details of function [estimateMissSBM()]
    smoothing_forward = function(control) {

      trace <- control$trace; control$trace <- FALSE
      control_fast <- control
      control_fast$maxIter <- 2
      sampling <- private$missSBM_fit[[1]]$fittedSampling$type
      useCov   <- private$missSBM_fit[[1]]$fittedSBM$nbCovariates > 0

      if (length(self$models) == 1) return(NULL)
      if (trace) cat("   Going forward ")
      vBlocks <- self$vBlocks[-length(self$vBlocks)]

      for (k in seq_along(vBlocks)) {
        if (trace) cat("+")

        ## current  imputed network
        base_net <- private$missSBM_fit[[k]]$imputedNetwork
        if (private$missSBM_fit[[k]]$fittedSBM$directed) base_net <- base_net %*% t(base_net)
        ## current clustering
        cl0 <- private$missSBM_fit[[k]]$fittedSBM$memberships

        cl_splitable <- (1:vBlocks[k])[tabulate(cl0, nbins = vBlocks[k]) >= 4]
        sd <- sapply(cl_splitable, function(k_) sd(base_net[cl0 == k_, cl0 == k_]))
        cl_splitable <- cl_splitable[sd > 0]

        if (length(cl_splitable) > 0) {
          cl_split <- vector("list", vBlocks[k])
          cl_split[cl_splitable] <- mclapply(cl_splitable, function(k_) {
            A <- base_net[cl0 == k_, cl0 == k_]
            A <- 1/(1+exp(-A/sd(A)))
            D <- 1/sqrt(rowSums(abs(A)))
            L <- sweep(sweep(A, 1, D, "*"), 2, D, "*")
            Un <- eigen(L, symmetric = TRUE)$vectors[, 1:2]
            Un <- sweep(Un, 1, sqrt(rowSums(Un^2)), "/")
            ClusterR::KMeans_rcpp(Un, 2, num_init = 10)$clusters
          }, mc.cores = control$cores)

          ## build list of candidate clustering after splits
          cl_candidates <- lapply(cl_splitable, function(k_)  {
            split <- cl_split[[k_]]
            split[cl_split[[k_]] == 1] <- k_
            split[cl_split[[k_]] == 2] <- vBlocks[k] + 1
            candidate <- cl0
            candidate[candidate == k_] <- split
            ## in case of empty classes, add randomly one guy in thoses classes
            candidate <- factor(candidate, levels = 1:(vBlocks[k] + 1))
            absent <- which(tabulate(candidate) == 0)
            swap <- base::sample(1:length(candidate), length(absent))
            candidate[swap] <- absent
            as.numeric(candidate) # relabeling to start from 1
          })

          icl_candidates <- mclapply(cl_candidates, function(cl_) {
            model <- missSBM_fit$new(private$partlyObservedNet, sampling, as.integer(cl_), useCov)
            model$doVEM(control_fast)
            model$ICL
          }, mc.cores = control$cores) %>% unlist()

          best_one <-  missSBM_fit$new(private$partlyObservedNet, sampling, cl_candidates[[which.min(icl_candidates)]], useCov)
          best_one$doVEM(control)

          if (best_one$ICL < private$missSBM_fit[[k + 1]]$ICL) {
            private$missSBM_fit[[k + 1]] <- best_one
          }

        }
      }

      if (trace) cat("\r                                                                                                    \r")
    },
    # method for performing backward smoothing of the ICL
    # control a list of parameters controlling the variational EM algorithm. See details of function [`estimate`]
    smoothing_backward = function(control) {

      trace <- control$trace; control$trace <- FALSE
      control_fast <- control
      control_fast$maxIter <- 2
      sampling    <- private$missSBM_fit[[1]]$fittedSampling$type
      useCov      <- private$missSBM_fit[[1]]$fittedSBM$nbCovariates > 0

      if (length(self$models) == 1) return(NULL)
      if (trace) cat("   Going backward ")
      vBlocks <- self$vBlocks
      for (k in seq(from = length(vBlocks), to = 2, by = -1) ) {
        if (trace) cat("+")
        cl0 <- factor(private$missSBM_fit[[k]]$fittedSBM$memberships, 1:vBlocks[k])
        absent <- which(tabulate(cl0) == 0)
        swap <- base::sample(1:length(cl0), length(absent))
        cl0[swap] <- absent

        ## build list of candidate clustering after merge
        cl_candidates <- lapply(combn(vBlocks[k], 2, simplify = FALSE), function(couple) {
          cl_merged <- cl0
          levels(cl_merged)[couple] <- couple[1]
          levels(cl_merged) <- as.character(1:(nlevels(cl0) - 1))
          as.integer(cl_merged)
        })

        icl_candidates <- mclapply(cl_candidates, function(cl_) {
          model <- missSBM_fit$new(private$partlyObservedNet, sampling, as.integer(cl_), useCov)
          model$doVEM(control_fast)
          model$ICL
        }, mc.cores = control$cores) %>% unlist()

        best_one <-  missSBM_fit$new(private$partlyObservedNet, sampling, cl_candidates[[which.min(icl_candidates)]], useCov)
        best_one$doVEM(control)

        if (best_one$ICL < private$missSBM_fit[[k - 1]]$ICL) {
          private$missSBM_fit[[k - 1]] <- best_one
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
    #' @param partlyObservedNet An object with class [`partlyObservedNetwork`].
    #' @param sampling The sampling design for the modelling of missing data: MAR designs ("dyad", "node") and NMAR designs ("double-standard", "block-dyad", "block-node" ,"degree")
    #' @param clusterInit Initial clustering: a list of vectors, each with size \code{ncol(adjacencyMatrix)}.
    #' @param cores integer for number of cores used. Default is 1.
    #' @param trace integer for verbosity (0, 1, 2). Default is 1. Useless when \code{cores} > 1
    #' @param useCov logical. If covariates are present in partlyObservedNet, should they be used for the inference or of the network sampling design, or just for the SBM inference? default is TRUE.
    initialize = function(partlyObservedNet, sampling, clusterInit, cores, trace, useCov) {

      if (trace) cat("\n")
      if (trace) cat("\n Adjusting Variational EM for Stochastic Block Model\n")
      if (trace) cat("\n\tImputation assumes a '", sampling,"' network-sampling process\n", sep = "")
      if (trace) cat("\n")

      stopifnot(inherits(partlyObservedNet, "partlyObservedNetwork"))
      private$partlyObservedNet <- partlyObservedNet
      private$missSBM_fit <-lapply(clusterInit,
        function(cl0) {
          if (trace) cat(" Initialization of model with", length(unique(cl0)), "blocks.", "\r")
          missSBM_fit$new(partlyObservedNet, sampling, cl0, useCov)
        }
      )
    },
    #' @description method to launch the estimation of the collection of models
    #' @param control a list of parameters controlling the variational EM algorithm. See details of function [estimateMissSBM()]
    estimate = function(control) {
      if (control$trace) cat("\n")
      private$missSBM_fit <- mclapply(private$missSBM_fit, function(model) {
        if (control$trace) cat(" Performing VEM inference for model with", model$fittedSBM$nbBlocks,"blocks.\r")
        model$doVEM(control)
        model
      }, mc.cores = control$cores)
      invisible(self)
    },
    #' @description method for performing smoothing of the ICL
    #' @param control a list of parameters controlling the smoothing, similar to those found in the regular function [estimateMissSBM()]
    smooth = function(control) {
      if (control$trace) cat("\n Smoothing ICL\n")
      prop_swap <- control$prop_swap
      if (length(prop_swap) == 1) prop_swap <- rep(prop_swap, control$iterates)
      for (i in 1:control$iterates) {
        control$prop_swap <- prop_swap[i]
        if (control$smoothing %in% c('forward' , 'both')) private$smoothing_forward(control)
        if (control$smoothing %in% c('backward', 'both')) private$smoothing_backward(control)
      }
    },
    #' @description show method for missSBM_collection
    show = function() {
      cat("--------------------------------------------------------\n")
      cat("COLLECTION OF", length(self$vBlocks), "SBM fits          \n")
      cat("========================================================\n")
      cat(" - Number of blocks considers: from ", min(self$vBlocks), " to ", max(self$vBlocks),"\n", sep = "")
      cat(" - Best model (smaller ICL): ", self$bestModel$fittedSBM$nbBlocks, "\n", sep = "")
      cat(" - Fields: $models, $ICL, $vBlocks, $bestModel, $optimizationStatus\n")
    },
    #' @description User friendly print method
    print = function() { self$show() }
  ),
  active = list(
    #' @field models a list of models
    models = function(value) (private$missSBM_fit),
    #' @field ICL the vector of Integrated Classification Criterion (ICL) associated to the models in the collection (the smaller, the better)
    ICL = function(value) {setNames(sapply(self$models, function(model) model$ICL), names(self$vBlocks))},
    #' @field bestModel the best model according to the ICL
    bestModel = function(value) {self$models[[which.min(self$ICL)]]},
    #' @field vBlocks a vector with the number of blocks
    vBlocks = function(value) {sapply(self$models, function(model) model$fittedSBM$nbBlocks)},
    #' @field optimizationStatus a data.frame summarizing the optimization process for all models
    optimizationStatus = function(value) {
      Reduce("rbind",
             lapply(self$models, function(model) {
               res <- model$monitoring
               res$nBlock <- model$fittedSBM$nbBlocks
               res
             })
      )
    }
  )
)
