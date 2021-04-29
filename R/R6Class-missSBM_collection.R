#' An R6 class to represent a collection of SBM fits with missing data
#'
#' @description The function [estimateMissSBM()] fits a collection of SBM with missing data for
#' a varying number of block. These models with class [`missSBM_fit`]  are stored in an instance
#' of an object with class [`missSBM_collection`], described here.
#'
#' Fields are accessed via active binding and cannot be changed by the user.
#'
#' This class comes with a set of R6 methods, some of them being useful for the user and exported
#' as S3 methods. See the documentation for [show()], [print()] and [smooth()], the latter being
#' used to smooth the ICL on a collection of model, as post-treatment.
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

      trace <- control$trace > 0; control$trace <- FALSE
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
            n <- ncol(A)
            cl <- rep(1L, n)
            unconnected <- which(rowSums(abs(A)) == 0)
            connected   <- setdiff(1:n, unconnected)
            A <- A[connected,connected]
            D <- 1/sqrt(rowSums(abs(A)))
            L <- sweep(sweep(A, 1, D, "*"), 2, D, "*")
            Un <- base::svd(L, nu = 2, nv = 0)$u
            Un <- sweep(Un, 1, sqrt(rowSums(Un^2)), "/")
            Un[is.nan(Un)] <- 0
            cl_ <- ClusterR::KMeans_rcpp(Un, 2, num_init = 10)$clusters
            cl[connected] <- cl_
            cl[unconnected] <- which.min(rowsum(D, cl_))
            cl
          }, mc.cores = control$cores)

          ## build list of candidate clustering after splits
          cl_candidates <- mclapply(cl_splitable, function(k_)  {
            split <- cl_split[[k_]]
            split[cl_split[[k_]] == 1] <- k_
            split[cl_split[[k_]] == 2] <- vBlocks[k] + 1
            candidate <- cl0
            candidate[candidate == k_] <- split
            as.numeric(as.factor(candidate)) # relabeling to start from 1
          }, mc.cores = control$cores)

          loglik_candidates <- mclapply(cl_candidates, function(cl_) {
            model <- missSBM_fit$new(private$partlyObservedNet, sampling, as.integer(cl_), useCov)
            model$doVEM(control_fast)
            model$loglik
          }, mc.cores = control$cores) %>% unlist()

          best_one <-  missSBM_fit$new(private$partlyObservedNet, sampling, cl_candidates[[which.max(loglik_candidates)]], useCov)
          best_one$doVEM(control)

          if (best_one$loglik > private$missSBM_fit[[k + 1]]$loglik) {
            private$missSBM_fit[[k + 1]] <- best_one
          }

        }
      }
      if (trace) cat("\r                                                                                                    \r")
    },
    # method for performing backward smoothing of the ICL
    # control a list of parameters controlling the variational EM algorithm. See details of function [`estimate`]
    smoothing_backward = function(control) {

      trace <- control$trace > 0; control$trace <- FALSE
      control_fast <- control
      control_fast$maxIter <- 2

      sampling    <- private$missSBM_fit[[1]]$fittedSampling$type
      useCov      <- private$missSBM_fit[[1]]$fittedSBM$nbCovariates > 0

      if (length(self$models) == 1) return(NULL)
      if (trace) cat("   Going backward ")
      vBlocks <- self$vBlocks
      for (k in seq(from = length(vBlocks), to = 2, by = -1) ) {
        if (vBlocks[k] >= 2 ) {
          if (trace) cat("+")
          cl0 <- factor(private$missSBM_fit[[k]]$fittedSBM$memberships)
          ## build list of candidate clustering after merge
          cl_candidates <- mclapply(combn(vBlocks[k], 2, simplify = FALSE), function(couple) {
            cl_merged <- cl0
            levels(cl_merged)[which(levels(cl_merged) == paste(couple[1]))] <- paste(couple[2])
            levels(cl_merged) <- as.character(1:(vBlocks[k] - 1))
            as.integer(cl_merged)
          }, mc.cores = control$cores)

          loglik_candidates <- mclapply(cl_candidates, function(cl_) {
            model <- missSBM_fit$new(private$partlyObservedNet, sampling, as.integer(cl_), useCov)
            model$doVEM(control_fast)
            model$loglik
          }, mc.cores = control$cores) %>% unlist()

          best_one <-  missSBM_fit$new(private$partlyObservedNet, sampling, cl_candidates[[which.max(loglik_candidates)]], useCov)
          best_one$doVEM(control)

          if (best_one$loglik > private$missSBM_fit[[k - 1]]$loglik) {
            private$missSBM_fit[[k - 1]] <- best_one
          }
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
      trace_main <- control$trace > 0
      control$trace <- ifelse (control$trace > 1, TRUE, FALSE)
      if (trace_main) cat("\n")
      private$missSBM_fit <- mclapply(private$missSBM_fit, function(model) {
        if (trace_main) cat(" Performing VEM inference for model with", model$fittedSBM$nbBlocks,"blocks.\r")
        model$doVEM(control)
        model
      }, mc.cores = control$cores)
      invisible(self)
    },
    #' @description method for performing smoothing of the ICL
    #' @param control a list of parameters controlling the smoothing. See details of regular function [smooth()]
    smooth = function(control) {
      if (control$trace > 0) control$trace <- TRUE else control$trace <- FALSE
      if (control$trace) cat("\n Smoothing ICL\n")
      for (i in 1:control$iterates) {
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

#' Smooth the path ICL in a collection of missSBM_fit models
#'
#' Apply a split and/or merge strategy of the clustering in a path of models in a collection
#' of SBM ordered by number of block. The goal is to find better initialization. This results
#' in a "smoothing" of the ICL, that should be close to concave.
#'
#' @param Robject an object with class missSBM_collection, i.e. an output from [estimateMissSBM()]
#' @param type character indicating what kind of ICL smoothing should be use among "forward", "backward" or "both". Default is "both".
#' @param control a list controlling the variational EM algorithm. See details.
#'
#' @details The list of parameters \code{control} controls the optimization process and the variational EM algorithm, with the following entries
#'  \itemize{
#'  \item{"iterates": }{integer for the number of iterations of smoothing. Default is 1.}
#'  \item{"threshold": }{V-EM algorithm stops stop when an optimization step changes the objective function or the parameters
#'         by less than threshold. Default is 1e-3.}
#'  \item{"maxIter": }{V-EM algorithm stops when the number of iteration exceeds maxIter.
#'        Default is 100 with no covariate, 50 otherwise.}
#'  \item{"fixPointIter": }{number of fix-point iterations in the V-E step.
#'        Default is 5 with no covariate, 2 otherwise.}
#'  \item{"cores": }{integer for number of cores used. Default is 1.}
#'  \item{"trace": }{integer for verbosity. Useless when \code{cores} > 1}
#' }
#' @return An invisible missSBM_collection, in which the ICL has been smoothed
#' @export
smooth <- function(Robject, type = c("both", "forward", "backward"), control = list()) {

  stopifnot(inherits(Robject, "missSBM_collection"))

  ## defaut control parameter for VEM, overwritten by user specification
  ctrl <- list(threshold = 1e-2, maxIter = 50, fixPointIter = 3, cores = 1, trace = 1, iterates = 1)
  ctrl[names(control)] <- control
  ctrl$smoothing <- match.arg(type)
  if(Sys.info()['sysname'] == "Windows") ctrl$cores <- 1

  ## Run the smoothing
  Robject$smooth(ctrl)

  invisible(Robject)
}
