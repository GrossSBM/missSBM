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
#' ## Uncomment to set parallel computing with future
#' ## future::plan("multicore", workers = 2)
#'
#' ## Sample 75% of dyads in  French political Blogosphere's network data
#' adjacencyMatrix <- missSBM::frenchblog2007 %>%
#'   igraph::delete.vertices(1:100) %>%
#'   igraph::as_adjacency_matrix() %>%
#'   missSBM::observeNetwork(sampling = "dyad", parameters = 0.75)
#' collection <- estimateMissSBM(adjacencyMatrix, 1:5, sampling = "dyad")
#' class(collection)
#'
#' @rdname missSBM_collection
#' @importFrom future.apply future_lapply
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
    control           = NULL, # default control list (set at construction by estimateMissSBM()),
                               # used by estimate()/polish()/explore() when not overridden per call

    # for each q, ask the model at q for trial-fitted split candidates
    # (missSBM_fit$candidates_split()), fully refit the best one, and keep it in place of the
    # q+1 model only if it strictly improves the ICL.
    # control: a list of parameters controlling the variational EM algorithm. See details of
    # function [estimateMissSBM()]
    explore_forward = function(control) {

      trace <- control$trace; control$trace <- FALSE

      if (length(self$models) == 1) return(NULL)
      if (trace) cat("   Going forward ")
      vBlocks <- self$vBlocks[-length(self$vBlocks)]
      for (k in seq_along(vBlocks)) {
        if (trace) cat("+")

        candidates <- private$missSBM_fit[[k]]$candidates_split(control)
        if (length(candidates) > 0) {
          best_one <- candidates[[which.min(sapply(candidates, function(m) m$ICL))]]
          best_one$doVEM(control)

          if (best_one$ICL < private$missSBM_fit[[k + 1]]$ICL) {
            private$missSBM_fit[[k + 1]] <- best_one
          }
        }
      }

      if (trace) cat("\r                                                                                                    \r")
    },
    # same as explore_forward(), but asking each model for merge candidates
    # (missSBM_fit$candidates_merge()) and propagating improvements to the q-1 neighbor.
    # control: a list of parameters controlling the variational EM algorithm. See details of
    # function [`estimate`]
    explore_backward = function(control) {

      trace <- control$trace; control$trace <- FALSE
      max_candidates <- control$maxMergeCandidates
      if (is.null(max_candidates)) max_candidates <- Inf

      if (length(self$models) == 1) return(NULL)
      if (trace) cat("   Going backward ")
      vBlocks <- self$vBlocks
      for (k in seq(from = length(vBlocks), to = 2, by = -1) ) {
        if (trace) cat("+")

        candidates <- private$missSBM_fit[[k]]$candidates_merge(control, max_candidates = max_candidates)
        if (length(candidates) > 0) {
          best_one <- candidates[[which.min(sapply(candidates, function(m) m$ICL))]]
          best_one$doVEM(control)

          if (best_one$ICL < private$missSBM_fit[[k - 1]]$ICL) {
            private$missSBM_fit[[k - 1]] <- best_one
          }
        }
      }
      if (trace) cat("\r                                                                                                    \r")
    },
    plot_icl = function() {
      qplot(self$vBlocks, self$ICL, geom = "line") + theme_bw() + geom_point() +
        labs(x = "#blocks", y = "Integrated Classification likelihood") + ggtitle("Model Selection")
    },
    plot_elbo = function() {
      elbo <- sapply(self$models, function(model) model$loglik)
      qplot(self$vBlocks, elbo, geom = "line") + theme_bw() + geom_point() +
        labs(x = "#blocks", y = "Evidence (Varitional) Lower Bound") + ggtitle("Model Selection")
    },
    plot_monitoring = function() {
      monitoring <- self$optimizationStatus
      monitoring$nBlock <- as.factor(monitoring$nBlock)
      ggplot(monitoring, aes(x = cumsum(iteration), y = elbo, color = nBlock)) + geom_line() + geom_point() +
        theme_bw() + labs(x = "# cumulated V-EM iterations", y = "Evidence (variational) Lower Bound") + ggtitle("Optimization monitoring")
    }
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description constructor for networkSampling
    #' @param partlyObservedNet An object with class [`partlyObservedNetwork`].
    #' @param sampling The sampling design for the modelling of missing data: MAR designs ("dyad", "node") and MNAR designs ("double-standard", "block-dyad", "block-node" ,"degree")
    #' @param clusterInit Initial clustering: a list of vectors, each with size \code{ncol(adjacencyMatrix)}.
    #' @param control a list of parameters controlling advanced features. Only 'trace' and 'useCov' are relevant here. See [estimateMissSBM()] for details.
    initialize = function(partlyObservedNet, sampling, clusterInit, control) {

      if (control$trace) cat("\n")
      if (control$trace) cat("\n Adjusting Variational EM for Stochastic Block Model\n")
      if (control$trace) cat("\n\tImputation assumes a '", sampling,"' network-sampling process\n", sep = "")
      if (control$trace) cat("\n")

      stopifnot(inherits(partlyObservedNet, "partlyObservedNetwork"))
      private$partlyObservedNet <- partlyObservedNet
      private$control           <- control
      if (control$trace) cat(" Initialization of", length(clusterInit), "model(s).", "\n")
      private$missSBM_fit <- lapply(clusterInit,
        function(cl0) {
          missSBM_fit$new(partlyObservedNet, sampling, cl0, control$useCov)
        }
      )
    },
    #' @description method to launch the estimation of the collection of models
    #' @param control optional list of parameters overriding the collection's stored control
    #'   (set at construction by [estimateMissSBM()], see its details for the full list).
    #'   Default \code{NULL} uses the stored control as-is.
    estimate = function(control = NULL) {
      if (is.null(control)) control <- private$control
      if (control$trace) cat(" Performing VEM inference\n")
      private$missSBM_fit <- future_lapply(private$missSBM_fit, function(model) {
        if (control$trace) cat(" \tModel with", model$fittedSBM$nbBlocks,"blocks.\r")
        control$trace <- FALSE
        model$doVEM(control)
        model
      }, future.seed = TRUE, future.scheduling = structure(TRUE, ordering = "random"))
      invisible(self)
    },
    #' @description method to node-swap-polish every model in the collection (see
    #'   [missSBM_fit]'s \code{polish()}); fixes individually misclassified nodes at each
    #'   model's own number of blocks, unlike \code{explore()} which searches across blocks.
    #' @param control optional list of parameters overriding the collection's stored control.
    #'   Default \code{NULL} uses the stored control as-is.
    polish = function(control = NULL) {
      if (is.null(control)) control <- private$control
      if (control$trace) cat(" Polishing (node-swap)\n")
      private$missSBM_fit <- future_lapply(private$missSBM_fit, function(model) {
        control$trace <- FALSE
        model$polish(control)
        model
      }, future.seed = TRUE, future.scheduling = structure(TRUE, ordering = "random"))
      invisible(self)
    },
    #' @description method for performing exploration of the ICL (split/merge search across
    #'   numbers of blocks, see [missSBM_fit]'s \code{candidates_split()}/\code{candidates_merge()}).
    #'   Uses the collection's stored control by default; \code{iterates} and \code{direction}
    #'   let the caller override just those two aspects for this call, without altering the
    #'   stored control -- handy to alternate \code{explore()}/\code{polish()} calls without
    #'   having to reconstruct a full control list each time.
    #' @param control optional list of parameters overriding the collection's stored control.
    #'   Default \code{NULL} uses the stored control as-is.
    #' @param iterates optional integer overriding \code{control$iterates} for this call only.
    #' @param direction optional character ("forward", "backward", "both" or "none") overriding
    #'   \code{control$exploration} for this call only.
    explore = function(control = NULL, iterates = NULL, direction = NULL) {
      if (is.null(control)) control <- private$control
      if (!is.null(iterates))  control$iterates    <- iterates
      if (!is.null(direction)) control$exploration <- direction
      if (control$iterates > 0 && control$exploration != "none") {
        if (control$trace) cat("\n Looking for better solutions\n")
        for (i in 1:control$iterates) {
          if (control$exploration %in% c('forward' , 'both')) {
            if (control$trace) cat(" Pass",i)
            private$explore_forward(control)
          }
          if (control$exploration %in% c('backward', 'both')) {
            if (control$trace) cat(" Pass",i)
            private$explore_backward(control)
          }
        }
      }
    },
    #' @description plot method for missSBM_collection
    #' @param type the type specifies the field to plot, either "icl", "elbo" or "monitoring". Default is "icl"
    plot = function(type = c("icl", "elbo", "monitoring")) {
      gg_obj <- switch(match.arg(type),
          "icl"        = private$plot_icl(),
          "elbo"       = private$plot_elbo(),
          "monitoring" = private$plot_monitoring()
      )
      gg_obj
    },
    #' @description show method for missSBM_collection
    show = function() {
      cat("--------------------------------------------------------\n")
      cat("COLLECTION OF", length(self$vBlocks), "SBM fits          \n")
      cat("========================================================\n")
      cat(" - Number of blocks considered: from ", min(self$vBlocks), " to ", max(self$vBlocks),"\n", sep = "")
      cat(" - Best model (smaller ICL): ", self$bestModel$fittedSBM$nbBlocks, "\n", sep = "")
      cat(" - Fields: $models, $ICL, $vBlocks, $bestModel, $optimizationStatus\n")
      cat(" - Method: $estimate(), $polish(), $explore(), $plot() \n")
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
    #' @field optimizationSettings the control list used by estimate()/polish()/explore() when
    #'   not overridden per call (set at construction by [estimateMissSBM()])
    optimizationSettings = function(value) {private$control},
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
