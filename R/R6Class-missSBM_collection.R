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
    forwardLimit      = NULL, # index into vBlocks explore_forward() won't grow past; NULL = no limit

    ## caps forwardLimit just before a persistent run of >= min_run degenerate models; adaptive,
    ## lifts the cap again once a later pass cures it
    flag_degenerate_tail = function(min_run) {
      degenerate <- self$degenerate
      run <- 0
      for (i in seq_along(degenerate)) {
        run <- if (degenerate[i]) run + 1 else 0
        if (run >= min_run) {
          new_limit <- i - min_run
          if (!identical(private$forwardLimit, new_limit)) {
            private$forwardLimit <- new_limit
            warning(
              "VEM keeps collapsing classes from nbBlocks = ", self$vBlocks[i - min_run + 1],
              " onward, even after repair() and exploration: the network does not appear to ",
              "support that many blocks. Forward (split) exploration will not grow past ",
              if (i > min_run) self$vBlocks[i - min_run] else "0",
              " blocks on further passes. See $occupiedBlocks/$degenerate.", call. = FALSE
            )
          }
          return(invisible(NULL))
        }
      }
      private$forwardLimit <- NULL # no persistent run (any more): lift a previous cap, if any
      invisible(NULL)
    },

    # for each q, ask the model at q for trial-fitted split candidates
    # (missSBM_fit$candidates_split()), fully refit the best one, and keep it in place of the
    # q+1 model only if it strictly improves the ICL.
    # control: a list of parameters controlling the variational EM algorithm. See details of
    # function [estimateMissSBM()]
    explore_forward = function(control) {

      trace <- control$trace; control$trace <- FALSE

      if (length(self$models) == 1) return(NULL)
      if (trace) cat("   Going forward ")
      last_k <- length(self$vBlocks) - 1
      if (!is.null(private$forwardLimit)) last_k <- min(last_k, private$forwardLimit)
      vBlocks <- self$vBlocks[seq_len(last_k)]
      for (k in seq_along(vBlocks)) {
        if (trace) cat("+")

        candidates <- private$missSBM_fit[[k]]$candidates_split(control)
        if (length(candidates) > 0) {
          best_one <- candidates[[which.min(sapply(candidates, function(m) m$ICL))]]
          best_one$doVEM(control)
          best_one$repair(control) # the full refit can itself collapse a component; try to recover it

          ## a gap in vBlocks (e.g. c(2, 3, 5)) means a split candidate's block count need not
          ## match its target slot -- reject it in that case rather than corrupt the slot
          expected_nbBlocks <- private$missSBM_fit[[k + 1]]$fittedSBM$nbBlocks
          if (!is_degenerate(best_one) && best_one$fittedSBM$nbBlocks == expected_nbBlocks &&
              best_one$ICL < private$missSBM_fit[[k + 1]]$ICL) {
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
          best_one$repair(control) # the full refit can itself collapse a component; try to recover it

          ## a gap in vBlocks (e.g. c(2, 3, 5)) means a merge candidate's block count need not
          ## match its target slot -- reject it in that case rather than corrupt the slot
          expected_nbBlocks <- private$missSBM_fit[[k - 1]]$fittedSBM$nbBlocks
          if (!is_degenerate(best_one) && best_one$fittedSBM$nbBlocks == expected_nbBlocks &&
              best_one$ICL < private$missSBM_fit[[k - 1]]$ICL) {
            private$missSBM_fit[[k - 1]] <- best_one
          }
        }
      }
      if (trace) cat("\r                                                                                                    \r")
    },
    plot_icl = function() {
      df <- data.frame(nBlock = self$vBlocks, ICL = self$ICL, collapsed = self$degenerate)
      ggplot(df, aes(x = nBlock, y = ICL)) +
        geom_line() + geom_point(aes(shape = collapsed), size = 2) + theme_bw() +
        scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4)) +
        labs(x = "#blocks", y = "Integrated Classification likelihood",
             shape = "collapsed classes") + ggtitle("Model Selection")
    },
    plot_elbo = function() {
      df <- data.frame(nBlock = self$vBlocks,
                        elbo = sapply(self$models, function(model) model$loglik),
                        collapsed = self$degenerate)
      ggplot(df, aes(x = nBlock, y = elbo)) +
        geom_line() + geom_point(aes(shape = collapsed), size = 2) + theme_bw() +
        scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4)) +
        labs(x = "#blocks", y = "Evidence (Varitional) Lower Bound",
             shape = "collapsed classes") + ggtitle("Model Selection")
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
        model$repair(control)
        model
      }, future.seed = TRUE, future.scheduling = structure(TRUE, ordering = "random"))
      invisible(self)
    },
    #' @description alternative to \code{estimate()}: fits each model in increasing order of
    #'   number of blocks, initializing \code{vBlocks[k]} by splitting (see [missSBM_fit]'s
    #'   \code{split()}/\code{candidates_split()}) the already-converged model at
    #'   \code{vBlocks[k-1]} instead of an independent, cold spectral clustering. Meant to reduce
    #'   VEM component collapse at higher numbers of blocks (see \code{$degenerate}), at the cost
    #'   of being sequential in the number of blocks (unlike \code{estimate()}, which fits every
    #'   model in parallel) -- can be slower in wall-clock time with many workers available.
    #'   Falls back to this slot's own cold-started clustering (built at construction, same as
    #'   \code{estimate()} would use) whenever nothing is splittable along the chain.
    #' @param control optional list of parameters overriding the collection's stored control.
    #'   Default \code{NULL} uses the stored control as-is.
    estimate_chain = function(control = NULL) {
      if (is.null(control)) control <- private$control
      if (control$trace) cat(" Performing chained VEM inference\n")
      ctrl <- control; ctrl$trace <- FALSE
      n <- length(private$missSBM_fit)

      ## smallest requested Q: fit its cold-started clustering directly
      if (control$trace) cat(" \tModel with", self$vBlocks[1], "blocks.\r")
      chained <- private$missSBM_fit[[1]]
      chained$doVEM(ctrl)
      chained$repair(ctrl)
      private$missSBM_fit[[1]] <- chained

      if (n > 1) for (k in 2:n) {
        if (control$trace) cat(" \tChaining to", self$vBlocks[k], "blocks.\r")
        gap <- self$vBlocks[k] - self$vBlocks[k - 1]
        stopifnot(gap >= 1) # vBlocks assumed strictly increasing, see estimateMissSBM()

        next_fit <- chained
        for (step in seq_len(gap)) {
          candidates <- next_fit$candidates_split(ctrl)
          if (length(candidates) == 0) {
            next_fit <- NULL # nothing (more) splittable along the chain: fall back below
            break
          }
          next_fit <- candidates[[which.min(sapply(candidates, function(m) m$ICL))]]
        }
        if (is.null(next_fit)) next_fit <- private$missSBM_fit[[k]] # this slot's cold start

        next_fit$doVEM(ctrl)
        next_fit$repair(ctrl)
        private$missSBM_fit[[k]] <- next_fit
        chained <- next_fit
      }

      if (control$trace) cat("\r                                                                                                    \r")
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
    #'   Uses the collection's stored control by default; \code{iterates} lets the caller override
    #'   it for this call only, without altering the stored control -- handy to alternate
    #'   \code{explore()}/\code{polish()} calls without having to reconstruct a full control list
    #'   each time. \code{iterates <= 0} is a no-op.
    #' @param control optional list of parameters overriding the collection's stored control.
    #'   Default \code{NULL} uses the stored control as-is.
    #' @param iterates optional integer overriding \code{control$iterates} for this call only.
    #' @param direction character ("forward", "backward", "both" or "none") controlling which
    #'   directions are searched. Default "both".
    explore = function(control = NULL, iterates = NULL, direction = "both") {
      if (is.null(control)) control <- private$control
      if (!is.null(iterates)) control$iterates <- iterates
      if (control$iterates > 0 && direction != "none") {
        if (control$trace) cat("\n Looking for better solutions\n")
        for (i in 1:control$iterates) {
          if (direction %in% c('forward' , 'both')) {
            if (control$trace) cat(" Pass",i)
            private$explore_forward(control)
          }
          if (direction %in% c('backward', 'both')) {
            if (control$trace) cat(" Pass",i)
            private$explore_backward(control)
          }
          if (isTRUE(control$stopOnDegenerate)) private$flag_degenerate_tail(control$maxConsecutiveDegenerate)
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
      degenerate <- self$degenerate
      if (any(degenerate)) {
        cat(" - Warning:", sum(degenerate), "model(s) have collapsed classes",
            "(see $occupiedBlocks/$degenerate)\n")
      }
      cat(" - Fields: $models, $ICL, $vBlocks, $occupiedBlocks, $degenerate, $bestModel, $optimizationStatus\n")
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
    #' @field bestModel the best model according to the ICL, restricted to models without
    #'   collapsed classes when at least one such model is available (see \code{$degenerate})
    bestModel = function(value) {
      icl <- self$ICL
      candidates <- which(!self$degenerate)
      if (length(candidates) == 0) candidates <- seq_along(icl)
      self$models[[candidates[which.min(icl[candidates])]]]
    },
    #' @field vBlocks a vector with the number of blocks
    vBlocks = function(value) {sapply(self$models, function(model) model$fittedSBM$nbBlocks)},
    #' @field occupiedBlocks a vector with the number of classes actually occupied in each model
    #'   (see [missSBM_fit]'s \code{occupiedBlocks})
    occupiedBlocks = function(value) {sapply(self$models, function(model) model$occupiedBlocks)},
    #' @field degenerate logical vector, \code{TRUE} for models with collapsed classes
    #'   (\code{occupiedBlocks < vBlocks}, see [missSBM_fit]'s \code{repair()})
    degenerate = function(value) {self$occupiedBlocks < self$vBlocks},
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
