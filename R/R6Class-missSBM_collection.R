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

    # method for performing forward exploration of the ICL
    # a list of parameters controlling the variational EM algorithm. See details of function [estimateMissSBM()]
    explore_forward = function(control) {

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
            kmeans(Un, 2, nstart = 10)$cluster ## ClusterR::KMeans_rcpp() fdails on matrix(c(1,-1,-1,-1,-1,1,1,1), 4, 2)
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
    # method for performing backward exploration of the ICL
    # control a list of parameters controlling the variational EM algorithm. See details of function [`estimate`]
    explore_backward = function(control) {

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
    #' @param control a list of parameters controlling advanced features. Only 'cores', 'trace' and 'useCov' are relevant here. See [estimateMissSBM()] for details.
    initialize = function(partlyObservedNet, sampling, clusterInit, control) {

      if (control$trace) cat("\n")
      if (control$trace) cat("\n Adjusting Variational EM for Stochastic Block Model\n")
      if (control$trace) cat("\n\tImputation assumes a '", sampling,"' network-sampling process\n", sep = "")
      if (control$trace) cat("\n")

      stopifnot(inherits(partlyObservedNet, "partlyObservedNetwork"))
      private$partlyObservedNet <- partlyObservedNet
      if (control$trace) cat(" Initialization of", length(clusterInit), "model(s).", "\n")
      private$missSBM_fit <- lapply(clusterInit,
        function(cl0) {
          missSBM_fit$new(partlyObservedNet, sampling, cl0, control$useCov)
        }
      )
    },
    #' @description method to launch the estimation of the collection of models
    #' @param control a list of parameters controlling the variational EM algorithm. See details of function [estimateMissSBM()]
    estimate = function(control) {
      if (control$trace) cat(" Performing VEM inference\n")
      private$missSBM_fit <- mclapply(private$missSBM_fit, function(model) {
        if (control$trace) cat(" \tModel with", model$fittedSBM$nbBlocks,"blocks.\r")
        control$trace <- FALSE
        model$doVEM(control)
        model
      }, mc.cores = control$cores)
      invisible(self)
    },
    #' @description method for performing exploration of the ICL
    #' @param control a list of parameters controlling the exploration, similar to those found in the regular function [estimateMissSBM()]
    explore = function(control) {
      if (control$iterates > 0) {
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
