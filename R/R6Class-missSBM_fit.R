#' An R6 class to represent an SBM fit with missing data
#'
#' @description The function [estimateMissSBM()] fits a collection of SBM for varying number of block.
#' Each fitted SBM is an instance of an R6 object with class [`missSBM_fit`], described here.
#'
#' Fields are accessed via active binding and cannot be changed by the user.
#'
#' This class comes with a set of R6 methods, some of them being useful for the user and exported
#' as S3 methods. See the documentation for  [show()], [print()], [fitted()], [predict()], [plot()].
#'
#' @examples
#' ## Sample 75% of dyads in  French political Blogosphere's network data
#' adjMatrix <- missSBM::frenchblog2007 %>%
#'   igraph::as_adjacency_matrix(sparse = FALSE) %>%
#'   missSBM::observeNetwork(sampling = "dyad", parameters = 0.75)
#' collection <- estimateMissSBM(adjMatrix, 3:5, sampling = "dyad")
#' my_missSBM_fit <- collection$bestModel
#' class(my_missSBM_fit)
#' plot(my_missSBM_fit, "imputed")
#'
#' @include R6Class-simpleSBM_fit.R
#' @include R6Class-networkSampling_fit.R
#' @export
missSBM_fit <-
  R6::R6Class(classname = "missSBM_fit",
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## fields for internal use (referring to mathematical notations)
  private = list(
    nu         = NULL, # imputed values (a sparse Matrix with entries only for imputed values, "dgCMatrix)
    sampling   = NULL, # fit of the current sampling model (object of class 'networkSampling_fit')
    SBM        = NULL, # fit of the current stochastic block model (object of class 'SBM_fit')
    optStatus  = NULL, # status of the optimization process

    ## kept so split()/merge() can build a sibling fit with one more/fewer block
    partlyObservedNet = NULL,
    netSampling       = NULL,
    useCov            = NULL,

    ## builds a new missSBM_fit from a candidate clustering, either replacing self's own fit
    ## (in_place = TRUE) or returning it as a new object -- shared by split() and merge()
    build_candidate = function(candidate_labels, in_place) {
      new_fit <- missSBM_fit$new(private$partlyObservedNet, private$netSampling, candidate_labels, private$useCov)
      if (in_place) {
        private$SBM      <- new_fit$fittedSBM
        private$sampling <- new_fit$fittedSampling
        private$nu       <- NULL
        return(invisible(self))
      }
      new_fit
    }
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description constructor for networkSampling
    #' @param partlyObservedNet An object with class [`partlyObservedNetwork`].
    #' @param netSampling The sampling design for the modelling of missing data: MAR designs ("dyad", "node") and MNAR designs ("double-standard", "block-dyad", "block-node" ,"degree")
    #' @param clusterInit Initial clustering: a vector with size \code{ncol(adjacencyMatrix)}, providing a user-defined clustering. The number of blocks is deduced from the number of levels in with \code{clusterInit}.
    #' @param useCov logical. If covariates are present in partlyObservedNet, should they be used for the inference or of the network sampling design, or just for the SBM inference? default is TRUE.
    initialize = function(partlyObservedNet, netSampling, clusterInit, useCov = TRUE) {

      ## Basic sanity checks
      stopifnot(netSampling %in% available_samplings)
      stopifnot(inherits(partlyObservedNet, "partlyObservedNetwork"))
      stopifnot(is.numeric(clusterInit))

      ## Initialize the SBM fit
      covariates <- array2list(partlyObservedNet$covarArray)
      if (!useCov) covariates <- list()

      if (length(covariates) == 0) {
        if (netSampling %in% c("double-standard", "block-node", "block-dyad")) {
          private$SBM <- SimpleSBM_fit_MNAR$new(partlyObservedNet, clusterInit)
        } else {
          private$SBM <- SimpleSBM_fit_noCov$new(partlyObservedNet, clusterInit)
        }
      } else {
        private$SBM <- SimpleSBM_fit_withCov$new(partlyObservedNet, clusterInit, covariates)
      }

      ## Initialize the sampling fit
      private$sampling <- switch(netSampling,
        "dyad"            = dyadSampling_fit$new(partlyObservedNet),
        "node"            = nodeSampling_fit$new(partlyObservedNet),
        "covar-dyad"      = covarDyadSampling_fit$new(partlyObservedNet),
        "covar-node"      = covarNodeSampling_fit$new(partlyObservedNet),
        "double-standard" = doubleStandardSampling_fit$new(partlyObservedNet),
        "block-node"      = blockNodeSampling_fit$new(partlyObservedNet, clustering_indicator(clusterInit)),
        "block-dyad"      = blockDyadSampling_fit$new(partlyObservedNet, clustering_indicator(clusterInit)),
        "degree"          = degreeSampling_fit$new(partlyObservedNet, clustering_indicator(clusterInit), private$SBM$connectParam$mean),
        "snowball"        = nodeSampling_fit$new(partlyObservedNet) # estimated sampling parameter not relevant
      )

      private$partlyObservedNet <- partlyObservedNet
      private$netSampling       <- netSampling
      private$useCov            <- useCov
    },
    #' @description a method to perform inference of the current missSBM fit with variational EM
    #' @param control a list of parameters controlling the variational EM algorithm. See details of function [estimateMissSBM()]
    doVEM = function(control = list(threshold = 1e-2, maxIter = 100, fixPointIter = 3, trace = TRUE)) {

      ## Starting the Variational EM algorithm
      if (control$trace) cat("\n Adjusting Variational EM for Stochastic Block Model\n")
      if (control$trace) cat("\n\tDyads are distributed according to a '",
                             ifelse(private$SBM$directed, "directed", "undirected"),"' SBM.\n", sep = "")
      if (control$trace) cat("\n\tImputation assumes a '", private$sampling$type,"' network-sampling process\n", sep = "")

      private$optStatus <- run_VEM(
        control    = control,
        init_stop  = FALSE,
        e_step     = function(fixPointIter) {
          for (i in seq.int(fixPointIter)) {
            # update the variational parameters for missing entries (a.k.a nu)
            private$nu <- private$sampling$update_imputation(private$SBM$imputation)
            # update the variational parameters for block memberships (a.k.a tau)
            private$SBM$update_blocks(log_lambda = private$sampling$log_lambda)
          }
        },
        m_step     = function() {
          # update the parameters of the SBM (a.k.a pi and theta)
          private$SBM$update_parameters(private$nu)
          # update the parameters of network sampling process (a.k.a psi)
          private$sampling$update_parameters(private$nu, private$SBM$probMemberships)
        },
        get_loglik = function() self$loglik,
        get_theta  = function() private$SBM$connectParam$mean,
        ## a lightweight snapshot of the fit (no cloning of the network data)
        snapshot   = function() list(SBM = private$SBM$get_state(), sampling = private$sampling$clone(), nu = private$nu),
        restore    = function(state) {
          private$SBM$set_state(state$SBM)
          private$sampling <- state$sampling
          private$nu       <- state$nu
        },
        reorder    = function() private$SBM$reorder()
      )

      invisible(private$optStatus)
    },
    #' @description clone of the current fit after splitting cluster \code{index} in two, via a
    #'   spectral bipartition of the sub-network it induces. Builds but does not fit the
    #'   candidate (see \code{candidates_split()}).
    #' @param index index (integer) of the cluster to split
    #' @param in_place replace \code{self}'s own fit (\code{TRUE}) or return a new object
    #'   (\code{FALSE}, the default)?
    #' @param base_net optional precomputed network to bipartition (as built internally at the
    #'   top of this method); lets \code{candidates_split()} avoid recomputing it once per
    #'   candidate.
    #' @return a new [`missSBM_fit`] with one more block, or \code{NULL} if \code{index} cannot
    #'   be split (its induced sub-network has zero variance)
    split = function(index, in_place = FALSE, base_net = NULL) {
      if (is.null(base_net)) {
        base_net <- self$imputedNetwork
        if (private$SBM$directed) base_net <- base_net %*% t(base_net)
      }
      cl0 <- private$SBM$memberships
      Q   <- private$SBM$nbBlocks

      A <- base_net[cl0 == index, cl0 == index]
      sdA <- sd(A)
      if (sdA == 0) return(NULL)

      A <- 1 / (1 + exp(-A / sdA))
      D <- 1 / sqrt(rowSums(abs(A)))
      L <- sweep(sweep(A, 1, D, "*"), 2, D, "*")
      Un <- eigen(L, symmetric = TRUE)$vectors[, 1:2]
      Un <- sweep(Un, 1, sqrt(rowSums(Un^2)), "/")
      bipartition <- kmeans_missSBM(Un, 2)

      split_labels <- bipartition
      split_labels[bipartition == 1] <- index
      split_labels[bipartition == 2] <- Q + 1
      candidate <- cl0
      candidate[candidate == index] <- split_labels
      candidate <- repair_empty_classes(candidate, Q + 1)

      private$build_candidate(candidate, in_place)
    },
    #' @description generate and cheaply trial-fit candidates obtained by splitting each
    #'   splittable cluster in two (see \code{split()}). A cluster is splittable if it has at
    #'   least 4 members and non-zero variance in its induced sub-network.
    #' @param control a list of VEM control parameters (see [estimateMissSBM()]); \code{maxIter}
    #'   is overridden by \code{trial_niter}
    #' @param trial_niter number of VEM iterations used for the trial fits. Default is 2.
    #' @return a list of trial-fitted [`missSBM_fit`] candidates (one per splittable cluster)
    candidates_split = function(control, trial_niter = 2) {
      base_net <- self$imputedNetwork
      if (private$SBM$directed) base_net <- base_net %*% t(base_net)
      cl0 <- private$SBM$memberships
      Q   <- private$SBM$nbBlocks

      cl_splitable <- (1:Q)[tabulate(cl0, nbins = Q) >= 4]
      sds <- sapply(cl_splitable, function(k_) sd(base_net[cl0 == k_, cl0 == k_]))
      cl_splitable <- cl_splitable[sds > 0]
      if (length(cl_splitable) == 0) return(list())

      control_fast <- control
      control_fast$maxIter <- trial_niter
      control_fast$trace   <- FALSE

      candidates <- future_lapply(cl_splitable, function(k_) {
        candidate <- self$split(k_, base_net = base_net)
        candidate$doVEM(control_fast)
        candidate
      }, future.seed = TRUE, future.scheduling = structure(TRUE, ordering = "random"))
      Filter(function(m) !is_degenerate(m), candidates)
    },
    #' @description clone of the current fit after merging clusters \code{indices[1]} and
    #'   \code{indices[2]} into one. Builds but does not fit the candidate (see
    #'   \code{candidates_merge()}).
    #' @param indices indices (couple of integers) of the clusters to merge
    #' @param in_place replace \code{self}'s own fit (\code{TRUE}) or return a new object
    #'   (\code{FALSE}, the default)?
    #' @return a new [`missSBM_fit`] with one fewer block
    merge = function(indices, in_place = FALSE) {
      Q  <- private$SBM$nbBlocks
      cl0 <- repair_empty_classes(private$SBM$memberships, Q)
      cl0 <- factor(cl0, 1:Q)

      cl_merged <- cl0
      levels(cl_merged)[sort(indices)] <- indices[1]
      levels(cl_merged) <- as.character(1:(nlevels(cl0) - 1))
      candidate <- as.integer(cl_merged)

      private$build_candidate(candidate, in_place)
    },
    #' @description generate and cheaply trial-fit candidates obtained by merging pairs of
    #'   clusters (see \code{merge()}). Beyond \code{max_candidates} pairs (quadratic in the
    #'   number of blocks), only the most similar-connectivity pairs are tried.
    #' @param control a list of VEM control parameters (see [estimateMissSBM()]); \code{maxIter}
    #'   is overridden by \code{trial_niter}
    #' @param max_candidates cap on the number of pairs tried. Default is 30.
    #' @param trial_niter number of VEM iterations used for the trial fits. Default is 2.
    #' @return a list of trial-fitted [`missSBM_fit`] candidates
    candidates_merge = function(control, max_candidates = 30, trial_niter = 2) {
      Q <- private$SBM$nbBlocks
      if (Q <= 1) return(list())

      pairs <- combn(Q, 2, simplify = FALSE)
      if (length(pairs) > max_candidates) {
        theta <- private$SBM$connectParam$mean
        score <- sapply(pairs, function(ij) {
          sqrt(sum((theta[ij[1], ] - theta[ij[2], ])^2) + sum((theta[, ij[1]] - theta[, ij[2]])^2))
        })
        pairs <- pairs[order(score)[1:max_candidates]]
      }

      control_fast <- control
      control_fast$maxIter <- trial_niter
      control_fast$trace   <- FALSE

      candidates <- future_lapply(pairs, function(couple) {
        candidate <- self$merge(couple)
        candidate$doVEM(control_fast)
        candidate
      }, future.seed = TRUE, future.scheduling = structure(TRUE, ordering = "random"))
      Filter(function(m) !is_degenerate(m), candidates)
    },
    #' @description show method for missSBM_fit
    show = function() {
      cat("missSBM-fit\n")
      cat("==================================================================\n")
      cat("Structure for storing an SBM fitted under missing data condition  \n")
      cat("==================================================================\n")
      cat("* Useful fields (first 2 special objects themselves with methods) \n")
      cat("  $fittedSBM (the adjusted stochastic block model)                \n")
      cat("  $fittedSampling (the estimated sampling process)                \n")
      cat("  $imputedNetwork (the adjacency matrix with imputed values)      \n")
      cat("  $monitoring, $ICL, $loglik, $vExpec, $penalty                  \n")
      cat("* S3 methods                                                      \n")
      cat("  plot, coef, fitted, predict, print                              \n")
    },
    #' @description User friendly print method
    print = function() { self$show() }
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## ACTIVE BINDING
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field fittedSBM the fitted SBM with class [`SimpleSBM_fit_noCov`], [`SimpleSBM_fit_withCov`] or
    #' [`SimpleSBM_fit_MNAR`] inheriting from class [`sbm::SimpleSBM_fit`]
    fittedSBM = function(value) {private$SBM},
    #' @field fittedSampling  the fitted sampling, inheriting from class [`networkSampling`] and corresponding fits
    fittedSampling = function(value) {private$sampling},
    #' @field imputedNetwork The network data as a matrix with NAs values imputed with the current model
    imputedNetwork = function(value) {
      res <- private$SBM$networkData + private$nu
      if (!private$SBM$directed) res <- res + t(res)
      res
    },
    #' @field monitoring a list carrying information about the optimization process
    monitoring     = function(value) {private$optStatus},
    #' @field entropyImputed the entropy of the distribution of the imputed dyads
    entropyImputed = function(value) {
      ## operate on the stored values only (private$nu is sparse, nonzero only at missing
      ## dyads): `xlogx(1 - private$nu)` would densify the whole N x N matrix and dispatch
      ## ifelse() through costly S4 machinery for every VEM iteration
      if (is.null(private$nu)) return(0)
      nu_x <- private$nu@x
      - sum(xlogx(nu_x) + xlogx(1 - nu_x))
    },
    #' @field entropy the entropy due to the distribution of the imputed dyads and of the clustering
    entropy = function(value) {private$SBM$entropy + self$entropyImputed},
    #' @field vExpec double: variational expectation of the complete log-likelihood
    vExpec  = function(value) {
      ## if(), not ifelse(): ifelse() evaluates both branches eagerly
      private$sampling$vExpec +
        if (private$sampling$type %in% c("block-dyad", "block-node", "double-standard"))
          private$SBM$vExpec else private$SBM$vExpec_corrected
    },
    #' @field penalty double, value of the penalty term in ICL
    penalty = function(value) {private$SBM$penalty + private$sampling$penalty},
    #' @field loglik double: approximation of the log-likelihood (variational lower bound) reached
    loglik  = function(value) {self$vExpec + self$entropy},
    #' @field ICL double: value of the integrated classification log-likelihood
    ICL    = function(value) {-2 * self$vExpec + self$penalty}
  )
)

## PUBLIC S3 METHODS FOR missSBMfit
## =========================================================================================

## Auxiliary function to check the given class of an object
is_missSBMfit <- function(Robject) {inherits(Robject, "missSBM_fit")}

#' Extract model fitted values from object  [`missSBM_fit`], return by [estimateMissSBM()]
#'
#' @name fitted.missSBM_fit
#'
#' @param object an R6 object with class [`missSBM_fit`]
#' @param ... additional parameters for S3 compatibility.
#'
#' @return A matrix of estimated probabilities of connection
#'
#' @importFrom stats fitted
#' @export
fitted.missSBM_fit <- function(object, ...) {
  stopifnot(is_missSBMfit(object))
  fitted(object$fittedSBM)
}

#' Prediction of a [`missSBM_fit`] (i.e. network with imputed missing dyads)
#'
#' @name predicted.missSBM_fit
#'
#' @param object an R6 object with class [`missSBM_fit`]
#' @param ... additional parameters for S3 compatibility.
#'
#' @return an adjacency matrix between pairs of nodes. Missing dyads are imputed with
#' their expected values, i.e. by there estimated probabilities of connection under the missing SBM.
#'
#' @importFrom stats predict
#' @export
predict.missSBM_fit <- function(object, ...) {
  stopifnot(is_missSBMfit(object))
  object$imputedNetwork
}

#' Summary method for a [`missSBM_fit`]
#'
#' @name summary.missSBM_fit
#'
#' @param object an R6 object with class [`missSBM_fit`]
#' @param ... additional parameters for S3 compatibility.
#'
#' @return a basic printing output
#'
#' @method summary missSBM_fit
#' @export
summary.missSBM_fit <- function(object, ...) {
  stopifnot(is_missSBMfit(object))
  object$show()
}

#' Visualization for an object [`missSBM_fit`]
#'
#' @description Plot function for the various fields of a [`missSBM_fit`]: the fitted
#' SBM (network or connectivity), and a plot monitoring the optimization.
#'
#' @return a ggplot object
#'
#' @name plot.missSBM_fit
#'
#' @param x an object with class [`missSBM_fit`]
#' @param type the type specifies the field to plot, either "imputed", "expected", "meso",  or "monitoring"
#' @param dimLabels : a list of two characters specifying the labels of the nodes. Default to \code{list(row= 'node',col = 'node')})
#' @param ... additional parameters for S3 compatibility. Not used
#' @export
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom sbm plotMyMatrix
plot.missSBM_fit <- function(x, type = c("imputed", "expected", "meso", "monitoring"), dimLabels = list(row= 'node',col = 'node'), ...) {
  stopifnot(is_missSBMfit(x))
  type <- match.arg(type)
  gg_obj <- switch(type,
    "expected"   = plotMyMatrix(x$fittedSBM$expectation, dimLabels, list(row = x$fittedSBM$memberships)),
    "meso"       = x$fittedSBM$plot("meso"),
    "imputed"    = plotMyMatrix(as.matrix(predict(x)), dimLabels,  list(row = x$fittedSBM$memberships)),
    "monitoring" = ggplot(x$monitoring, aes(x = .data$iteration, y = .data$elbo)) + geom_line() + theme_bw()
  )
  if (type != "meso") gg_obj ## return ggobject unless igraph is invoked
}

#' Extract model coefficients
#'
#' @description Extracts model coefficients from objects [`missSBM_fit`] returned by [estimateMissSBM()]
#'
#' @name coef.missSBM_fit
#'
#' @param object an R6 object with class [`missSBM_fit`]
#' @param type type of parameter that should be extracted. Either "mixture" (default), "connectivity", "covariates" or "sampling"
#' @param ... additional parameters for S3 compatibility. Not used
#' @return A vector or matrix of coefficients extracted from the missSBM_fit model.
#'
#' @export
coef.missSBM_fit <- function(object, type = c("mixture", "connectivity", "covariates", "sampling"), ...) {
  stopifnot(is_missSBMfit(object))
  switch(match.arg(type),
         mixture      = object$fittedSBM$blockProp,
         connectivity = object$fittedSBM$connectParam,
         covariates   = object$fittedSBM$covarParam,
         sampling     = object$fittedSampling$parameters)
}
