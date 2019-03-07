#' Inference of Stochastic Block Model from sampled data
#'
#' Perform variational inference of Stochastic Block Model from sampled adjacency matrix
#'
#' @param adjacencyMatrix The adjacency matrix of the network
#' @param vBlocks The vector of number of blocks considered in the collection
#' @param sampling The sampling design for missing data modeling : "dyad", "double_standard", "node", "snowball", "degree", "block"
#' by default "undirected" is choosen
#' @param clusterInit character in "hierarchical" or "spectral" for initialization
#' @param trace logical, control the verbosity. Default to \code{TRUE}.
#' @param mc.cores integer, the number of cores to use when multiply model are fitted
#' @param smoothing character indicating what kind of ICL smoothing should be use among "none", "forward", "backward" or "both"
#' @param iter_both integer for the number of iteration in case of foward-backward (aka both) smoothing
#' @param control_VEM a list controlling the variational EM algorithm. See details.
#' @param Robject an object with class \code{missSBMcollection}
#' @return Returns an S3 object with class \code{missSBMcollection}, which is a list with all models estimated for all Q in vBlocks. \code{missSBMcollection} owns a couple of S3 methods: \code{is.missSBMcollection} to test the class of the object, a method \code{ICL} to extract the values of the Integrated Classification Criteria for each model, a method \code{getBestModel} which extract from the list the best model (and object of class \code{missSBM-fit}) according to the ICL, and a method \code{optimizationStatus} to monitor the objective function a convergence of the VEM algorithm.
#' @references [1] Tabouy, P. Barbillon, J. Chiquet. Variationnal inference of Stochastic Block Model from sampled data (2017). arXiv:1707.04141.
#' @seealso \code{\link{sampleNetwork}} and \code{\link{simulateSBM}} and \code{\link{missingSBM_fit}}.
#' @examples
#' ### A SBM model : ###
#' N <- 300
#' Q <- 3
#' alpha <- rep(1,Q)/Q     # mixture parameter
#' pi <- diag(.45,Q) + .05 # connectivity matrix
#' directed <- FALSE
#' ## Simulation of an Bernoulli non-directed SBM
#' mySBM <- simulateSBM(N, alpha, pi, directed)
#'
#' ## Sampling of the data : ##
#' samplingParameters <- .5 # the sampling rate
#' sampling <- "dyad"       # the sampling design
#' sampledNet <-
#'    sampleNetwork(
#'      mySBM$adjMatrix,
#'      sampling,
#'      samplingParameters
#'    )
#'
#' ## Inference :
#' vBlocks <- 1:5 # number of classes
#' sbm <-
#'    missSBM(
#'       sampledNet$adjacencyMatrix,
#'       vBlocks,
#'       sampling
#'    )
#' @import R6 parallel
#' @include smoother_SBM.R
#' @export
missSBM <- function(
  adjacencyMatrix,
  vBlocks,
  sampling,
  clusterInit = "hierarchical",
  trace     = TRUE,
  smoothing = c("none", "forward", "backward", "both"),
  mc.cores = 1,
  iter_both = 1,
  control_VEM = list()) {

  ## some sanity checks
  try(
    !all.equal(unique(as.numeric(adjacencyMatrix[!is.na(adjacencyMatrix)])), c(0,1)),
    stop("Only binary graphs are supported.")
  )

  sampledNet <- sampledNetwork$new(adjacencyMatrix)
  if (trace) cat("\n")
  if (trace) cat("\n Adjusting Variational EM for Stochastic Block Model\n")
  if (trace) cat("\n\tImputation assumes a '", sampling,"' network-sampling process\n", sep = "")
  if (trace) cat("\n")
  models <- mclapply(vBlocks,
    function(nBlocks) {
    if (trace) cat(" Initialization of model with", nBlocks,"blocks.", "\r")
      if (is.list(clusterInit)) {
        missingSBM_fit$new(sampledNet, nBlocks, sampling, clusterInit[[nBlocks]])
      } else {
        missingSBM_fit$new(sampledNet, nBlocks, sampling, clusterInit)
      }
    }, mc.cores = mc.cores
  )

  ## defaut control parameter for VEM, overwritten by user specification
  control <- list(threshold = 1e-4, maxIter = 200, fixPointIter = 5, trace = FALSE)
  control[names(control_VEM)] <- control_VEM
  cat("\n")
  mclapply(models,
    function(model) {
      if (trace) cat(" Performing VEM inference for model with", model$fittedSBM$nBlocks,"blocks.\r")
      model$doVEM(control)
    }, mc.cores = mc.cores
  )

  smoothing <- match.arg(smoothing)
  if (smoothing != "none") {
    if (trace) cat("\n Smoothing ICL\n")
    smoothing_fn <- switch(smoothing,
                           "forward"  = smoothingForward ,
                           "backward" = smoothingBackward,
                           "both"     = smoothingForBackWard
    )
    if (!is.character(clusterInit)) {
      split_fn <- init_hierarchical
    } else {
      split_fn <- switch(clusterInit,
                         "spectral" = init_spectral,
                         "hierarchical" = init_hierarchical,
                         init_hierarchical)
    }
    control$trace <- FALSE # forcing no trace while smoothing
    models <- smoothing_fn(models, vBlocks, sampledNet, sampling, split_fn, mc.cores, iter_both, control)

  }

  structure(setNames(models, vBlocks), class = "missSBMcollection")
}

