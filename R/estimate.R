#' Inference of an SBM with missing data
#'
#' Perform variational inference of a Stochastic Block Model from a sampled adjacency matrix
#'
#' @param adjacencyMatrix The adjacency matrix of the network
#' @param vBlocks The vector of number of blocks considered in the collection
#' @param sampling The sampling design for missing data modeling : "dyad", "double_standard", "node", "snowball", "degree", "block" by default "undirected" is choosen
#' @param covarMatrix An optional matrix of covariates with dimension N x M (M covariates per node).
#' @param covarSimilarity An optional R x R -> R function  to compute similarity between node covariates. Default is #'
#' @param clusterInit Initial method for clustering: either a character in "hierarchical", "spectral" or "kmeans", or a list with \code{length(vBlocks)} vectors, each with size \code{ncol(adjacencyMatrix)} providing a user-defined clustering
#' @param trace logical, control the verbosity. Default to \code{TRUE}.
#' @param mc.cores integer, the number of cores to use when multiply model are fitted
#' @param smoothing character indicating what kind of ICL smoothing should be use among "none", "forward", "backward" or "both"
#' @param iter_both integer for the number of iteration in case of foward-backward (aka both) smoothing
#' @param control_VEM a list controlling the variational EM algorithm. See details.
#' @param Robject an object with class \code{missSBMcollection}
#' @return Returns an S3 object with class \code{missSBMcollection}, which is a list with all models estimated for all Q in vBlocks. \code{missSBMcollection} owns a couple of S3 methods: \code{is.missSBMcollection} to test the class of the object, a method \code{ICL} to extract the values of the Integrated Classification Criteria for each model, a method \code{getBestModel} which extract from the list the best model (and object of class \code{missSBM-fit}) according to the ICL, and a method \code{optimizationStatus} to monitor the objective function a convergence of the VEM algorithm.
#' @seealso \code{\link{sample}}, \code{\link{simulate}} and \code{\link{missingSBM_fit}}.
#' @examples
#' ## SBM parameters
#' directed <- FALSE
#' N <- 300 # number of nodes
#' Q <- 3   # number of clusters
#' alpha <- rep(1,Q)/Q     # mixture parameter
#' pi <- diag(.45,Q) + .05 # connectivity matrix
#'
#' ## simulate a SBM without covariates
#' sbm <- missSBM::simulate(N, alpha, pi, directed)
#'
#' ## Sample network data
#' samplingParameters <- .5 # the sampling rate
#' sampling <- "dyad"       # the sampling design
#' sampledNet <- missSBM::sample(sbm$adjMatrix, sampling, samplingParameters)
#'
#' ## Inference :
#' vBlocks <- 1:5 # number of classes
#' sbm <- missSBM::estimate(sampledNet$adjMatrix, vBlocks, sampling)
#' @import R6 parallel
#' @include utils_smoothing.R
#' @export
estimate <- function(
  adjacencyMatrix,
  vBlocks,
  sampling,
  clusterInit = ifelse(is.null(covarMatrix), "hierarchical", "spectral"),
  covarMatrix = NULL,
  covarSimilarity = l1_similarity,
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

  ## defaut control parameter for VEM, overwritten by user specification
  control <- list(threshold = 1e-4, maxIter = 200, fixPointIter = 5, trace = FALSE)
  control[names(control_VEM)] <- control_VEM

  ## Instatntiate the collection of missingSBM_fit
  myCollection <- missSBM_collection(adjacencyMatrix, vBlocks, sampling, clusterInit, covarMatrix, covarArray, trace)

  myCollection$estimate(control, mc.cores, trace)

  smoothing <- match.arg(smoothing)
  if (smoothing != "none") {
    myCollection$smooth_ICL(smoothing, control,  clusterInit, trace)
  }

  structure(setNames(models, vBlocks), class = "missSBMcollection")
}

