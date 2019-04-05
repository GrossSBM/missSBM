#' Inference of an SBM with missing data
#'
#' Perform variational inference of a Stochastic Block Model from a sampled adjacency matrix
#'
#' @param adjacencyMatrix The adjacency matrix of the network
#' @param vBlocks The vector of number of blocks considered in the collection
#' @param sampling The sampling design for missing data modeling : "dyad", "node", "double-standard", "block-dyad", "block-node" ,"degree
#' @param covarMatrix An optional matrix of covariates with dimension N x M (M covariates per node).
#' @param covarSimilarity An optional R x R -> R function  to compute similarity between node covariates. Default is #'
#' @param clusterInit Initial method for clustering: either a character in "hierarchical", "spectral" or "kmeans", or a list with \code{length(vBlocks)} vectors, each with size \code{ncol(adjacencyMatrix)} providing a user-defined clustering
#' @param trace logical, control the verbosity. Default to \code{TRUE}.
#' @param cores integer, the number of cores to use when multiply model are fitted
#' @param control_VEM a list controlling the variational EM algorithm. See details.
#' @return Returns an R6 object with class \code{\link{missSBM_collection}}.
#' @seealso \code{\link{sample}}, \code{\link{simulate}}, \code{\link{missSBM_collection}} and \code{\link{missSBM_fit}}.
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
#' sampledNet <- missSBM::sample(sbm$adjacencyMatrix, sampling, samplingParameters)
#'
#' ## Inference :
#' vBlocks <- 1:5 # number of classes
#' collection <- missSBM::estimate(sampledNet$adjacencyMatrix, vBlocks, sampling, trace = FALSE)
#' collection$ICL
#' @import R6 parallel
#' @export
estimate <- function(adjacencyMatrix, vBlocks, sampling,
  clusterInit = ifelse(is.null(covarMatrix), "hierarchical", "spectral"),
  covarMatrix = NULL,
  covarSimilarity = l1_similarity,
  trace = TRUE, cores = 1, control_VEM = list()) {

  ## defaut control parameter for VEM, overwritten by user specification
  control <- list(threshold = 1e-4, maxIter = 200, fixPointIter = 5, trace = FALSE)
  control[names(control_VEM)] <- control_VEM

  ## Instantiate the collection of missSBM_fit
  myCollection <- missSBM_collection$new(
      adjMatrix       = adjacencyMatrix,
      vBlocks         = vBlocks,
      sampling        = sampling,
      clusterInit     = clusterInit,
      covarMatrix     = covarMatrix,
      covarSimilarity = covarSimilarity,
      cores           = cores,
      trace           = trace
  )

  ## Launch estimation of each missSBM_fit
  myCollection$estimate(control, cores, trace)

  ## Return the collection of optimized missSBM_fit
  myCollection
}

