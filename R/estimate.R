#' Inference of an SBM with missing data
#'
#' Perform variational inference of a Stochastic Block Model from a sampled adjacency matrix
#'
#' @param adjacencyMatrix The adjacency matrix of the network
#' @param vBlocks The vector of number of blocks considered in the collection
#' @param sampling The sampling design for missing data modeling : "dyad", "node", "double-standard", "block-dyad", "block-node" ,"degree
#' @param covariates A list with M entries (the M covariates). If the covariates are node-centred, each entry of \code{covariates}
#' must be a size-N vector;  if the covariates are dyad-centred, each entry of \code{covariates} must be N x N matrix.
#' @param similarity An optional R x R -> R function to compute similarities between node covariates. Default is \code{l1_similarity}, that is, -abs(x-y).
#' Only relevent when the covariates are node-centered (i.e. \code{covariates} is a list of size-N vectors).
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
  clusterInit = ifelse(is.null(covariates), "hierarchical", "spectral"),
  covariates = NULL,
  similarity = l1_similarity,
  trace = TRUE, cores = 1, control_VEM = list()) {

  ## defaut control parameter for VEM, overwritten by user specification
  control <- list(threshold = 1e-4, maxIter = 200, fixPointIter = 5, trace = FALSE)
  control[names(control_VEM)] <- control_VEM

  if (!is.null(covariates)) {
    stopifnot(sampling %in% available_samplings_covariates)
    # Conversion of covariates to an array
    covariates <- simplify2array(covariates)
    # if a list of vector (covariates node-centered), will be a matrix
    # and thus must be node centered
    if (is.matrix(covariates)) {
      stopifnot(sampling == "node")
      covarMatrix <- covariates
      covarArray  <- getCovarArray(covarMatrix, similarity)
    }
    # if a list of matrix (covariates dyad-centered), will be a 3-dimensional array
    # and thus must be dyad centered
    if (length(dim(covariates)) == 3) {
      stopifnot(sampling  == "dyad")
      covarMatrix <- NULL
      covarArray  <- covariates
    }
  } else {
    covarMatrix <- NULL
    covarArray  <- NULL
  }

  ## Instantiate the collection of missSBM_fit
  myCollection <- missSBM_collection$new(
      adjMatrix    = adjacencyMatrix,
      vBlocks      = vBlocks,
      sampling     = sampling,
      clusterInit  = clusterInit,
      covarMatrix  = covarMatrix,
      covarArray   = covarArray,
      cores        = cores,
      trace        = trace
  )

  ## Launch estimation of each missSBM_fit
  myCollection$estimate(control, cores, trace)

  ## Return the collection of optimized missSBM_fit
  myCollection
}

