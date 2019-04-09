#' Inference of an SBM with missing data
#'
#' Perform variational inference of a Stochastic Block Model from a sampled adjacency matrix
#'
#' @param sampledNet An object with class \code{\link{sampledNetwork}}, typically obtained wthe the function \code{\link{prepare_data}} or \code{\link{sample}}.
#' @param vBlocks The vector of number of blocks considered in the collection
#' @param sampling The sampling design for missing data modeling : "dyad", "node", "double-standard", "block-dyad", "block-node" ,"degree
#' @param covariates A list with M entries (the M covariates). If the covariates are node-centred, each entry of \code{covariates}
#' must be a size-N vector;  if the covariates are dyad-centred, each entry of \code{covariates} must be N x N matrix.
#' @param similarity An optional R x R -> R function to compute similarities between node covariates. Default is \code{l1_similarity}, that is, -abs(x-y).
#' Only relevent when the covariates are node-centered (i.e. \code{covariates} is a list of size-N vectors).
#' @param clusterInit Initial method for clustering: either a character in "hierarchical", "spectral" or "kmeans", or a list with \code{length(vBlocks)} vectors, each with size \code{ncol(adjacencyMatrix)} providing a user-defined clustering
#' @param control a list controlling the variational EM algorithm. See details.
#' @return Returns an R6 object with class \code{\link{missSBM_collection}}.
#'
#' @details The list of parameters \code{control} controls the optimziation process and the variational EM algorithm, with the following entries
#'  \itemize{
#'  \item{"threshold"}{stop when an optimization step changes the objective function by less than threshold. Default is 1e-4.}
#'  \item{"maxiter"}{V-EM algorithm stops when the number of iteration exceeds maxIter. Default is 200}
#'  \item{"fixPointIter"}{number of fix-point iteration for the Variational E step. Default is 5.}
#'  \item{"cores"}{integer for number of cores used. Default is 1.}
#'  \item{"trace"}{integer for verbosity (0, 1, 2). Useless when \code{cores} > 1}
#' }
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
#' collection <- missSBM::estimate(sampledNet, vBlocks, sampling)
#' collection$ICL
#' @import R6 parallel
#' @export
estimate <- function(sampledNet, vBlocks, sampling, clusterInit = "hierarchical", control = list()) {

  ## defaut control parameter for VEM, overwritten by user specification
  ctrl <- list(threshold = 1e-4, maxIter = 200, fixPointIter = 5, trace = 1, cores = 1)
  ctrl[names(control)] <- control

  ## Instantiate the collection of missSBM_fit
  myCollection <- missSBM_collection$new(
      sampledNet  = sampledNet,
      vBlocks     = vBlocks,
      sampling    = sampling,
      clusterInit = clusterInit,
      cores       = ctrl$cores,
      trace       = (ctrl$trace > 0)
  )

  ## Launch estimation of each missSBM_fit
  myCollection$estimate(ctrl)

  ## Return the collection of optimized missSBM_fit
  myCollection
}



#' Prepare network data for estimation with missing data
#'
#' This function put together the adjacency matrix of a network and an optional list of covariates
#' into a single \code{\link{sampledNetwork}} object, ready to use for inference with the \code{\link{estimate}}
#' function of the missSBM package.
#'
#' @param adjacencyMatrix The adjacency matrix of the network (NAs allowed)
#' @param covariates An optional list with M entries (the M covariates). If the covariates are node-centred, each entry of \code{covariates}
#' must be a size-N vector;  if the covariates are dyad-centred, each entry of \code{covariates} must be N x N matrix.
#' @param similarity An optional R x R -> R function to compute similarities between node covariates. Default is \code{l1_similarity}, that is, -abs(x-y).
#' Only relevent when the \code{covariates} is a list of size-N vectors.
#' @return Returns an R6 object with class \code{\link{sampledNetwork}}.
#'
#' @seealso \code{\link{estimate}} and \code{\link{sampledNetwork}}.
#' @importFrom igraph as_adj
#' @examples
#' data(war_graphs)
#' adj_beligerent <- war_graphs$beligerent %>% igraph::as_adj(sparse = FALSE)
#' sampledNet_war_nocov <- prepare_data(adj_beligerent)
#' military_power <- igraph::get.vertex.attribute(war_graphs$beligerent)$military_power
#' sampledNet_war_withcov <- prepare_data(adj_beligerent, list(military_power = military_power))
#' @export
prepare_data <- function(adjacencyMatrix, covariates = NULL, similarity = missSBM:::l1_similarity) {

  covar <- format_covariates(covariates, similarity)

  sampledNet <- sampledNetwork$new(adjacencyMatrix, covar$Matrix, covar$Array)

  sampledNet
}
