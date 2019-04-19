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
#' data(war)
#' adj_beligerent <- war$beligerent %>% igraph::as_adj(sparse = FALSE)
#' sampledNet_war_nocov <- prepare_data(adj_beligerent)
#' sampledNet_war_withcov <- prepare_data(adj_beligerent, list(military_power = war$beligerent$power))
#' @export
prepare_data <- function(adjacencyMatrix, covariates = NULL, similarity = missSBM:::l1_similarity) {

  covar <- format_covariates(covariates, similarity)

  sampledNet <- sampledNetwork$new(adjacencyMatrix, covar$Matrix, covar$Array)

  sampledNet
}
