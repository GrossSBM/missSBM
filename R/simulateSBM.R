#' Simulation of a Stochastic Block Model
#'
#' Generates realizations (blocks and adjancency) of a Stochastic Block model
#'
#' @param N The number of nodes
#' @param alpha The mixture parameters
#' @param pi The connectivity matrix (probabilities inter and intra clusters)
#' @param directed Boolean variable to indicate whether the network is directed or not,
#' by default "undirected" is choosen
#' @param covarParam An optional vector of parameters associated with the covariates, with size M
#' @param covarMatrix An optional matrix of covariates with dimension N x M (M covariates per node).
#' @param covarSimilarity An optional R x R -> R function  to compute similarity between node covariates. Default is \code{l1_similarity}, that is, -abs(x-y).
#' @return a object with class \code{SBM-Class}
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
#' @export
simulateSBM <- function(N, alpha, pi, directed = FALSE, covarMatrix = NULL, covarParam = NULL, covarSimilarity=l1_similarity){
  mySBM <- SBM_sampler$new(directed, N, alpha, pi, covarMatrix, covarParam, covarSimilarity)
  mySBM$rBlocks()
  mySBM$rAdjMatrix()
  mySBM
}
