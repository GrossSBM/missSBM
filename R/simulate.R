#' Simulation of an SBM
#'
#' Generates a realization (blocks and adjancency matrix) of a Stochastic Block model
#'
#' @param N The number of nodes
#' @param alpha The mixture parameters
#' @param pi The connectivity matrix (either probabilities of inter and intra clusters )
#' @param directed Boolean variable to indicate whether the network is directed or not. Default to \code{FALSE}.
#' @param covarParam An optional vector of parameters associated with the covariates, with size M
#' @param covarMatrix An optional matrix of covariates with dimension N x M (M covariates per node).
#' @param covarSimilarity An optional R x R -> R function  to compute similarity between node covariates. Default is \code{l1_similarity}, that is, -abs(x-y).
#' @return an object with class \code{SBM_sampler}
#' @seealso The class \code{\link{SBM_sampler}}
#' @examples
#' ## SBM parameters
#' directed <- FALSE
#' N <- 300 # number of nodes
#' Q <- 3   # number of clusters
#' alpha <- rep(1,Q)/Q     # mixture parameter
#' pi <- diag(.45,Q) + .05 # connectivity matrix
#' gamma <- log(pi/(1-pi)) # logit transform fo the model with covariates
#' M <- 2 # two Gaussian covariates
#' covarMatrix <- matrix(rnorm(N*M,mean = 0, sd = 1), N, M)
#' covarParam  <- rnorm(M, -1, 1)
#'
#' ## draw a SBM without covariates
#' sbm <- missSBM::simulate(N, alpha, pi, directed)
#'
#' ## draw a SBM model with covariates
#' sbm_cov <- missSBM::simulate(N, alpha, gamma, directed, covarMatrix, covarParam)
#'
#' \dontrun{
#' par(mfrow = c(1,2))
#' plot(sbm)
#' plot(sbm_cov)
#' }
#'
#' @export
simulate <- function(N, alpha, pi, directed = FALSE, covarMatrix = NULL, covarParam = NULL, covarSimilarity=l1_similarity){
  mySBM <-
    SBM_sampler$new(
      directed     = directed,
      nNodes       = N,
      mixtureParam = alpha,
      connectParam = pi,
      covarParam   = covarParam,
      covarArray   = getCovarArray(covarMatrix, covarSimilarity)
    )
  mySBM$rBlocks()
  mySBM$rAdjMatrix()
  mySBM
}
