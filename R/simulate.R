#' Simulation of an SBM
#'
#' Generates a realization (blocks and adjancency matrix) of a Stochastic Block model
#'
#' @param nNodes The number of nodes
#' @param mixtureParam The mixture parameters
#' @param connectParam The connectivity matrix (inter/intra clusters probabilities. provided on a logit scale for a model with covariates)
#' @param directed Boolean variable to indicate whether the network is directed or not. Default to \code{FALSE}.
#' @param covariates A list with M entries (the M covariates). If the covariates are node-centred, each entry of \code{covariates}
#' must be a size-N vector;  if the covariates are dyad-centred, each entry of \code{covariates} must be N x N matrix.
#' @param covarParam An optional vector of parameters associated with the covariates, with size M
#' @param similarity An optional R x R -> R function to compute similarities between node covariates. Default is \code{l1_similarity}, that is, -abs(x-y). Only relevent when covariates is a list of size-N vectors.
#' @return an object with class \code{SBM_sampler}
#' @seealso The class \code{\link{SBM_sampler}}
#' @examples
#' ## SBM parameters
#' directed <- FALSE
#' N <- 300 # number of nodes
#' Q <- 3   # number of clusters
#' M <- 2 # two Gaussian covariates
#' alpha <- rep(1,Q)/Q     # mixture parameters
#' pi <- diag(.45,Q) + .05 # connectivity matrix
#' eta <- rnorm(M, -1, 1)  # covaraite parametes
#' gamma <- log(pi/(1-pi)) # logit transform of pi for the model with covariates
#' X   <- replicate(M, rnorm(N,mean = 0, sd = 1), simplify = FALSE)
#' phi <- replicate(M, matrix(rnorm(N * N ,mean = 0, sd = 1), N, N), simplify = FALSE)

#' ## draw a SBM without covariates
#' sbm <- missSBM::simulate(N, alpha, pi, directed)
#'
#' ## draw a SBM model with node-centred covariates
#' sbm_cov1 <- missSBM::simulate(N, alpha, gamma, directed, X, eta)
#'
#' ## draw a SBM model with dyad-centered covariates
#' sbm_cov2 <- missSBM::simulate(N, alpha, gamma, directed, X, eta)
#'
#' \dontrun{
#' par(mfrow = c(1,3))
#' plot(sbm)
#' plot(sbm_cov1)
#' plot(sbm_cov2)
#' }
#'
#' @export
simulate <- function(nNodes, mixtureParam, connectParam, directed = FALSE, covariates = NULL, covarParam = NULL, similarity = l1_similarity) {

  ## Conversion of covariates to an array
  if (!is.null(covariates)) {
    covariates <- simplify2array(covariates)
    if (is.matrix(covariates))
      covariates <- getCovarArray(covariates, similarity)
  }

  ## Instantiation of the SBM sampler
  mySBM <-
    SBM_sampler$new(
      directed     = directed,
      nNodes       = nNodes,
      mixtureParam = mixtureParam,
      connectParam = connectParam,
      covarParam   = covarParam,
      covarArray   = covariates
    )

  ## draw blocks
  mySBM$rBlocks()

  ## draw adjacency matrix
  mySBM$rAdjMatrix()

  ## send the sampled SBM object
  mySBM
}
