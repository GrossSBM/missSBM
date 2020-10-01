#' Simulation of an SBM
#'
#' Generates a realization (blocks and adjacency matrix) of a Stochastic Block model
#'
#' @param nbNodes The number of nodes
#' @param blockProp The mixture parameters
#' @param connectParam The connectivity matrix (inter/intra clusters probabilities. provided on a logit scale for a model with covariates)
#' @param directed Boolean variable to indicate whether the network is directed or not. Default to \code{FALSE}.
#' @param covariates A list with M entries (the M covariates). Each entry of the list must be an N x N matrix.
#' @param covarParam An optional vector of parameters associated with the covariates, with size M
#' @return an object with class \code{SBM_sampler}
#' @seealso The class \code{\link{SBM_sampler}}
#' @examples
#' ## SBM parameters
#' directed <- FALSE
#' N <- 300 # number of nodes
#' Q <- 3   # number of clusters
#' M <- 2 # two Gaussian covariates
#' alpha <- rep(1, Q)/Q     # mixture parameters
#' theta <- diag(.45, Q) + .05 # connectivity matrix
#' eta   <- rnorm(M, -1, 1)  # covariate parametes
#' gamma <- log(theta/(1-theta)) # logit transform of theta for the model with covariates
#' X <- replicate(M, matrix(rnorm(N * N ,mean = 0, sd = 1), N, N), simplify = FALSE)
#'
#' ## draw a SBM without covariates
#' sbm <- missSBM::simulate(N, alpha, theta, directed)
#' coef(sbm, "connectivity")
#'
#' ## draw a SBM model with node-centered covariates
#' sbm_cov <- missSBM::simulate(N, alpha, gamma, directed, X, eta)
#' coef(sbm_cov, "covariates")
#'
#' old_param <- par(mfrow = c(1,2))
#' plot(sbm)
#' plot(sbm_cov)
#' par(old_param)
#'
#' @export
simulate <- function(nbNodes, blockProp, connectParam, directed = FALSE, covariates = list(), covarParam = numeric(length(covariates))) {

  ## Instantiation of the SBM sampler
  mySBM <-
    SBM_sampler$new(
      directed     = directed,
      nbNodes      = nbNodes,
      blockProp    = blockProp,
      connectParam = connectParam,
      covarParam   = covarParam,
      covarList    = covariates
    )

  ## draw blocks
  mySBM$rMemberships()

  ## draw adjacency matrix
  mySBM$rAdjacency()

  ## send the sampled SBM object
  mySBM
}
