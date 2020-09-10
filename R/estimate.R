#' Estimation of SBMs with missing data
#'
#' Variational inference from sampled network data on a collection of
#' Stochastic Block Models indexed by block number.
#'
#' @param sampledNet An object with class \code{\link{sampledNetwork}}, typically obtained with
#' the function \code{\link{prepare_data}} (real-word data) or \code{\link{sample}} (simulation).
#' @param vBlocks The vector of number of blocks considered in the collection
#' @param sampling The sampling design for the modelling of missing data: MAR designs ("dyad", "node")
#' and NMAR designs ("double-standard", "block-dyad", "block-node" ,"degree")
#' @param clusterInit Initial method for clustering: either a character in "hierarchical", "spectral"
#' or "kmeans", or a list with \code{length(vBlocks)} vectors, each with size \code{ncol(adjacencyMatrix)},
#' providing a user-defined clustering. Default is "hierarchical".
#' @param useCovariates logical. If covariates are present in sampledNet, should they be used for the inference or of the network sampling design, or just for the SBM inference? default is TRUE.
#' @param control a list of parameters controlling the variational EM algorithm. See details.
#' @return Returns an R6 object with class \code{\link{missSBM_collection}}.
#'
#' @details The list of parameters \code{control} essentially tunes the optimization process and the
#' variational EM algorithm, with the following parameters
#'  \itemize{
#'  \item{"threshold"}{stop when an optimization step changes the objective function by less than threshold. Default is 1e-4.}
#'  \item{"maxIter"}{V-EM algorithm stops when the number of iteration exceeds maxIter. Default is 200}
#'  \item{"fixPointIter"}{number of fix-point iterations in the Variational E step. Default is 5.}
#'  \item{"cores"}{integer for number of cores used. Default is 1.}
#'  \item{"trace"}{integer for verbosity (0, 1, 2). Default is 1. Useless when \code{cores} > 1}
#' }
#'
#' @details The different sampling designs are split into two families in which we find dyad-centered and
#' node-centered samplings. See <doi:10.1080/01621459.2018.1562934> for complete description.
#' \itemize{
#' \item Missing at Random (MAR)
#'   \itemize{
#'     \item{"dyad": parameter = p and \deqn{p = P(Dyad (i,j) is sampled)}}
#'     \item{"node": parameter = p and \deqn{p = P(Node i is sampled)}}
#'     \item{"covar-dyad": parameter = beta in R^M and \deqn{P(Dyad (i,j) is sampled) = logistic(parameter' covarArray (i,j, ))}}
#'     \item{"covar-node": parameter = nu in R^M and \deqn{P(Node i is sampled)  = logistic(parameter' covarMatrix (i,)}}
#'   }
#' \item Not Missing At Random (NMAR)
#'   \itemize{
#'     \item{"double-standard": parameter = (p0,p1) and \deqn{p0 = P(Dyad (i,j) is sampled | the dyad is equal to 0)=}, p1 = P(Dyad (i,j) is sampled | the dyad is equal to 1)}
#'     \item{"block-node": parameter = c(p(1),...,p(Q)) and \deqn{p(q) = P(Node i is sampled | node i is in cluster q)}}
#'     \item{"block-dyad": parameter = c(p(1,1),...,p(Q,Q)) and \deqn{p(q,l) = P(Edge (i,j) is sampled | node i is in cluster q and node j is in cluster l)}}
#'     \item{"degree": parameter = c(a,b) and \deqn{logit(a+b*Degree(i)) = P(Node i is sampled | Degree(i))}}
#'   }
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
#' coef(collection$bestModel$fittedSBM, "connectivity")
#'
#' myModel <- collection$bestModel
#' plot(myModel, "monitoring")
#' coef(myModel, "sampling")
#' coef(myModel, "connectivity")
#' head(predict(myModel))
#' head(fitted(myModel))
#' @import R6 parallel
#' @export
estimate <- function(sampledNet, vBlocks, sampling, clusterInit = "hierarchical", useCovariates = TRUE, control = list()) {

  ## Sanity checks
  stopifnot(sampling %in% available_samplings)

  ## If no covariates, you don't have to use them
  if (is.null(sampledNet$covarArray)) useCovariates <- FALSE

  ## Defaut control parameters for VEM, overwritten by user specification
  if (useCovariates) {
    stopifnot(sampling %in% available_samplings_covariates)
    ctrl <- list(threshold = 1e-3, maxIter = 50, fixPointIter = 2, trace = 1, cores = 1)
  } else {
    ctrl <- list(threshold = 1e-3, maxIter = 100, fixPointIter = 5, trace = 1, cores = 1)
  }
  ctrl[names(control)] <- control

  ## Instantiate the collection of missSBM_fit
  myCollection <- missSBM_collection$new(
      sampledNet  = sampledNet,
      vBlocks     = vBlocks,
      sampling    = sampling,
      clusterInit = clusterInit,
      cores       = ctrl$cores,
      trace       = (ctrl$trace > 0),
      useCov      = useCovariates
  )

  ## Launch estimation of each missSBM_fit
  myCollection$estimate(ctrl)

  ## Return the collection of optimized missSBM_fit
  myCollection
}

