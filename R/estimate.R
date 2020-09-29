#' Estimation of SBMs with missing data
#'
#' Variational inference from sampled network data on a collection of
#' Stochastic Block Models indexed by block number.
#'
#' @param adjacencyMatrix The N x N adjacency matrix of the network to sample. If \code{adjacencyMatrix} is symmetric,
#' we assume an undirected network with no loop; otherwise the network is assumed directed.
#' @param vBlocks The vector of number of blocks considered in the collection
#' @param sampling The sampling design for the modelling of missing data: MAR designs ("dyad", "node","covar-dyad","covar-node","snowball")
#' and NMAR designs ("double-standard", "block-dyad", "block-node" ,"degree"). See details.
#' @param covariates A list with M entries (the M covariates). If the covariates are node-centered, each entry of \code{covariates}
#' must be a size-N vector;  if the covariates are dyad-centered, each entry of \code{covariates} must be N x N matrix.
#' @param control a list of parameters controlling the variational EM algorithm. See details.
#' @return Returns an R6 object with class \code{\link{missSBM_collection}}.
#'
#' @details The list of parameters \code{control} tunes more advanced features, such as the
#' initialization, how covariates are handled in the model, and the variational EM algorithm:
#'  \itemize{
#'  \item{"useCovSBM"}{logical. If covariates is not null, should they be used for the
#'         for the SBM inference (or just for the sampling)? Default is TRUE.}
#'  \item{"clusterInit"}{Initial method for clustering: either a character in "hierarchical", "spectral"
#'         or "kmeans", or a list with \code{length(vBlocks)} vectors, each with size
#'         \code{ncol(adjacencyMatrix)},  providing a user-defined clustering. Default is "hierarchical".}
#'  \item{"similarity"}{An R x R -> R function to compute similarities between node covariates. Default is \code{l1_similarity}, that is, -abs(x-y).}
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
#'     \item{"snowball": parameter estimation is not relevant for this sampling}
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
#' ## Sampling parameters
#' samplingParameters <- .5 # the sampling rate
#' sampling  <- "dyad"       # the sampling design
#'
#' ## simulate a SBM without covariate
#' sbm <- missSBM::simulate(N, alpha, pi, directed)
#'
#' ## Sample some dyads data + Infer SBM with missing data
#' collection <-
#'    missSBM::sample(sbm$adjacencyMatrix, sampling, samplingParameters) %>%
#'    missSBM::estimate(vBlocks = 1:5, sampling = sampling)
#' collection$ICL
#' coef(collection$bestModel$fittedSBM, "connectivity")
#'
#' myModel <- collection$bestModel
#' plot(myModel, "monitoring")
#' coef(myModel, "sampling")
#' coef(myModel, "connectivity")
#' head(predict(myModel))
#' head(fitted(myModel))
#'
#' @import R6 parallel
#' @export

estimate <- function(adjacencyMatrix, vBlocks, sampling, covariates = NULL, control = list()) {

  ## Sanity checks
  stopifnot(sampling %in% available_samplings)
  stopifnot(is.numeric(vBlocks))
  stopifnot(is.character(sampling))

  ## If no covariate is provided, you cannot ask for using them
  if (is.null(covariates)) control$useCovSBM <- FALSE
  ## If nothing specified by the user, use covariates by default
  else if (is.null(control$useCovSBM)) control$useCovSBM <- TRUE

  ## Defaut control parameters overwritten by user specification
  ctrl <- list(threshold = 1e-3, trace = 1, cores = 1, clusterInit = "hierarchical")
  if (control$useCovSBM) {
    stopifnot(sampling %in% available_samplings_covariates)
    ctrl <- c(ctrl,list(maxIter = 50, fixPointIter = 2, similarity = l1_similarity))
  } else {
    ctrl <- c(ctrl, list(maxIter = 100, fixPointIter = 5))
  }
  ctrl[names(control)] <- control

  ## Prepare network data for estimation with missing data
  sampledNet <- sampledNetwork$new(adjacencyMatrix, covariates, ctrl$similarity)

  ## Instantiate the collection of missSBM_fit
  myCollection <- missSBM_collection$new(
      sampledNet  = sampledNet,
      vBlocks     = vBlocks,
      sampling    = sampling,
      clusterInit = ctrl$clusterInit,
      cores       = ctrl$cores,
      trace       = (ctrl$trace > 0),
      useCov      = ctrl$useCovSBM
  )

  ## Launch estimation of each missSBM_fit
  myCollection$estimate(ctrl)

  ## Return the collection of optimized missSBM_fit
  myCollection
}

