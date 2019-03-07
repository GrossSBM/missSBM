#' @title Simulation of a Stochastic Block Model
#'
#' @description \code{simulateSBM} is a function that generates a matrix (the adjacency matrix of a network) under the SBM
#'
#' @param N The number of nodes
#' @param alpha The mixture parameters
#' @param pi The connectivity matrix (probabilities inter and intra clusters)
#' @param directed Boolean variable to indicate whether the network is directed or not,
#' by default "undirected" is choosen
#' @param covarParam An optional vector of parameters associated with the covariates, with size M
#' @param covarMatrix An optional matrix of covariates with dimension N x M (M covariates per node).
#' @param covarSimilarity An optional R x R -> R function  to compute similarity between node covariates. Default is \code{l1_similarity}, that is, -abs(x-y).
#' @return \code{simulateSBM} returns a vector with clusters of nodes and a matrix (the adjacency matrix of the network)
#' @references [1] Tabouy, P. Barbillon, J. Chiquet. Variationnal inference of Stochastic Block Model from sampled data (2017). arXiv:1707.04141.
#' @seealso \code{\link{inferSBM}} and \code{\link{samplingSBM}}
#' @details The emission law can be :\itemize{\item{Bernoulli:
#' \deqn{P(Y[i,j] = 1 | Zi = q, Zj = l) = p(Zi,Zj)}}
#' \item{Poisson:
#' \deqn{P(Y[i,j] = k | Zi = q, Zj = l) = (\lambda(Z_i,Z_j)^k/(k!)) * exp(-\lambda(Z_i,Z_j))}}
#' }
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
simulateSBM <- function(N, alpha, pi, directed = FALSE, covariates = NULL, covarParam = NULL, covarSimilarity=l1_similarity){
  mySBM <- SBM_sampler$new(directed, N, alpha, pi, covariates, covarParam, covarSimilarity)
  mySBM$rBlocks()
  mySBM$rAdjMatrix()
  mySBM
}

#' @title Sampling of a network
#'
#' @description \code{samplingSBM} is a function that sample a matrix (the adjacency matrix of a network) under the SBM
#'
#' @param adjacencyMatrix The adjacency matrix of the network
#' @param sampling The sampling design used to sample the adjacency matrix
#' @param parameters The sampling parameters adapted to each sampling
#' @param clusters Clusters membership vector of the nodes, only necessary for class sampling, by default equal to NULL
#' @return \code{samplingSBM} returns a matrix (the sampled adjacency matrix of the network given in parameter)
#' @references [1] Tabouy, P. Barbillon, J. Chiquet. Variationnal inference of Stochastic Block Model from sampled data (2017). arXiv:1707.04141.
#' @seealso \code{\link{inferSBM}} and \code{\link{samplingSBM}}
#' @details The differents sampling designs are splitted into two families in which we find dyad-centered and node-centered samplings, for
#' more details see (\cite{1}) :\itemize{\item Missing At Random (MAR) \itemize{\item{"dyad": parameter = p
#' \deqn{p = P(Dyad (i,j) is sampled)}}
#' \item{"node": parameter = p and
#' \deqn{p = P(Node i is sampled)}}
#' \item{"snowball" (one step):
#' like the MARNode sampling plus we sample neighbours of nodes sampled at the first batch}
#' }
#' \item Not Missing At Random (NMAR) \itemize{ \item{"double_standard": parameter = (p0,p1) and
#' \deqn{p0 = P(Dyad (i,j) is sampled | the dyad is equal to 0)=}, p1 = P(Dyad (i,j) is sampled | the dyad is equal to 1)}
#' \item{"block_dyad": parameter = c(p(1,1),...,p(Q,Q)) and
#' \deqn{p(q,l) = P(Edge (i,j) is sampled | node i is in cluster q and node j is in cluster l)}}
#' \item{"degree": parameter = c(a,b) and
#' \deqn{logit(a+b*Degree(i)) = P(Node i is sampled | Degree(i))}}
#' \item{"block": parameter = c(p(1),...,p(Q)) and
#' \deqn{p(q) = P(Node i is sampled | node i is in cluster q)}}
#' }}
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
#' ## Sampling of the data : ##
#' samplingParameters <- .5 # the sampling rate
#' sampling <- "dyad"       # the sampling design
#' sampledNetwork <-
#'    samplingSBM(
#'      mySBM$adjMatrix,
#'      sampling,
#'      samplingParameters
#'    )
#' @export
samplingSBM <- function(adjacencyMatrix, sampling, parameters, clusters = NULL, covarMatrix = NULL, covarSimilarity = l1_similarity){

  stopifnot(sampling %in% available_samplings)
  if (!is.null(covarMatrix)) stopifnot(sampling %in% available_samplings_covariates)

  N <- ncol(adjacencyMatrix)
  directed <- !isSymmetric(adjacencyMatrix)

  ## instantiate the sampler
  mySampler <-
    switch(sampling,
      "dyad"            = simpleDyadSampler$new(parameters, N, directed, getCovarArray(covarMatrix, covarSimilarity)),
      "node"            = simpleNodeSampler$new(parameters, N, directed, covarMatrix),
      "double-standard" = doubleStandardSampler$new(parameters, adjacencyMatrix, directed),
      "block-dyad"      = blockDyadSampler$new(parameters, N, directed, clusters),
      "block-node"      = blockNodeSampler$new(parameters, N, directed, clusters),
      "degree"          = degreeSampler$new(parameters, rowSums(adjacencyMatrix), directed),
  )
  ## draw a sampling matrix R
  mySampler$rSamplingMatrix()

  ## turn this matrix to a sampled Network object
  adjacencyMatrix[mySampler$samplingMatrix == 0] <- NA
  sampledNet <- sampledNetwork$new(adjacencyMatrix)
  sampledNet
}

#' @title Inference of Stochastic Block Model from sampled data
#'
#' @description \code{inferSBM} is a function that makes variationnal inference of Stochastic Block Model from sampled adjacency matrix
#'
#' @param adjacencyMatrix The adjacency matrix of the network
#' @param vBlocks The vector of number of blocks considered in the collection
#' @param sampling The sampling design for missing data modeling : "dyad", "double_standard", "node", "snowball", "degree", "block"
#' by default "undirected" is choosen
#' @param clusterInit character in "hierarchical" or "spectral" for initialization
#' @param trace logical, control the verbosity. Default to \code{TRUE}.
#' @param mc.cores integer, the number of cores to use when multiply model are fitted
#' @param smoothing character indicating what kind of ICL smoothing should be use among "none", "forward", "backward" or "both"
#' @param iter_both integer for the number of iteration in case of foward-backward (aka both) smoothing
#' @param control_VEM a list controlling the variational EM algorithm. See details.
#' @param Robject an object with class \code{missSBMcollection}
#' @return \code{inferSBM} returns an S3 object with class \code{missSBMcollection}, which is a list with all models estimated for all Q in vBlocks. \code{missSBMcollection} owns a couple of S3 methods: \code{is.missSBMcollection} to test the class of the object, a method \code{ICL} to extract the values of the Integrated Classification Criteria for each model, a method \code{getBestModel} which extract from the list the best model (and object of class \code{missSBM-fit}) according to the ICL, and a method \code{optimizationStatus} to monitor the objective function a convergence of the VEM algorithm.
#' @references [1] Tabouy, P. Barbillon, J. Chiquet. Variationnal inference of Stochastic Block Model from sampled data (2017). arXiv:1707.04141.
#' @seealso \code{\link{samplingSBM}} and \code{\link{simulateSBM}} and \code{\link{missingSBM_fit}}.
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
#' ## Sampling of the data : ##
#' samplingParameters <- .5 # the sampling rate
#' sampling <- "dyad"       # the sampling design
#' sampledNetwork <-
#'    samplingSBM(
#'      mySBM$adjMatrix,
#'      sampling,
#'      samplingParameters
#'    )
#'
#' ## Inference :
#' vBlocks <- 1:5 # number of classes
#' sbm <-
#'    inferSBM(
#'       sampledNetwork$adjMatrix,
#'       vBlocks,
#'       sampling
#'    )
#' @import R6 parallel
#' @include smoother_SBM.R
#' @export
inferSBM <- function(
  adjacencyMatrix,
  vBlocks,
  sampling,
  clusterInit = "hierarchical",
  trace     = TRUE,
  smoothing = c("none", "forward", "backward", "both"),
  mc.cores = 1,
  iter_both = 1,
  control_VEM = list()) {

  ## some sanity checks
  try(
    !all.equal(unique(as.numeric(adjacencyMatrix[!is.na(adjacencyMatrix)])), c(0,1)),
    stop("Only binary graphs are supported.")
  )

  sampledNet <- sampledNetwork$new(adjacencyMatrix)
  if (trace) cat("\n")
  if (trace) cat("\n Adjusting Variational EM for Stochastic Block Model\n")
  if (trace) cat("\n\tImputation assumes a '", sampling,"' network-sampling process\n", sep = "")
  if (trace) cat("\n")
  models <- mclapply(vBlocks,
    function(nBlocks) {
    if (trace) cat(" Initialization of model with", nBlocks,"blocks.", "\r")
      if (is.list(clusterInit)) {
        missingSBM_fit$new(sampledNet, nBlocks, sampling, clusterInit[[nBlocks]])
      } else {
        missingSBM_fit$new(sampledNet, nBlocks, sampling, clusterInit)
      }
    }, mc.cores = mc.cores
  )

  ## defaut control parameter for VEM, overwritten by user specification
  control <- list(threshold = 1e-4, maxIter = 200, fixPointIter = 5, trace = FALSE)
  control[names(control_VEM)] <- control_VEM
  cat("\n")
  mclapply(models,
    function(model) {
      if (trace) cat(" Performing VEM inference for model with", model$fittedSBM$nBlocks,"blocks.\r")
      model$doVEM(control)
    }, mc.cores = mc.cores
  )

  smoothing <- match.arg(smoothing)
  if (smoothing != "none") {
    if (trace) cat("\n Smoothing ICL\n")
    smoothing_fn <- switch(smoothing,
                           "forward"  = smoothingForward ,
                           "backward" = smoothingBackward,
                           "both"     = smoothingForBackWard
    )
    if (!is.character(clusterInit)) {
      split_fn <- init_hierarchical
    } else {
      split_fn <- switch(clusterInit,
                         "spectral" = init_spectral,
                         "hierarchical" = init_hierarchical,
                         init_hierarchical)
    }
    control$trace <- FALSE # forcing no trace while smoothing
    models <- smoothing_fn(models, vBlocks, sampledNet, sampling, split_fn, mc.cores, iter_both, control)

  }

  structure(setNames(models, vBlocks), class = "missSBMcollection")
}

#' @rdname inferSBM
#' @export
is.missSBMcollection <- function(Robject) {
  inherits(Robject, "missSBMcollection")
}

#' @rdname inferSBM
#' @export
ICL <- function(Robject) { UseMethod("ICL", Robject) }

#' @rdname inferSBM
#' @importFrom stats setNames
#' @export
ICL.missSBMcollection <- function(Robject) {
  stopifnot(is.missSBMcollection(Robject))
  setNames(sapply(Robject, function(model) model$vICL), names(Robject))
}

#' @rdname inferSBM
#' @export
optimizationStatus <- function(Robject) { UseMethod("optimizationStatus", Robject) }

#' @rdname inferSBM
#' @export
optimizationStatus.missSBMcollection <- function(Robject) {
  stopifnot(is.missSBMcollection(Robject))
  Reduce("rbind",
    lapply(Robject,
      function(model) {
        res <- model$monitoring
        res$nBlock <- model$fittedSBM$nBlocks
        res
    })
  )
}

#' @rdname inferSBM
#' @export
getBestModel <- function(Robject) {UseMethod("getBestModel", Robject)}

#' @rdname inferSBM
#' @export
getBestModel.missSBMcollection <- function(Robject) {
  stopifnot(is.missSBMcollection(Robject))
  Robject[[which.min(ICL(Robject))]]
}

