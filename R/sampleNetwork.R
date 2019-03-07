#' Sampling of a network
#'
#' Samples a matrix (the adjacency matrix of a network) under the SBM
#'
#' @param adjMatrix The adjacency matrix of the network
#' @param sampling The sampling design used to sample the adjacency matrix
#' @param parameters The sampling parameters adapted to each sampling
#' @param clusters Clusters membership vector of the nodes, only necessary for class sampling, by default equal to
#' @param covarMatrix An optional matrix of covariates with dimension N x M (M covariates per node).
#' @param covarSimilarity An optional R x R -> R function  to compute similarity between node covariates. Default is
#' @return \code{sampleNetwork} returns a matrix (the sampled adjacency matrix of the network given in parameter)
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
#'    sampleNetwork(
#'      mySBM$adjMatrix,
#'      sampling,
#'      samplingParameters
#'    )
#' @export
sampleNetwork <- function(adjMatrix, sampling, parameters, clusters = NULL, covarMatrix = NULL, covarSimilarity = l1_similarity){

  stopifnot(sampling %in% available_samplings)
  if (!is.null(covarMatrix)) stopifnot(sampling %in% available_samplings_covariates)

  N <- ncol(adjMatrix)
  directed <- !isSymmetric(adjMatrix)

  ## instantiate the sampler
  mySampler <-
    switch(sampling,
      "dyad"            = simpleDyadSampler$new(parameters, N, directed, getCovarArray(covarMatrix, covarSimilarity)),
      "node"            = simpleNodeSampler$new(parameters, N, directed, covarMatrix),
      "double-standard" = doubleStandardSampler$new(parameters, adjMatrix, directed),
      "block-dyad"      = blockDyadSampler$new(parameters, N, directed, clusters),
      "block-node"      = blockNodeSampler$new(parameters, N, directed, clusters),
      "degree"          = degreeSampler$new(parameters, rowSums(adjMatrix), directed),
  )
  ## draw a sampling matrix R
  mySampler$rSamplingMatrix()

  ## turn this matrix to a sampled Network object
  adjMatrix[mySampler$samplingMatrix == 0] <- NA
  sampledNet <- sampledNetwork$new(adjMatrix)
  sampledNet
}
