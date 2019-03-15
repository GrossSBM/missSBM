#' Sampling of a network data set
#'
#' Samples a matrix (the adjacency matrix of a network), possible generated under a stochastic block model
#'
#' @param adjacencyMatrix The N x N adjacency matrix of the network to sample
#' @param sampling The sampling design used to sample the adjacency matrix
#' @param parameters The sampling parameters adapted to each sampling
#' @param clusters An optional clustering membership vector of the nodes, only necessary for class sampling
#' @param covarMatrix An optional matrix of covariates with dimension N x M (M covariates per node).
#' @param covarSimilarity An optional R x R -> R function to compute the similarity between node covariates. The default internal function missSBM:::l1_smilarity is based on the oppositite of the absolute difference between two vector of covariates.
#'
#' @return \code{sampleNetwork} returns an object with class \code{\link{sampledNetwork}} containing all the useful information about the sampling.
#' @seealso The class \code{\link{sampledNetwork}}
#'
#' @details The differents sampling designs are splitted into two families in which we find dyad-centered and node-centered samplings. See \cite{1} for details.
#' \itemize{
#' \item Missing At Random (MAR)
#'   \itemize{
#'     \item{"dyad": parameter = p \deqn{p = P(Dyad (i,j) is sampled)}}
#'     \item{"node": parameter = p and \deqn{p = P(Node i is sampled)}}
#'   }
#' \item Not Missing At Random (NMAR)
#'   \itemize{
#'     \item{"double-standard": parameter = (p0,p1) and \deqn{p0 = P(Dyad (i,j) is sampled | the dyad is equal to 0)=}, p1 = P(Dyad (i,j) is sampled | the dyad is equal to 1)}
#'     \item{"block-node": parameter = c(p(1),...,p(Q)) and \deqn{p(q) = P(Node i is sampled | node i is in cluster q)}}
#'     \item{"block-dyad": parameter = c(p(1,1),...,p(Q,Q)) and \deqn{p(q,l) = P(Edge (i,j) is sampled | node i is in cluster q and node j is in cluster l)}}
#'     \item{"degree": parameter = c(a,b) and \deqn{logit(a+b*Degree(i)) = P(Node i is sampled | Degree(i))}}
#'   }
#' }
#' @examples
#' ## SBM parameters
#' directed <- FALSE
#' N <- 300 # number of nodes
#' Q <- 3   # number of clusters
#' alpha <- rep(1,Q)/Q     # mixture parameter
#' pi <- diag(.45,Q) + .05 # connectivity matrix
#'
#' ## draw a SBM without covariates
#' sbm <- simulateSBM(N, alpha, pi, directed)
#'
#' ## Sampling of the network data
#'
#' # some sampling design and their associated parameters
#' sampling_parameters <- list(
#'    "dyad" = .3,
#'    "node" = .3,
#'    "double-standard" = c(0.4, 0.8),
#'    "block-node" = c(.3, .8, .5),
#'    "block-dyad" = pi,
#'    "degree" = c(.01, .01)
#'  )
#'
#' sampled_networks <- list()
#'
#' for (sampling in names(sampling_parameters)) {
#'   sampled_networks[[sampling]] <-
#'      sampleNetwork(
#'        adjacencyMatrix = sbm$adjMatrix,
#'        sampling        = sampling,
#'        parameters      = sampling_parameters[[sampling]],
#'        cluster         = sbm$memberships
#'      )
#' }
#' \dontrun{
#' par(mfrow = c(2,3))
#' for (sampling in names(sampling_parameters)) {
#'   plot(sampled_networks[[sampling]] , main = paste(sampling, "sampling"))
#' }
#' par(mfrow = c(1,1))
#' }
#' @export
sampleNetwork <- function(adjacencyMatrix, sampling, parameters, clusters = NULL, covarMatrix = NULL, covarSimilarity = l1_similarity){

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
      "degree"          = degreeSampler$new(parameters, rowSums(adjacencyMatrix), directed)
  )
  ## draw a sampling matrix R
  mySampler$rSamplingMatrix()

  ## turn this matrix to a sampled Network object
  adjacencyMatrix[mySampler$samplingMatrix == 0] <- NA
  sampledNet <- sampledNetwork$new(adjacencyMatrix)
  sampledNet
}
