#' Sampling of network data
#'
#' Samples a matrix (the adjacency matrix of a network), possible generated under a stochastic block model
#'
#' @param adjacencyMatrix The N x N adjacency matrix of the network to sample. If \code{adjacencyMatrix} is symmetric,
#' we assume an undirected network with no loop; otherwise the network is assumed directed.
#' @param sampling The sampling design used to sample the adjacency matrix
#' @param parameters The sampling parameters adapted to each sampling
#' @param clusters An optional clustering membership vector of the nodes, only necessary for class sampling
#' @param covariates A list with M entries (the M covariates). If the covariates are node-centred, each entry of \code{covariates}
#' must be a size-N vector;  if the covariates are dyad-centred, each entry of \code{covariates} must be N x N matrix.
#' @param similarity An optional R x R -> R function to compute similarities between node covariates. Default is \code{l1_similarity}, that is, -abs(x-y).
#' Only relevent when the covariates are node-centered (i.e. \code{covariates} is a list of size-N vectors).
#'
#' @return an object with class \code{\link{sampledNetwork}} containing all the useful information about the sampling.
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
#' ## simulate a SBM without covariates
#' sbm <- missSBM::simulate(N, alpha, pi, directed)
#'
#' ## Sample network data
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
#'      missSBM::sample(
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
sample <- function(adjacencyMatrix, sampling, parameters, clusters = NULL, covariates = NULL, similarity = l1_similarity) {

  ## Turn the covariates to an array if not null
  if (!is.null(covariates)) {
    stopifnot(sampling %in% available_samplings_covariates)
    # Conversion of covariates to an array
    covariates <- simplify2array(covariates)
    # if a list of vector (covariates node-centered), will be a matrix
    if (is.matrix(covariates)) stopifnot(sampling == "node")
    # if a list of matrix (covariates dyad-centered), will be a 3-dimensional array
    if (length(dim(covariates)) == 3) stopifnot(sampling  == "dyad")
  } else {
    stopifnot(sampling %in% available_samplings)
  }

  nNodes   <- ncol(adjacencyMatrix)
  directed <- !isSymmetric(adjacencyMatrix)

  ## SAMPLER INSTANTIATION
  mySampler <-
    switch(sampling,
      "dyad"       = simpleDyadSampler$new(
        parameters = parameters, nNodes = nNodes, directed = directed, covarArray  = covariates),
      "node"       = simpleNodeSampler$new(
        parameters = parameters, nNodes = nNodes, directed = directed, covarMatrix = covariates),
      "double-standard" = doubleStandardSampler$new(
        parameters = parameters, adjMatrix = adjacencyMatrix, directed = directed),
      "block-dyad" = blockDyadSampler$new(
        parameters = parameters, nNodes = nNodes, directed = directed, clusters = clusters),
      "block-node" = blockNodeSampler$new(
        parameters = parameters, nNodes = nNodes, directed = directed, clusters = clusters),
      "degree"     = degreeSampler$new(
        parameters = parameters, degrees = rowSums(adjacencyMatrix), directed = directed)
  )

  ## draw a sampling matrix R
  mySampler$rSamplingMatrix()

  ## turn this matrix to a sampled Network object
  adjacencyMatrix[mySampler$samplingMatrix == 0] <- NA
  sampledNet <- sampledNetwork$new(adjacencyMatrix)
  sampledNet
}
