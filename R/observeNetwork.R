#' Observe a network partially according to a given sampling design
#'
#' This function samples observations in an adjacency matrix according to a given network sampling design.
#'
#' @param adjacencyMatrix The N x N adjacency matrix of the network to sample.
#' @param sampling The sampling design used to observe the adjacency matrix, see details.
#' @param parameters The sampling parameters (adapted to each sampling, see details).
#' @param clusters An optional clustering membership vector of the nodes. Only necessary for block samplings.
#' @param covariates An optional list with M entries (the M covariates). If the covariates are node-centered,
#'     each entry of \code{covariates}. must be a size-N vector;  if the covariates are dyad-centered, each entry
#'     of \code{covariates} must be N x N matrix.
#' @param similarity An optional function to compute similarities between node covariates. Default is
#' \code{missSBM:::l1_similarity}, that is, -abs(x-y). Only relevant when the covariates are node-centered.
#' @param intercept An optional intercept term to be added in case of the presence of covariates. Default is 0.
#'
#' @inherit estimateMissSBM details
#'
#' @return an adjacency matrix with the same dimension as the input, yet with additional NAs.
#'
#' @examples
#' ## SBM parameters
#' N <- 300 # number of nodes
#' Q <- 3   # number of clusters
#' pi <- rep(1,Q)/Q     # block proportion
#' theta <- list(mean = diag(.45,Q) + .05 ) # connectivity matrix
#'
#' ## simulate an unidrected binary SBM without covariate
#' sbm <- sbm::sampleSimpleSBM(N, pi, theta)
#'
#' ## Sample network data
#'
#' # some sampling design and their associated parameters
#' sampling_parameters <- list(
#'    "dyad" = .3,
#'    "node" = .3,
#'    "double-standard" = c(0.4, 0.8),
#'    "block-node" = c(.3, .8, .5),
#'    "block-dyad" = theta$mean,
#'    "degree" = c(.01, .01),
#'    "snowball" = c(2,.1)
#'  )
#'
#' sampled_networks <- list()
#'
#' for (sampling in names(sampling_parameters)) {
#'   sampled_networks[[sampling]] <-
#'      missSBM::observeNetwork(
#'        adjacencyMatrix = sbm$netMatrix,
#'        sampling        = sampling,
#'        parameters      = sampling_parameters[[sampling]],
#'        cluster         = sbm$memberships
#'      )
#' }
#' @export
observeNetwork <- function(adjacencyMatrix, sampling, parameters, clusters = NULL, covariates = NULL, similarity = l1_similarity, intercept = 0) {

### TEMPORARY FIX
  adjacencyMatrix[is.na(adjacencyMatrix)] <- 0

  ## Sanity check
  stopifnot(sampling %in% available_samplings)

  ## general network parameters
  nbNodes   <- ncol(adjacencyMatrix)
  directed <- !isSymmetric(adjacencyMatrix)

  ## Prepare the covariates
  covar <- format_covariates(covariates, similarity)
  if (!is.null(covar$Array)) stopifnot(sampling %in% available_samplings_covariates)

  ## instantiate the sampler
  mySampler <-
    switch(sampling,
      "dyad"       = simpleDyadSampler$new(
        parameters = parameters, nbNodes = nbNodes, directed = directed),
      "node"       = simpleNodeSampler$new(
        parameters = parameters, nbNodes = nbNodes, directed = directed),
      "covar-dyad" = simpleDyadSampler$new(
        parameters = parameters, nbNodes = nbNodes, directed = directed, covarArray  = covar$Array, intercept = intercept),
      "covar-node" = simpleNodeSampler$new(
        parameters = parameters, nbNodes = nbNodes, directed = directed, covarMatrix = covar$Matrix, intercept = intercept),
      "double-standard" = doubleStandardSampler$new(
        parameters = parameters, adjMatrix = adjacencyMatrix, directed = directed),
      "block-dyad" = blockDyadSampler$new(
        parameters = parameters, nbNodes = nbNodes, directed = directed, clusters = clusters),
      "block-node" = blockNodeSampler$new(
        parameters = parameters, nbNodes = nbNodes, directed = directed, clusters = clusters),
      "degree"     = degreeSampler$new(
        parameters = parameters, degrees = rowSums(adjacencyMatrix), directed = directed),
      "snowball" = snowballSampler$new(
        parameters = parameters, adjacencyMatrix = adjacencyMatrix ,directed=directed)
  )

  ## draw a sampling matrix R
  mySampler$rSamplingMatrix()

  ## turn this matrix to a sampled Network object
  adjacencyMatrix[mySampler$samplingMatrix == 0] <- NA
  adjacencyMatrix
}
