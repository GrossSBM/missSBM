#' @title Simulation of a Stochastic Block Model
#'
#' @description \code{simulateSBM} is a function that generates a matrix (the adjacency matrix of a network) under the SBM
#'
#' @param n The number of nodes
#' @param alpha The mixture parameters
#' @param pi The connectivity matrix (probabilities inter and intra clusters)
#' @param directed Boolean variable to indicate whether the network is directed or not,
#' by default "undirected" is choosen
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
#' n <- 300
#' Q <- 3
#' alpha <- rep(1,Q)/Q                                      # mixture parameter
#' pi <- diag(.45,Q) + .05                                  # connectivity matrix
#' directed <- FALSE
#' mySBM <- simulateSBM(n, alpha, pi, directed)             # simulation of ad Bernoulli non-directed SBM
#'
#'### Results : ###
#' clusters <-  mySBM$clusters                              # clusters
#' adjacencyMatrix <- mySBM$adjacencyMatrix                 # the adjacency matrix
#'
#' @export
simulateSBM <- function(n, alpha, pi, directed=FALSE){
  mySBM <- SBM_sampler$new(directed, n, alpha, pi)
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
#' n <- 300
#' Q <- 3
#' alpha <- rep(1,Q)/Q                                  # mixture parameter
#' pi <- diag(.45,Q) + .05                              # connectivity matrix
#' directed <- FALSE                                    # if the network is directed or not
#' mySBM <- simulateSBM(n, alpha, pi, directed)         # simulation of ad Bernoulli non-directed SBM
#'
#'### Results : ###
#' adjacencyMatrix <- mySBM$adjacencyMatrix             # the adjacency matrix
#'
#' ## Sampling of the data : ##
#' samplingParameters <- .5                             # the sampling rate
#' sampling <- "dyad"                                   # the sampling design
#' sampledNetwork <- samplingSBM(adjacencyMatrix, sampling, samplingParameters)
#'
#' @export
samplingSBM <- function(adjacencyMatrix, sampling, parameters, clusters = NULL){

  if (!(sampling %in% available_samplings))
    stop("This sampling is not available!")

  if (sampling == "block" & is.null(clusters))
    stop("For class sampling you must give clusters !")

  mySampling <- networkSampling_sampler$new(sampling, parameters)
  mySampling$rSampling(adjacencyMatrix, clusters)
}

#' @title Inference of Stochastic Block Model from sampled data
#'
#' @description \code{inferSBM} is a function that makes variationnal inference of Stochastic Block Model from sampled adjacency matrix
#'
#' @param sampledNetwork The sampled network data (a square matrix)
#' @param vBlocks The vector of number of blocks considered in the collection
#' @param sampling The sampling design for missing data modeling : "dyad", "double_standard", "node", "snowball", "degree", "block"
#' by default "undirected" is choosen
#' @param plot Summary of the output of the algorithm, by default TRUE is choosen
#' @return \code{inferSBM} returns a list with the best model choosen following the ICL criterion, a list with all models estimated for all Q in vBlocks
#' and a vector with ICL calculated for all Q in vBlocks
#' @references [1] Tabouy, P. Barbillon, J. Chiquet. Variationnal inference of Stochastic Block Model from sampled data (2017). arXiv:1707.04141.
#' @seealso \code{\link{samplingSBM}} and \code{\link{simulateSBM}} and \code{\link{SBM_collection}}.
#' @examples
#' ### A SBM model : ###
#' n <- 300
#' Q <- 3
#' alpha <- rep(1,Q)/Q                                      # mixture parameter
#' pi <- diag(.45,Q) + .05                                  # connectivity matrix
#' directed <- FALSE
#' mySBM <- simulateSBM(n, alpha, pi, directed)             # simulation of ad Bernoulli non-directed SBM
#'
#'### Results : ###
#' adjacencyMatrix <- mySBM$adjacencyMatrix                                           # the adjacency matrix
#'
#' ## Sampling of the data : ##
#' samplingParameters <- .5                                                           # the sampling rate
#' sampling <- "dyad"                                                                 # the sampling design
#' sampledAdjMatrix <- samplingSBM(adjacencyMatrix, sampling, samplingParameters)     # the sampled adjacency matrix
#'
#' ## Inference :
#' vBlocks <- 1:5                                                                     # number of classes
#' sbm <- inferSBM(sampledAdjMatrix$adjacencyMatrix, vBlocks, sampling)
#'
#' @import R6 parallel
#' @include smoother_SBM.R
#' @export
inferSBM <- function(adjacencyMatrix, vBlocks, sampling, clusterInit = "spectral",
                     smoothing = c("none", "forward", "backward", "both"), mc.cores = 2, control_VEM = list()){

  sampledNet <- sampledNetwork$new(adjacencyMatrix)
  cat("\n")
  cat("\n Adjusting Variational EM for Stochastic Block Model\n")
  cat("\n\tImputation assumes a '", sampling,"' network-sampling process\n", sep = "")
  cat("\n")
  models <- lapply(vBlocks,
    function(nBlocks) {
    cat(" Initialization of model with", nBlocks,"blocks.", "\r")
      missingSBM_fit$new(sampledNet, nBlocks, sampling, clusterInit)
    }
  )

  ## defaut control parameter for VEM, overwritten by user specification
  control <- list(threshold = 1e-4, maxIter = 200, fixPointIter = 5, trace = FALSE)
  control[names(control_VEM)] <- control_VEM
  cat("\n")
  res_optim <- do.call(rbind, lapply(models,
    function(model) {
      cat(" Performing VEM inference for model with", model$fittedSBM$nBlocks,"blocks.\r")
      res <- model$doVEM(control)
      res$nBlocks <- model$fittedSBM$nBlocks
      res$iteration <- 1:nrow(res)
      res
    }
  ))

  smoothing <- match.arg(smoothing)
  if (smoothing != "none") {
    cat("\n Smoothing ICL\n")
    smoothing_fn <- switch(smoothing,
      "forward"  = smoothingForward_SpCl  ,
      "backward" = smoothingBackward      ,
      "both"     = smoothingForBackWard_SpCl
    )
    control_VEM$trace <- FALSE # forcing no trace while smoothing
    models <- smoothing_fn(models, vBlocks, sampledNet, sampling, mc.cores, control_VEM)
  }

  return(list(models = models, monitor = res_optim))
}

