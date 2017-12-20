#' @title Inference of Stochastic Block Model from sampled data
#'
#' @description \code{missingSBM} is a function that makes variationnal inference of Stochastic Block Model from sampled adjacency matrix
#'
#'
#' @param sampledNetwork The sampled network data (a square matrix)
#' @param vBlocks The vector of number of blocks considered in the collection
#' @param sampling The sampling design for missing data modeling : MAREdge, doubleStandard, MARNode, snowball, starDegree, class
#' @param family The emission law of the adjacency matrix : Bernoulli or Poisson
#' @param directed Boolean variable to indicate whether the network is directed or not,
#' by default "undirected" is choosen
#' @return \code{missingSBM} returns a \code{\link{SBM_collection}} object.
#' @author T. Tabouy
#' @references [1] Tabouy et al., Variationnal inference of Stochastic Block Model from sampled data (2017). arXiv:1707.04141.
#' @seealso \code{\link{SBM_collection}}.
#' @examples
#' ## A SBM model : ##
#' n <- 300
#' Q <- 3
#' alpha <- rep(1,Q)/Q                                                                # mixture parameter
#' pi <- diag(.45,Q) + .05                                                            # connectivity matrix
#' family <- "Bernoulli"                                                              # the emmission law
#' directed <- FALSE                                                                  # if the network is directed or not
#' mySBM <- SBM_BernoulliUndirected$new(n, alpha, pi)                                 # the model object
#' SBMdata <- mySBM$rSBM()                                                            # simulation of the complete data
#'
#' ## Sampling of the data : ##
#' sampling_rate <- .5                                                                # the sampling rate
#' sampling <- "MAREdge"                                                              # the sampling design
#' mySampled  <- sampling_randomPairMAR$new(n, sampling_rate, directed)               # the sampling object
#' sample     <- mySampled$rSampling(SBMdata$adjacencyMatrix)                         # simulation of a sampling matrix
#'
#' ## Inference :
#' sampledNetwork <- sample$adjacencyMatrix                                           # the adjacency matrix
#' vBlocks <- 1:5                                                                     # number of classes
#' sbm <- missingSBM(sample$adjacencyMatrix, vBlocks, sampling, family, directed)     # the inference
#'
#' ## Results/Estimation :
#' plot(sbm$vICLs)                                                                    # the ICL criterion to select Q
#' sbm$getBestModel()                                                                 # estimation of Q
#' apply(sbm$models[[Q]]$blockVarParam, 1, which.max)                                 # clustering
#' sbm$models[[Q]]$SBM$connectParam                                                   # connectivity matrix
#' sbm$models[[Q]]$SBM$mixtureParam                                                   # mixture parameters
#'
#' @export
missingSBM <- function(sampledNetwork, vBlocks, sampling, family, directed){

  library(R6)
  library(parallel)
  library(mclust)
  library(igraph)

  if(!directed){
    if(!isSymmetric(sampledNetwork)){
      stop("The adjacency matrix is not symmetric !")
    }
  }
  if(family == "Poisson"){
    if(!(sampling %in% c("MAREdge", "MARNode"))){
      stop("This sampling for Poisson emission law is not available")
    }
  }
  return(SBM_collection$new(sampledNetwork, vBlocks, sampling, family, directed))
}
