#' @title Simulation of a Stochastic Block Model
#'
#' @description \code{simulateSBM} is a function that generates a matrix (the adjacency matrix of a network) under the SBM
#'
#' @param n The number of nodes
#' @param alpha The mixture parameters
#' @param pi The connectivity matrix (probabilities inter and intra clusters)
#' @param family The emission law of the adjacency matrix : Bernoulli or Poisson
#' @param directed Boolean variable to indicate whether the network is directed or not,
#' by default "undirected" is choosen
#' @return \code{simulateSBM} returns a vector with clusters of nodes and a matrix (the adjacency matrix of the network)
#' @author T. Tabouy
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
#' alpha <- rep(1,Q)/Q                                                                # mixture parameter
#' pi <- diag(.45,Q) + .05                                                            # connectivity matrix
#' family <- "Bernoulli"                                                              # the emmission law
#' directed <- FALSE                                                                  # if the network is directed or not
#' mySBM <- simulateSBM(n, alpha, pi, family, directed)                               # simulation of ad Bernoulli non-directed SBM
#'
#'### Results : ###
#' clusters <-  mySBM$clusters                                                        # clusters
#' adjacencyMatrix <- mySBM$adjacencyMatrix                                           # the adjacency matrix
#'
#'
#' @export
simulateSBM <- function(n, alpha, pi, family="Bernoulli", directed=FALSE){

  mySBM <- SBM$new(family, directed, n, alpha, pi)
  mySBM$rBlocks()
  mySBM$rAdjMatrix()

  return(list(clusters = mySBM$clusters, adjacencyMatrix = mySBM$adjacencyMatrix))
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
#' @author T. Tabouy
#' @references [1] Tabouy, P. Barbillon, J. Chiquet. Variationnal inference of Stochastic Block Model from sampled data (2017). arXiv:1707.04141.
#' @seealso \code{\link{inferSBM}} and \code{\link{samplingSBM}}
#' @details The differents sampling designs are splitted into two families in which we find dyad-centered and node-centered sampling, for
#' more details see (\cite{1}) :\itemize{\item Missing At Random (MAR) \itemize{\item{MAREdge: parameter = p
#' \deqn{p = P(Dyad (i,j) is sampled)}}
#' \item{MARnode: parameter = p and
#' \deqn{p = P(Node i is sampled)}}
#' \item{snowball (one step):
#' like the MARNode sampling plus we sample neighbours of nodes sampled at the first batch}
#' }
#' \item Not Missing At Random (NMAR) \itemize{ \item{doubleStandard: parameter = (p0,p1) and
#' \deqn{p0 = P(Dyad (i,j) is sampled | the dyad is equal to 0)=}, p1 = P(Dyad (i,j) is sampled | the dyad is equal to 1)}
#' \item{starDegree: parameter = c(a,b) and
#' \deqn{logit(a+b*Degree(i)) = P(Node i is sampled | Degree(i))}}
#' \item{class: parameter = c(p(1),...,p(Q)) and
#' \deqn{p(q) = P(Node i is sampled | node i is in cluster q)}}
#' }}
#' @examples
#' ### A SBM model : ###
#' n <- 300
#' Q <- 3
#' alpha <- rep(1,Q)/Q                                                                # mixture parameter
#' pi <- diag(.45,Q) + .05                                                            # connectivity matrix
#' family <- "Bernoulli"                                                              # the emmission law
#' directed <- FALSE                                                                  # if the network is directed or not
#' mySBM <- simulateSBM(n, alpha, pi, family, directed)                               # simulation of ad Bernoulli non-directed SBM
#'
#'### Results : ###
#' adjacencyMatrix <- mySBM$adjacencyMatrix                                           # the adjacency matrix
#'
#' ## Sampling of the data : ##
#' samplingParameters <- .5                                                           # the sampling rate
#' sampling <- "MAREdge"                                                              # the sampling design
#' sampledAdjMatrix <- samplingSBM(adjacencyMatrix, sampling, samplingParameters)     # the sampled adjacency matrix
#'
#'
#' @export
samplingSBM <- function(adjacencyMatrix, sampling, parameters, clusters = NULL){

  ## TODO: postponed all the checks to the Class sampling_model
  family <- ifelse(length(tabulate(adjacencyMatrix)) == 1, "Bernoulli", "Poisson")

  if (!(sampling %in% available_samplings)) stop("This sampling is not in the list !")
  if (!(family == "Bernoulli" | sampling %in% c("edge", "node"))) stop("This sampling for Poisson emission law is not available !")

  N <- nrow(adjacencyMatrix)

  if (sampling == "block"){
    if (is.null(clusters)) stop("For class sampling you must give clusters !")
    Q <- nlevels(as.factor(clusters))
    if (!(length(clusters) == N)) stop(paste("The parameter clusters must have a length equal to", N,"!"))
  }

  if (!switch(sampling,
        "double_standard" = ifelse(length(parameters) == 2, TRUE, FALSE),
        "degree"          = ifelse(length(parameters) == 2, TRUE, FALSE),
        "block"           = ifelse(length(parameters) == Q, TRUE, FALSE),
        "dyad"            = ifelse(length(parameters) == 1, TRUE, FALSE),
        "node"            = ifelse(length(parameters) == 1, TRUE, FALSE),
        "snowball"        = ifelse(length(parameters) == 1, TRUE, FALSE))) {
    stop("Sampling parameters have not good length")
  }

  if (sampling != "degree") {
    if(any(parameters < 0) | any(parameters > 1)){
      stop("Sampling parameters must be probabilities (i.e between 0 and 1)")
    }
  }

  mySampling <- sampling_model$new(sampling, parameters)

  if(sampling == "block"){
    mySampling$rSampling(adjacencyMatrix, clusters)
  } else {
    mySampling$rSampling(adjacencyMatrix)
  }
  mySampling$sampledNetwok$adjacencyMatrix
}

#' @title Inference of Stochastic Block Model from sampled data
#'
#' @description \code{inferSBM} is a function that makes variationnal inference of Stochastic Block Model from sampled adjacency matrix
#'
#' @param sampledNetwork The sampled network data (a square matrix)
#' @param vBlocks The vector of number of blocks considered in the collection
#' @param sampling The sampling design for missing data modeling : MAREdge, doubleStandard, MARNode, snowball, starDegree, class
#' by default "undirected" is choosen
#' @param plot Summary of the output of the algorithm, by default TRUE is choosen
#' @return \code{inferSBM} returns a list with the best model choosen following the ICL criterion, a list with all models estimated for all Q in vBlocks
#' and a vector with ICL calculated for all Q in vBlocks
#' @author T. Tabouy
#' @references [1] Tabouy, P. Barbillon, J. Chiquet. Variationnal inference of Stochastic Block Model from sampled data (2017). arXiv:1707.04141.
#' @seealso \code{\link{samplingSBM}} and \code{\link{simulateSBM}} and \code{\link{SBM_collection}}.
#' @examples
#' ### A SBM model : ###
#' n <- 300
#' Q <- 3
#' alpha <- rep(1,Q)/Q                                                                # mixture parameter
#' pi <- diag(.45,Q) + .05                                                            # connectivity matrix
#' family <- "Bernoulli"                                                              # the emmission law
#' directed <- FALSE                                                                  # if the network is directed or not
#' mySBM <- simulateSBM(n, alpha, pi, family, directed)                               # simulation of ad Bernoulli non-directed SBM
#'
#'### Results : ###
#' adjacencyMatrix <- mySBM$adjacencyMatrix                                           # the adjacency matrix
#'
#' ## Sampling of the data : ##
#' samplingParameters <- .5                                                           # the sampling rate
#' sampling <- "MAREdge"                                                              # the sampling design
#' sampledAdjMatrix <- samplingSBM(adjacencyMatrix, sampling, samplingParameters)     # the sampled adjacency matrix
#'
#' ## Inference :
#' vBlocks <- 1:5                                                                     # number of classes
#' sbm <- inferSBM(sampledAdjMatrix, vBlocks, sampling, family, directed)             # the inference
#'
#' @import R6
#' @import igraph
#' @export
inferSBM <- function(sampledNetwork, vBlocks, sampling, plot = TRUE){

  if (isSymmetric(sampledNetwork)) directed <- FALSE else directed <- TRUE
  if (length(table(sampledNetwork)) > 2) family <- "Poisson" else family <- "Bernoulli"

  if (!(sampling %in% available_samplings)) stop("This sampling is not in the list !")
  if (!(family == "Bernoulli" | sampling %in% c("MAREdge", "MARNode"))) stop("This sampling for Poisson emission law is not available !")
  if (is.null(vBlocks)) stop(" The parameter vBlocks must be a least of length 1 !")

  collection <- SBM_collection$new(sampledNetwork, vBlocks, sampling, family, directed)



  collectionList <- lapply(collection$models, function(x){
    return(list(Q = x$SBM$nBlocks, alpha = x$SBM$mixtureParam, pi = x$SBM$connectParam, clusters = apply(x$blockVarParam, 1, which.max), samplingParameters = x$sampling$missingParam))
  })

  if (length(vBlocks) > 1 & min(vBlocks) == 1) {
    bestModel <- collectionList[[collection$getBestModel()]]
  }
  if (min(vBlocks) > 1 & length(vBlocks) != 1) {
    bestModel <- collectionList[[collection$getBestModel() - min(vBlocks) + 1]]
  }
  if (length(vBlocks) == 1) {
    bestModel <- collectionList[[1]]
  }

  if (plot) {
    mode <- ifelse(directed, "directed", "undirected")
    sampAdjMatZeros <- sampledNetwork; sampAdjMatZeros[is.na(sampAdjMatZeros)] <- 0

    G1 <- graph_from_adjacency_matrix(bestModel$pi, mode = mode, weighted = TRUE, diag = TRUE)
    G2 <- graph_from_adjacency_matrix(sampAdjMatZeros, mode = mode, weighted = TRUE, diag = TRUE)
    cl <- bestModel$clusters
    par(mfrow=c(2,2))
    plot.igraph(G1,vertex.size=table(cl),edge.width=E(G1)$weight*10, main="connectivity matrix", vertex.color=1:20)
    image.NA(sampledNetwork[order(cl), order(cl)], axes=FALSE, na.color = "red", main="clustered network + NA")
    plot(vBlocks,collection$vICLs, type = 'l', main = "Integrated complete Likekihood (ICL)", xlab = "Q", ylab = "ICL")
    plot.igraph(G2, vertex.color = cl, main="Network structure")
  }

  return(list(bestModel = bestModel, models = collectionList, ICL = collection$vICLs))
}

image.NA <- function(z,  zlim=c(0,1), col=c("white", "midnightblue"), na.color='red', outside.below.color='black', outside.above.color='white',...)
{
  zstep <- (zlim[2] - zlim[1]) / length(col); # step in the color palette
  newz.below.outside <- zlim[1] - 2 * zstep # new z for values below zlim
  newz.above.outside <- zlim[2] + zstep # new z for values above zlim
  newz.na <- zlim[2] + 2 * zstep # new z for NA

  z[which(z<zlim[1])] <- newz.below.outside # we affect newz.below.outside
  z[which(z>zlim[2])] <- newz.above.outside # we affect newz.above.outside
  z[which(is.na(z>zlim[2]))] <- newz.na # same for newz.na

  zlim[1] <- zlim[1] - 2 * zstep # extend lower limit to include below value
  zlim[2] <- zlim[2] + 2 * zstep # extend top limit to include the two new values above and na

  col <- c(outside.below.color, col[1], col, outside.above.color, na.color) #correct by including col[1] at bottom of range

  par(mar=c(2.1,8.1,3.1,3.1))
  image(z[nrow(z):1,],  zlim=zlim, col=col, xaxt="n", yaxt="n", main="clustered network + NA") # we finally call image(...)
  box()
  # par(mar=c(5.1,4.1,4.1,2.1))
}




