#' an SBM model
#'
#' @field nNodes          number of nodes
#' @field nBlocks         number of blocks
#' @field blockProportion vector of block proportion (a.k.a. alpha)
#' @field modelParameters vector of model parameters (a.k.a. theta)
#'
#' @importFrom R6 R6Class
#' @export
SBM <-
R6Class(classname = "SBM",
  public = list(
    ## fields
    nNodes         = NULL, # number of nodes
    nBlocks        = NULL, # number of blocks
    mixtureParam   = NULL, # vector of block parameters (a.k.a. alpha)
    connectParam   = NULL, # vector of model parameters (a.k.a. theta)
    ## methods
    initialize = function(nNodes=NA, mixtureParam=NA, connectParam=NA) {
      self$nNodes       <- nNodes
      self$mixtureParam <- mixtureParam
      self$connectParam <- connectParam
      self$nBlocks      <- length(mixtureParam)
    },
    ## a function to generate a vector of clusters indicators
    rBlocks = function() {
      return(t(rmultinom(self$nNodes, size = 1, prob = self$mixtureParam)))
    },
    ## a function to generate a matrix of block indicators
    rSBM = function() {
      return(list(blocks = self$rBlocks(), adjacencyMatrix = NA))
    }
  )
)

#'
#' @field completeLogLik  a function to compute the completed log-likelihood
#'
#' @export
SBM_BernoulliUndirected <-
R6Class(classname = "SBM_BernoulliUndirected",
  inherit = SBM,
  public = list(
    initialize = function(nNodes=NA, mixtureParam=NA, connectParam=NA) {
      super$initialize(nNodes, mixtureParam, connectParam)
    },
    completeLogLik = function(completedNetwork, blockIndicators, sampledNetwork) {
      if(missing(sampledNetwork)){
        network     <- completedNetwork
        network.bar <- 1 - network ; diag(network.bar) <- 0
      } else if(missing(completedNetwork)){
        network     <- sampledNetwork$adjacencyMatrix * sampledNetwork$samplingMatrix
        network.bar <- (1 - network) * sampledNetwork$samplingMatrix ; diag(network.bar) <- 0
      }
      return(sum(blockIndicators %*% log(self$mixtureParam)) +
        .5 * sum( network.bar *(blockIndicators %*% log(self$connectParam) %*% t(blockIndicators)) +
                    network * (blockIndicators %*% log(1-self$connectParam) %*% t(blockIndicators))))
    },
    rSBM = function() {
      blocks <- super$rSBM()$blocks
      adjacencyMatrix <- matrix(rbinom(self$nNodes^2,1, blocks %*% self$connectParam %*% t(blocks)),self$nNodes)
      adjacencyMatrix <- adjacencyMatrix * lower.tri(adjacencyMatrix) + t(adjacencyMatrix * lower.tri(adjacencyMatrix))
      diag(adjacencyMatrix) <- 0
      return(list(blocks = blocks, adjacencyMatrix = adjacencyMatrix))
    }
    }
  )
)

#' @export
SBM_BernoulliUndirected.fit <-
  R6Class(classname = "SBM_BernoulliUndirected.fit",
          inherit = SBM_BernoulliUndirected,
          public = list(
            initialize = function(nNodes=NA, mixtureParam=NA, connectParam=NA) {
              super$initialize(nNodes, mixtureParam, connectParam)
            },
            maximization = function(SBM, completedNetwork, blockVarParam) {
              SBM$connectParam <- (t(blockVarParam)%*% completedNetwork %*%blockVarParam) / (t(blockVarParam)%*%((1-diag(self$nNodes)))%*%blockVarParam)
              SBM$mixtureParam <-  colMeans(blockVarParam)
              return(SBM)
            },
            maximization_MAR = function(SBM, completedNetwork, blockVarParam, samplingMatrix) {
              SBM$connectParam <- (t(blockVarParam)%*% completedNetwork %*%blockVarParam) / (t(blockVarParam)%*%((1-diag(self$nNodes))*samplingMatrix)%*%blockVarParam)
              SBM$mixtureParam <- colMeans(blockVarParam)
              return(SBM)
            }, 
            fixpoint = function(SBM, blockVarParam, completedNetwork) {
              completedNetwork.bar <- 1 - completedNetwork; diag(completedNetwork.bar) <- 0
              blockVarParam_new <- exp(sweep(completedNetwork %*% blockVarParam %*% t(log(SBM$connectParam)) + 
                                               completedNetwork.bar %*% blockVarParam %*% t(log(1-SBM$connectParam)),2,log(SBM$mixtureParam),"+"))
              num               <- rowSums(blockVarParam_new)
              blockVarParam_new <- blockVarParam_new/num
            }
          )
  )



#' @export
SBM_BernoulliDirected <-
  R6Class(classname = "SBM_BernoulliDirected",
          inherit = SBM,
          public = list(
            initialize = function(nNodes=NA, mixtureParam=NA, connectParam=NA) {
              super$initialize(nNodes, mixtureParam, connectParam)
            },
            completeLogLik = function(completedNetwork, blockIndicators, sampledNetwork) {
              if(missing(sampledNetwork)){
                network     <- completedNetwork
                network.bar <- 1 - network ; diag(network.bar) <- 0
              } else if(missing(completedNetwork)){
                network     <- sampledNetwork$adjacencyMatrix * sampledNetwork$samplingMatrix
                network.bar <- (1 - network) * sampledNetwork$samplingMatrix ; diag(network.bar) <- 0
              }
              return(sum(blockIndicators %*% log(self$mixtureParam)) +
                sum( network.bar *(blockIndicators %*% log(self$connectParam) %*% t(blockIndicators)) +
                            network * (blockIndicators %*% log(1-self$connectParam) %*% t(blockIndicators))))
            },
            rSBM = function() {
              blocks <- super$rSBM()$blocks
              adjacencyMatrix <- matrix(rbinom(self$nNodes^2,1, blocks %*% self$connectParam %*% t(blocks)),self$nNodes)
              diag(adjacencyMatrix) <- 0
              return(list(blocks = blocks, adjacencyMatrix = adjacencyMatrix))
            }
          )
  )

# mySBM <- SBM_BernoulliDirected$new(20, c(1/2, 1/2), matrix(c(.2, .05, .05, .2),2,2))
# mySBM$rBlocks()
# image(mySBM$rSBM()$adjacencyMatrix)

### Pour le cas MAR :
# sum(blockIndicators%*%log(mixtureParam)) + sum(completedNetwork*(blockIndicators%*%log(connectParam)%*%t(blockIndicators))) -
#   sum(R*log(factorial(completedNetwork))*(blockIndicators%*%matrix(1,Q,Q)%*%t(blockIndicators))) -
#   sum((R*L)*blockIndicators%*%connectParam%*%t(blockIndicators))


#' @export
SBM_PoissonDirected <-
  R6Class(classname = "SBM_PoissonDirected",
          inherit = SBM,
          public = list(
            initialize = function(nNodes=NA, mixtureParam=NA, connectParam=NA) {
              super$initialize(nNodes, mixtureParam, connectParam)
            },
            completeLogLik = function(completedNetwork, blockIndicators, sampledNetwork) {
              if (missing(sampledNetwork)) {
                loop.bar    <- matrix(1,self$nNodes,self$nNodes) ; diag(loop.bar) <- 0
                return(sum(blockIndicators%*%log(self$mixtureParam)) + 
                  sum(completedNetwork*(blockIndicators%*%log(self$connectParam)%*%t(blockIndicators))) - 
                  sum(log(factorial(completedNetwork))*(blockIndicators%*%matrix(1,Q,Q)%*%t(blockIndicators))) -
                  sum(loop.bar*blockIndicators%*%connectParam%*%t(blockIndicators)))
              } else if (missing(completedNetwork)) {
                network     <- sampledNetwork$adjacencyMatrix * sampledNetwork$samplingMatrix
                loop.bar    <- matrix(1,self$nNodes,self$nNodes) ; diag(loop.bar) <- 0
                return(sum(blockIndicators%*%log(self$mixtureParam)) + 
                  sum(network*(blockIndicators%*%log(self$connectParam)%*%t(blockIndicators))) - 
                  sum(sampledNetwork$samplingMatrix * log(factorial(network))*(blockIndicators%*%matrix(1,Q,Q)%*%t(blockIndicators))) -
                  sum((sampledNetwork$samplingMatrix * loop.bar)*blockIndicators%*%connectParam%*%t(blockIndicators)))
              }
            },
            rSBM = function() {
              blocks <- super$rSBM()$blocks
              adjacencyMatrix <- matrix(rpois(self$nNodes^2, blocks %*% self$connectParam %*% t(blocks)) ,self$nNodes)
              diag(adjacencyMatrix) <- 0
              return(list(blocks = blocks, adjacencyMatrix = adjacencyMatrix))
            }
          )
  )

# mySBM <- SBM_PoissonDirected$new(20, c(1/2, 1/2), matrix(c(1,2,3,4),2,2))
# mySBM$rBlocks()
# image(mySBM$rSBM()$adjacencyMatrix)

#' @export
SBM_PoissonUndirected <-
  R6Class(classname = "SBM_PoissonUndirected",
          inherit = SBM,
          public = list(
            initialize = function(nNodes=NA, mixtureParam=NA, connectParam=NA) {
              super$initialize(nNodes, mixtureParam, connectParam)
            },
            completeLogLik = function(completedNetwork, blockIndicators) {
              if (missing(sampledNetwork)) {
                loop.bar    <- matrix(1,self$nNodes,self$nNodes) ; diag(loop.bar) <- 0
                return(sum(blockIndicators%*%log(self$mixtureParam)) + 
                         .5 * (sum(completedNetwork*(blockIndicators%*%log(self$connectParam)%*%t(blockIndicators))) - 
                         sum(log(factorial(completedNetwork))*(blockIndicators%*%matrix(1,Q,Q)%*%t(blockIndicators))) -
                         sum(loop.bar*blockIndicators%*%connectParam%*%t(blockIndicators))))
              } else if (missing(completedNetwork)) {
                network     <- sampledNetwork$adjacencyMatrix * sampledNetwork$samplingMatrix
                loop.bar    <- matrix(1,self$nNodes,self$nNodes) ; diag(loop.bar) <- 0
                return(sum(blockIndicators%*%log(self$mixtureParam)) + 
                         .5 * (sum(network*(blockIndicators%*%log(self$connectParam)%*%t(blockIndicators))) - 
                         sum(sampledNetwork$samplingMatrix * log(factorial(network))*(blockIndicators%*%matrix(1,Q,Q)%*%t(blockIndicators))) -
                         sum((sampledNetwork$samplingMatrix * loop.bar)*blockIndicators%*%connectParam%*%t(blockIndicators))))
              }
            },
            rSBM = function() {
              blocks <- super$rSBM()$blocks
              adjacencyMatrix <- matrix(rpois(self$nNodes^2, blocks %*% self$connectParam %*% t(blocks)) ,self$nNodes)
              adjacencyMatrix <- adjacencyMatrix * lower.tri(adjacencyMatrix) + t(adjacencyMatrix * lower.tri(adjacencyMatrix))
              diag(adjacencyMatrix) <- 0
              return(list(blocks = blocks, adjacencyMatrix = adjacencyMatrix))
            }
          )
  )

# mySBM <- SBM_PoissonUndirected$new(20, c(1/2, 1/2), matrix(c(1,2,3,4),2,2))
# mySBM$rBlocks()
# image(mySBM$rSBM()$adjacencyMatrix)


