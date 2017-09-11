#' an SBM model
#'
#' @field nNodes          number of nodes
#' @field nBlocks         number of blocks
#' @field blockProportion vector of block proportion (a.k.a. alpha)
#' @field modelParameters vector of model parameters (a.k.a. theta)
#'
#' @importFrom R6 R6Class
#' @export
sampling <-
  R6Class(classname = "sampling",
          public = list(
            ## fields
            nNodes         = NULL, # number of nodes
            missingParam   = NULL, # vector of missing parameters (a.k.a. alpha)
            directed       = FALSE, #
            ## methods
            initialize = function(nNodes=NA, missingParam=NA, directed = FALSE) {
              self$nNodes         <- nNodes
              self$missingParam   <- missingParam
              self$directed       <- directed
            }
          )
  )

#' @export
sampling_doubleStandard <-
  R6Class(classname = "sampling_doublestandard",
          inherit = sampling,
          public = list(
            initialize = function(nNodes, missingParam, directed = FALSE) {
              super$initialize(nNodes, missingParam, directed)
            },
            rSampling = function(adjMatrix) {
              samplingMatrix <- matrix(0, self$nNodes, self$nNodes)
              if(!self$directed){
                areOne  <- (adjMatrix == 1) & upper.tri(adjMatrix)
                areZero <- (adjMatrix == 0) & upper.tri(adjMatrix)
                
                samplingMatrix[areOne]  <- runif(sum(areOne))  < self$missingParam[2]
                samplingMatrix[areZero] <- runif(sum(areZero)) < self$missingParam[1]
                samplingMatrix <- t(samplingMatrix) | samplingMatrix
                diag(samplingMatrix) <- 1
              } else {
                areOne  <- (adjMatrix == 1)
                areZero <- (adjMatrix == 0)
                
                samplingMatrix[areOne]  <- runif(sum(areOne))  < self$missingParam[2]
                samplingMatrix[areZero] <- runif(sum(areZero)) < self$missingParam[1]
                diag(samplingMatrix) <- 1
              }
              
              sampAdjMatrix  <- adjMatrix
              sampAdjMatrix[which(samplingMatrix == 0)] <- NA
              return(sampledNetwork$new(sampAdjMatrix, self$directed))
            },
            samplingLogLik = function(sampledNetwork) {
              obsEdges  <- which(!is.na(sampledNetwork$adjacencyMatrix) & (upper.tri(sampledNetwork$adjacencyMatrix) | lower.tri(sampledNetwork$adjacencyMatrix)), arr.ind = TRUE)
              missEdges <- which(is.na(sampledNetwork$adjacencyMatrix) , arr.ind = TRUE)
              ll       <- sum(log(self$missingParam[2]) * sampledNetwork$adjacencyMatrix[obsEdges] + log(self$missingParam[1]) * (1-sampledNetwork$adjacencyMatrix[obsEdges])) +
                        sum(log(1-self$missingParam[2]) * sampledNetwork$adjacencyMatrix[missEdges] + log(1-self$missingParam[1]) * (1-sampledNetwork$adjacencyMatrix[missEdges]))
              if(!self$directed){
                return(ll/2)
              } else {
                return(ll)
              }
            }
          )
  )

# Undirected :
# mySBM <- sampling_doubleStandard$new(10, c(1/4, 1/2))
# Y <- matrix(round(runif(100, 0,1)),10,10)
# samp <- mySBM$rSampling(Y)
# ll <- mySBM$samplingLogLik(Y)

# # Directed :
# mySBM <- sampling_doubleStandard$new(10, c(1/4, 1/2), directed = TRUE)
# samp <- mySBM$rSampling(matrix(1,10,10))
# ll <- mySBM$samplingLogLik(matrix(1,10,10))mySBM <- sampling_doubleStandard$new(10, c(1/4, 1/2))
Y    <- matrix(round(runif(100, 0,1)),10,10)
samp <- mySBM$rSampling(Y)
ll   <- mySBM$samplingLogLik(Y)


#' @export
sampling_class <-
  R6Class(classname = "sampling_class",
          inherit = sampling,
          public = list(
            initialize = function(nNodes, missingParam, directed = FALSE) {
              super$initialize(nNodes, missingParam, directed)
            },
            rSampling = function(adjMatrix, blockIndicators) {
              samplingMatrix <- matrix(0, self$nNodes, self$nNodes)
              sampProb <- runif(self$nNodes) < self$missingParam[apply(blockIndicators, 1, which.max)]
              obsNodes <- which(runif(self$nNodes) < sampProb)
              
              samplingMatrix <- matrix(0,self$nNodes,self$nNodes) ; samplingMatrix[obsNodes,] <- 1
              diag(samplingMatrix) <- 1
              samplingMatrix <- (t(samplingMatrix) | samplingMatrix)*1
              
              sampAdjMatrix  <- adjMatrix
              sampAdjMatrix[which(samplingMatrix == 0)] <- NA
              return(sampledNetwork$new(sampAdjMatrix, self$directed))
            },
            samplingLogLik = function(sampledNetwork, blockIndicators) {
              samplingVector <- rep(0, self$nNodes); samplingVector[which(!is.na(rowSums(sampledNetwork$adjacencyMatrix)))] <- 1
              obsEdges       <- which(!is.na(sampledNetwork$adjacencyMatrix) & (upper.tri(sampledNetwork$adjacencyMatrix) | lower.tri(sampledNetwork$adjacencyMatrix)), arr.ind = TRUE)
              missEdges      <- which(is.na(sampledNetwork$adjacencyMatrix) , arr.ind = TRUE)
              return(sum(t(samplingVector) %*% blockIndicators %*% log(self$missingParam) + t(1-samplingVector) %*% blockIndicators %*% log(1-self$missingParam)))
            }
          )
  )

# mySBM <- sampling_class$new(10, c(1/4, 1/2))
# Y     <- matrix(round(runif(100, 0,1)),10,10)
# s     <- round(runif(10, 1,2)); Z <- matrix(0, 10, 2); Z[cbind(1:10, s)] <- 1
# samp  <- mySBM$rSampling(Y, Z)
# ll    <- mySBM$samplingLogLik(Y, Z)

#' @export
sampling_degree <-
  R6Class(classname = "sampling_degree",
          inherit = sampling,
          public = list(
            initialize = function(nNodes, missingParam, directed = FALSE) {
              super$initialize(nNodes, missingParam, directed)
            },
            rSampling = function(adjMatrix) {
              samplingMatrix <- matrix(0, self$nNodes, self$nNodes)
              sampProb            <- self$missingParam[1]+self$missingParam[2]*rowSums(AdjMatrix)
              samprob             <- 1/(1+exp(-sampProb))
              obsNodes            <- which(runif(self$nNodes) < sampProb)
              
              samplingMatrix <- matrix(0,self$nNodes,self$nNodes) ; samplingMatrix[obsNodes,] <- 1
              diag(samplingMatrix) <- 1
              samplingMatrix <- (t(samplingMatrix) | samplingMatrix)*1

              sampAdjMatrix  <- adjMatrix ; sampAdjMatrix[which(samplingMatrix == 0)] <- NA
              return(sampledNetwork$new(sampAdjMatrix, self$directed))
            },
            samplingLogLik = function(sampledNetwork) {
              samplingVector <- rep(0, self$nNodes); samplingVector[which(!is.na(rowSums(sampledNetwork$adjacencyMatrix)))] <- 1
              sampProb       <- self$missingParam[1]+self$missingParam[2]*rowSums(sampledNetwork$adjacencyMatrix)
              samprob        <- 1/(1+exp(-sampProb))
              return(log((sampProb^samplingVector)%*%((1-sampProb)^(1-samplingVector))) )
            }
          )
  )


# mySBM <- sampling_degree$new(10, runif(10))
# Y     <- matrix(round(runif(100, 0,1)),10,10)
# samp  <- mySBM$rSampling(Y)
# ll    <- mySBM$samplingLogLik(Y)


#' @export
sampling_randomPairMAR <-
  R6Class(classname = "sampling_randomPairMAR",
          inherit = sampling,
          public = list(
            initialize = function(nNodes, missingParam, directed = FALSE) {
              super$initialize(nNodes, missingParam, directed)
            },
            rSampling = function(adjMatrix) {
              samplingMatrix <- matrix(0, self$nNodes, self$nNodes)
              
              if(!self$directed){
                edgeSamp <- sample(which(lower.tri(adjMatrix)), floor((self$nNodes*(self$nNodes-1)/2)*self$missingParam))
              } else {
                edgeSamp <- sample(which(lower.tri(adjMatrix) | upper.tri(adjMatrix)), floor((self$nNodes^2 - self$nNodes)*self$missingParam))
              }
              
              samplingMatrix <- matrix(0,self$nNodes,self$nNodes)
              samplingMatrix[edgeSamp] <- 1
              samplingMatrix <- t(samplingMatrix) | samplingMatrix ; diag(samplingMatrix) <- 1
              
              sampAdjMatrix  <- adjMatrix ; sampAdjMatrix[which(samplingMatrix == 0)] <- NA
              return(sampledNetwork$new(sampAdjMatrix, self$directed))
            },
            samplingLogLik = function(sampledNetwork) {
              obsEdges  <- which(!is.na(sampledNetwork$adjacencyMatrix) & (upper.tri(sampledNetwork$adjacencyMatrix) | lower.tri(sampledNetwork$adjacencyMatrix)), arr.ind = TRUE)
              missEdges <- which(is.na(sampledNetwork$adjacencyMatrix) , arr.ind = TRUE)
              logPsi    <- ifelse (self$missingParam < .Machine$double.eps, 0, log(self$missingParam))
              log1mPsi  <- ifelse (self$missingParam > 1-.Machine$double.eps, 0, log(1-self$missingParam))
              ll        <- length(obsEdges)*logPsi + length(missEdges)*log1mPsi
              if(!self$directed){
                return(ll/2)
              } else {
                return(ll)
              }
            }
          )
  )

# mySBM <- sampling_randomPairMAR$new(10, 0.5, directed = TRUE)
# Y     <- matrix(round(runif(100, 0,1)),10,10)
# samp  <- mySBM$rSampling(Y)
# ll    <- mySBM$samplingLogLik(Y)


#' @export
sampling_randomNodesMAR <-
  R6Class(classname = "sampling_randomNodesMAR",
          inherit = sampling,
          public = list(
            initialize = function(nNodes, missingParam, directed = FALSE) {
              super$initialize(nNodes, missingParam, directed)
            },
            rSampling = function(adjMatrix) {
              samplingMatrix <- matrix(0, self$nNodes, self$nNodes)
              obsNodes            <- which(runif(self$nNodes) < rep(self$missingParam, self$nNodes))
              
              samplingMatrix <- matrix(0,self$nNodes,self$nNodes) ; samplingMatrix[obsNodes,] <- 1
              diag(samplingMatrix) <- 1
              samplingMatrix <- (t(samplingMatrix) | samplingMatrix)*1
              
              sampAdjMatrix  <- adjMatrix ; sampAdjMatrix[which(samplingMatrix == 0)] <- NA
              return(sampledNetwork$new(sampAdjMatrix, self$directed))
            },
            samplingLogLik = function(sampledNetwork) {
              samplingVector <- rep(0, self$nNodes); samplingVector[which(!is.na(rowSums(sampledNetwork$adjacencyMatrix)))] <- 1
              sampProb       <- rep(self$missingParam, self$nNodes)
              logPsi         <- ifelse (sampProb < .Machine$double.eps, 0, log(sampProb))
              log1mPsi       <- ifelse (sampProb > 1-.Machine$double.eps, 0, log(1-sampProb))
              return(logPsi*sum(samplingVector) + log1mPsi * sum(1-samplingVector))
            }
          )
  )


# mySBM <- sampling_randomNodesMAR$new(10, 0.5)
# Y     <- matrix(round(runif(100, 0,1)),10,10)
# samp  <- mySBM$rSampling(Y)
# ll    <- mySBM$samplingLogLik(Y)

#' @export
sampling_snowball <-
  R6Class(classname = "sampling_snowball",
          inherit = sampling,
          public = list(
            initialize = function(nNodes, missingParam, directed = FALSE) {
              super$initialize(nNodes, missingParam, directed)
            },
            rSampling = function(adjMatrix) {
              samplingMatrix <- matrix(0, self$nNodes, self$nNodes)
              obsNodes            <- which(runif(self$nNodes) < self$missingParam)
              
              samplingMatrix <- matrix(0,self$nNodes,self$nNodes) ; samplingMatrix[obsNodes,] <- 1
              diag(samplingMatrix) <- 1
              samplingMatrix <- (t(samplingMatrix) | samplingMatrix)*1
              
              sampAdjMatrix  <- adjMatrix ; sampAdjMatrix[which(samplingMatrix == 0)] <- NA
              return(sampledNetwork$new(sampAdjMatrix, self$directed))
            },
            samplingLogLik = function(sampledNetwork) {
              samplingVector <- rep(0, self$nNodes); samplingVector[which(!is.na(rowSums(sampledNetwork$adjacencyMatrix)))] <- 1
              return(log((self$missingParam^samplingVector)%*%((1-self$missingParam)^(1-samplingVector))))
            }
          )
  )

# mySBM <- sampling_snowball$new(10, runif(10))
# Y     <- matrix(round(runif(100, 0,1)),10,10)
# samp  <- mySBM$rSampling(Y)
# ll    <- mySBM$samplingLogLik(Y)
