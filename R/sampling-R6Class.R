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
  R6Class(classname = "sampling_doubleStandard",
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
            samplingLogLik = function(sampledNetwork, completedNetwork, blockVarParam) {
              ll <- sum(log(self$missingParam[2]) * completedNetwork[sampledNetwork$observedDyads] + log(self$missingParam[1]) * (1-completedNetwork[sampledNetwork$observedDyads])) +
                        sum(log(1-self$missingParam[2]) * completedNetwork[sampledNetwork$missingDyads] + log(1-self$missingParam[1]) * (1-completedNetwork[sampledNetwork$missingDyads]))
              if(!self$directed){
                return(ll/2)
              } else {
                return(ll)
              }
            },
            updatePsi = function(completedNetwork, sampledNetwork, blockVarparam, taylorVarParam) {
              num   <- c(sum(1-completedNetwork[sampledNetwork$observedDyads])-self$nNodes, sum(completedNetwork[sampledNetwork$observedDyads]))
              denom <- c(sum(1-completedNetwork)-self$nNodes, sum(completedNetwork))
              return(num/denom)
            },
            penality = function(nBlocks) {
              return((2 + nBlocks*(nBlocks+1)/2)*log(self$nNodes*(self$nNodes-1)/2) + (nBlocks-1)*log(self$nNodes))
            }
          )
  )

# # Undirected :
# # mySBM <- sampling_doubleStandard$new(10, c(1/4, 1/2))
# # Y <- matrix(round(runif(100, 0,1)),10,10)
# # samp <- mySBM$rSampling(Y)
# # ll <- mySBM$samplingLogLik(Y)
#
# # # Directed :
# # mySBM <- sampling_doubleStandard$new(10, c(1/4, 1/2), directed = TRUE)
# # samp <- mySBM$rSampling(matrix(1,10,10))
# # ll <- mySBM$samplingLogLik(matrix(1,10,10))mySBM <- sampling_doubleStandard$new(10, c(1/4, 1/2))
# # Y    <- matrix(round(runif(100, 0,1)),10,10)
# # samp <- mySBM$rSampling(Y)
# # ll   <- mySBM$samplingLogLik(Y)


#' @export
sampling_class <-
  R6Class(classname = "sampling_class",
          inherit = sampling,
          public = list(
            initialize = function(nNodes, missingParam, directed = FALSE) {
              super$initialize(nNodes, missingParam, directed)
            },
            rSampling = function(adjMatrix, blockVarParam) {
              samplingMatrix <- matrix(0, self$nNodes, self$nNodes)
              sampProb <- runif(self$nNodes) < self$missingParam[apply(blockVarParam, 1, which.max)]
              obsNodes <- which(runif(self$nNodes) < sampProb)

              samplingMatrix <- matrix(0,self$nNodes,self$nNodes) ; samplingMatrix[obsNodes,] <- 1
              diag(samplingMatrix) <- 1
              samplingMatrix <- (t(samplingMatrix) | samplingMatrix)*1

              sampAdjMatrix  <- adjMatrix
              sampAdjMatrix[which(samplingMatrix == 0)] <- NA
              return(sampledNetwork$new(sampAdjMatrix, self$directed))
            },
            samplingLogLik = function(sampledNetwork, completedNetwork, blockVarParam) {
              return(sum(t(sampledNetwork$samplingVector) %*% blockVarParam %*% log(self$missingParam) + t(1-sampledNetwork$samplingVector) %*% blockVarParam %*% log(1-self$missingParam)))
            },
            updatePsi = function(completedNetwork, sampledNetwork, blockVarParam, taylorVarParam) {
              return(colSums(blockVarParam*sampledNetwork$samplingVector)/colSums(blockVarParam))
            },
            penality = function(nBlocks) {
              if(directed){
                return((nBlocks^2)*log(self$nNodes*(self$nNodes-1)) + 2*(nBlocks-1)*log(self$nNodes))
              } else {
                return((nBlocks*(nBlocks+1)/2)*log(self$nNodes*(self$nNodes-1)/2) + 2*(nBlocks-1)*log(self$nNodes))
              }
            }
          )
  )


#' @export
sampling_starDegree <-
  R6Class(classname = "sampling_starDegree",
          inherit = sampling,
          public = list(
            initialize = function(nNodes, missingParam, directed = FALSE) {
              super$initialize(nNodes, missingParam, directed)
            },
            rSampling = function(adjMatrix) {
              samplingMatrix      <- matrix(0, self$nNodes, self$nNodes)
              sampProb            <- self$missingParam[1]+self$missingParam[2]*rowSums(adjMatrix)
              samprob             <- 1/(1+exp(-sampProb))
              obsNodes            <- which(runif(self$nNodes) < sampProb)

              samplingMatrix <- matrix(0,self$nNodes,self$nNodes) ; samplingMatrix[obsNodes,] <- 1
              diag(samplingMatrix) <- 1
              samplingMatrix <- (t(samplingMatrix) | samplingMatrix)*1

              sampAdjMatrix  <- adjMatrix ; sampAdjMatrix[which(samplingMatrix == 0)] <- NA
              return(sampledNetwork$new(sampAdjMatrix, self$directed))
            },
            samplingLogLik = function(sampledNetwork, completedNetwork, blockVarParam) {
              sampProb       <- self$missingParam[1]+self$missingParam[2]*rowSums(completedNetwork)
              samprob        <- 1/(1+exp(-sampProb))
              return(log((sampProb^sampledNetwork$samplingVector)%*%((1-sampProb)^(1-sampledNetwork$samplingVector))) )
            },
            updatePsi = function(completedNetwork, sampledNetwork, blockVarParam, taylorVarParam) {
              networkWithZeros     <- completedNetwork
              networkWithZeros[sampledNetwork$missingDyads] <- 0
              Dtilde <- rowSums(completedNetwork)
              Dchap  <- rowSums((completedNetwork-networkWithZeros)*(1-(completedNetwork-networkWithZeros))) + Dtilde^2
              Nmiss  <- length(which(sampledNetwork$samplingVector == 0))
              b1     <- ( (((2*sum(private$g(taylorVarParam)*Dtilde))*(-length(Nmiss) + 0.5*sampledNetwork$nNodes)))/(sum(private$g(taylorVarParam))) - (-sum(Dtilde[Nmiss]) + sum(Dtilde)*0.5) )
              b2     <- ( 2*sum(private$g(taylorVarParam)*Dchap) - (((2*sum(private$g(taylorVarParam)*Dtilde))^2 ))/(sum(private$g(taylorVarParam))))
              b      <- b1/b2
              a      <- -(b*(2*sum(private$g(taylorVarParam)*Dtilde)) + (-length(Nmiss) + 0.5*sampledNetwork$nNodes))/(sum(private$g(taylorVarParam)))
              psi    <- c(a,b)
              return(psi)
            },
            penality = function(nBlocks) {
              if(directed){
                return((nBlocks^2)*log(self$nNodes*(self$nNodes-1)) + 2*(nBlocks-1)*log(self$nNodes))
              } else {
                return((nBlocks*(nBlocks+1)/2)*log(self$nNodes*(self$nNodes-1)/2) + 2*(nBlocks-1)*log(self$nNodes))
              }
            }
          ),
          private = list(
            g = function(x){
              return(-(1/(1+exp(-x)) - 0.5)/(0.5*x))
            }
          )
  )



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
              if(!self$directed){ samplingMatrix <- t(samplingMatrix) | samplingMatrix }
              diag(samplingMatrix) <- 1

              sampAdjMatrix  <- adjMatrix ; sampAdjMatrix[which(samplingMatrix == 0)] <- NA
              return(sampledNetwork$new(sampAdjMatrix, self$directed))
            },
            samplingLogLik = function(sampledNetwork, completedNetwork, blockVarParam) {
              logPsi    <- ifelse (self$missingParam < .Machine$double.eps, 0, log(self$missingParam))
              log1mPsi  <- ifelse (self$missingParam > 1-.Machine$double.eps, 0, log(1-self$missingParam))
              ll        <- (length(sampledNetwork$observedDyads)-sampledNetwork$nNodes)*logPsi + length(sampledNetwork$missingDyads)*log1mPsi
              if(!self$directed){
                return(ll/2)
              } else {
                return(ll)
              }
            },
            updatePsi = function(completedNetwork, sampledNetwork, blockVarParam, taylorVarParam) {
              if(!sampledNetwork$directed){
                return(length(sampledNetwork$observedDyads)/(2*sampledNetwork$nDyads))
              } else {
                return(length(sampledNetwork$observedDyads)/sampledNetwork$nDyads)
              }
            },
            penality = function(nBlocks) {
              if(self$directed){
                return((1 + (nBlocks^2))*log(self$nNodes*(self$nNodes-1)) + (nBlocks-1)*log(self$nNodes))
              } else {
                return((1 + nBlocks*(nBlocks+1)/2)*log(self$nNodes*(self$nNodes-1)/2) + (nBlocks-1)*log(self$nNodes))
              }
            },
            penalityPoisson = function(nBlocks, samplingMatrix) {
              nObsDyads <- length(which(!is.na(samplingMatrix)))-nrow(samplingMatrix)
              if(self$directed){
                return((nBlocks^2)*log(nObsDyads) + (nBlocks-1)*log(self$nNodes))
              } else {
                return((nBlocks*(nBlocks+1)/2)*log((nObsDyads)/2) + (nBlocks-1)*log(nObsDyads))
              }
            }
          )
  )


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
              obsNodes       <- which(runif(self$nNodes) < rep(self$missingParam, self$nNodes))

              samplingMatrix <- matrix(0,self$nNodes,self$nNodes) ; samplingMatrix[obsNodes,] <- 1
              diag(samplingMatrix) <- 1
              samplingMatrix <- (t(samplingMatrix) | samplingMatrix)*1

              sampAdjMatrix  <- adjMatrix ; sampAdjMatrix[which(samplingMatrix == 0)] <- NA
              return(sampledNetwork$new(sampAdjMatrix, self$directed))
            },
            samplingLogLik = function(sampledNetwork, completedNetwork, blockVarParam) {
              logPsi         <- ifelse (self$missingParam < .Machine$double.eps, 0, log(self$missingParam))
              log1mPsi       <- ifelse (self$missingParam > 1-.Machine$double.eps, 0, log(1-self$missingParam))
              return(logPsi*sum(sampledNetwork$samplingVector) + log1mPsi * sum(1-sampledNetwork$samplingVector))
            },
            updatePsi = function(completedNetwork, sampledNetwork, blockVarParam, taylorVarParam) {
              return(mean(colSums(blockVarParam*sampledNetwork$samplingVector)/colSums(blockVarParam)))
            },
            penality = function(nBlocks) {
              if(self$directed){
                return((nBlocks^2)*log(self$nNodes*(self$nNodes-1)) + nBlocks*log(self$nNodes))
              } else {
                return(nBlocks*(nBlocks+1)/2*log(self$nNodes*(self$nNodes-1)/2) + nBlocks*log(self$nNodes))
              }
            },
            penalityPoisson = function(nBlocks, samplingMatrix) {
              nObsDyads <- length(which(is.na(samplingMatrix)))
              if(self$directed){
                return((nBlocks^2)*log(nObsDyads*(nObsDyads-1)) + nBlocks*log(self$nNodes))
              } else {
                return((nBlocks*(nBlocks+1)/2)*log(nObsDyads*(nObsDyads-1)/2) + nBlocks*log(nObsDyads))
              }
            }
          )
  )


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
            samplingLogLik = function(sampledNetwork, completedNetwork) {
              return(log((self$missingParam^sampledNetwork$samplingVector)%*%((1-self$missingParam)^(1-sampledNetwork$samplingVector))))
            },
            penality = function(nBlocks) {
              if(directed){
                return((nBlocks^2)*log(self$nNodes*(self$nNodes-1)) + nBlocks*log(self$nNodes))
              }
              else {
                return(nBlocks*(nBlocks+1)/2*log(self$nNodes*(self$nNodes-1)/2) + nBlocks*log(self$nNodes))
              }
            }
          )
  )

