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
R6::R6Class(classname = "SBM",
  public = list(
    ## fields
    nNodes         = NULL, # number of nodes
    nBlocks        = NULL, # number of blocks
    mixtureParam   = NULL, # vector of block parameters (a.k.a. alpha)
    connectParam   = NULL, # vector of model parameters (a.k.a. theta)
    ## methods
    initialize = function(nNodes, mixtureParam, connectParam=NA) {
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
  ),
  private = list(
    zero  = .Machine$double.eps
  )
)

#'
#' @field completeLogLik  a function to compute the completed log-likelihood
#'
#' @export
SBM_BernoulliUndirected <-
R6::R6Class(classname = "SBM_BernoulliUndirected",
  inherit = SBM,
  public = list(
    initialize = function(nNodes=NA, mixtureParam=NA, connectParam=NA) {
      super$initialize(nNodes, mixtureParam, connectParam)
    },
    completeLogLik_MAR = function(blockIndicators, sampledNetwork) {
      network                 <- sampledNetwork$adjacencyMatrix
      network[is.na(network)] <- 0
      # network[is.na(sampledNetwork$adjacencyMatrix)] <- (blockIndicators %*% self$connectParam %*% t(blockIndicators))[is.na(sampledNetwork$adjacencyMatrix)]
      network.bar             <- (1 - network); diag(network.bar) <- 0
      return(sum(blockIndicators %*% log(self$mixtureParam)) +
        .5 * sum( network *(blockIndicators %*% log(self$connectParam) %*% t(blockIndicators)) +
                    network.bar * (blockIndicators %*% log(1-self$connectParam) %*% t(blockIndicators))))
    },
    completeLogLik = function(completedNetwork, blockIndicators) {
      network     <- completedNetwork
      network.bar <- 1 - network ; diag(network.bar) <- 0
      return(sum(blockIndicators %*% log(self$mixtureParam)) +
               .5 * sum( network *(blockIndicators %*% log(self$connectParam) %*% t(blockIndicators)) +
                             network.bar *(blockIndicators %*% log(1-self$connectParam) %*% t(blockIndicators))))
    },
    rSBM = function() {
      blocks <- super$rSBM()$blocks
      adjacencyMatrix <- matrix(rbinom(self$nNodes^2,1, blocks %*% self$connectParam %*% t(blocks)),self$nNodes)
      adjacencyMatrix <- adjacencyMatrix * lower.tri(adjacencyMatrix) + t(adjacencyMatrix * lower.tri(adjacencyMatrix))
      diag(adjacencyMatrix) <- 0
      return(list(blocks = blocks, adjacencyMatrix = as.matrix(adjacencyMatrix)))
    }
  )
)

#' @export
SBM_BernoulliUndirected.fit <-
  R6::R6Class(classname = "SBM_BernoulliUndirected.fit",
          inherit = SBM_BernoulliUndirected,
          public = list(
            initialize = function(nNodes=NA, mixtureParam=NA, connectParam=NA) {
              super$initialize(nNodes, mixtureParam, connectParam)
            },
            maximization = function(SBM, completedNetwork, blockVarParam) {
              SBM$connectParam <- (t(blockVarParam)%*% completedNetwork %*%blockVarParam) / (t(blockVarParam)%*%((1-diag(self$nNodes)))%*%blockVarParam)
              SBM$connectParam[is.nan(SBM$connectParam)] <- private$zero ; SBM$connectParam[SBM$connectParam > 1-private$zero] <- 1-private$zero ; SBM$connectParam[SBM$connectParam < private$zero] <- private$zero
              SBM$mixtureParam <-  colMeans(blockVarParam)
              return(SBM)
            },
            maximization_MAR = function(SBM, completedNetwork, blockVarParam, samplingMatrix) {
              SBM$connectParam <- (t(blockVarParam)%*% (completedNetwork*samplingMatrix) %*%blockVarParam) / (t(blockVarParam)%*%((1-diag(self$nNodes))*samplingMatrix)%*%blockVarParam)
              SBM$connectParam[is.nan(SBM$connectParam)] <- private$zero ; SBM$connectParam[SBM$connectParam > 1-private$zero] <- 1-private$zero ; SBM$connectParam[SBM$connectParam < private$zero] <- private$zero
              SBM$mixtureParam <- colMeans(blockVarParam)
              return(SBM)
            },
            fixPoint = function(SBM, blockVarParam, completedNetwork) {
              completedNetwork.bar <- 1 - completedNetwork; diag(completedNetwork.bar) <- 0
              blockVarParam.new    <- exp(sweep(completedNetwork %*% blockVarParam %*% t(log(SBM$connectParam)) +
                                               completedNetwork.bar %*% blockVarParam %*% t(log(1-SBM$connectParam)),2,log(SBM$mixtureParam),"+"))
              num                  <- rowSums(blockVarParam.new)
              blockVarParam.new    <- blockVarParam.new/num
              blockVarParam.new[is.nan(blockVarParam.new)] <- 0.5
              return(blockVarParam.new)
            },
            fixPoint_MAR = function(SBM, blockVarParam, completedNetwork, samplingMatrix) {
              completedNetwork.bar <- (1 - completedNetwork)*samplingMatrix; diag(completedNetwork.bar) <- 0
              blockVarParam.new    <- exp(sweep(completedNetwork %*% blockVarParam %*% t(log(SBM$connectParam)) +
                                                  completedNetwork.bar %*% blockVarParam %*% t(log(1-SBM$connectParam)),2,log(SBM$mixtureParam),"+"))
              num                  <- rowSums(blockVarParam.new)
              blockVarParam.new    <- blockVarParam.new/num
              blockVarParam.new[is.nan(blockVarParam.new)] <- 0.5
              return(blockVarParam.new)
            },
            updateNu = function(SBM, sampling, sampledNetwork, blockVarParam, completedNetwork, taylorVarParam) {
              PI                   <- log(SBM$connectParam) - log(1-SBM$connectParam)
              completedNetwork.new <- completedNetwork
              networkWithZeros     <- completedNetwork
              networkWithZeros[sampledNetwork$missingDyads] <- 0
              pap <- matrix(rowSums(completedNetwork), nrow = self$nNodes, ncol = self$nNodes, byrow = FALSE) - (completedNetwork-networkWithZeros)
              eph <- switch(class(sampling)[1],
                            "sampling_doubleStandard" = log(1-sampling$missingParam[2]) - log(1-sampling$missingParam[1]) + blockVarParam %*% PI %*% t(blockVarParam),
                            "sampling_class"          = blockVarParam %*% PI %*% t(blockVarParam),
                            "sampling_starDegree"     = blockVarParam %*% PI %*% t(blockVarParam)  - sampling$missingParam[2] + 2*private$g(taylorVarParam)*(sampling$missingParam[1]*sampling$missingParam[2] + (sampling$missingParam[2]^2)*(1+ pap))  + t(2*private$g(taylorVarParam)*(sampling$missingParam[1]*sampling$missingParam[2] + (sampling$missingParam[2]^2)*(1+ pap))))
              eph <- 1/(1+exp(-eph))
              completedNetwork.new[sampledNetwork$missingDyads] <- eph[sampledNetwork$missingDyads]
              return(completedNetwork.new)
            },
            updateKsi = function(sampling, completedNetwork, sampledNetwork){
              networkWithZeros <- completedNetwork
              networkWithZeros[sampledNetwork$missingDyads] <- 0
              Dtilde   <- rowSums(completedNetwork)
              Dchap    <- rowSums((completedNetwork-networkWithZeros)*(1-(completedNetwork-networkWithZeros))) + Dtilde^2
              ksi      <- sqrt(sampling$missingParam[1]^2 + (sampling$missingParam[2]^2)*Dchap + 2*sampling$missingParam[1]*sampling$missingParam[2]*Dtilde)
            }
          ),
          private = list(
            g = function(x){
              return(-(1/(1+exp(-x)) - 0.5)/(0.5*x))
            }
          )
  )


#' @export
SBM_BernoulliDirected <-
  R6::R6Class(classname = "SBM_BernoulliDirected",
          inherit = SBM,
          public = list(
            initialize = function(nNodes=NA, mixtureParam=NA, connectParam=NA) {
              super$initialize(nNodes, mixtureParam, connectParam)
            },
            completeLogLik_MAR = function(blockIndicators, sampledNetwork) {
              network                 <- sampledNetwork$adjacencyMatrix
              network[is.na(sampledNetwork$adjacencyMatrix)] <- (blockIndicators %*% self$connectParam %*% t(blockIndicators))[is.na(sampledNetwork$adjacencyMatrix)]
              network.bar             <- (1 - network); diag(network.bar) <- 0
              return(sum(blockIndicators %*% log(self$mixtureParam)) +
                       sum( network *(blockIndicators %*% log(self$connectParam) %*% t(blockIndicators)) +
                                   network.bar * (blockIndicators %*% log(1-self$connectParam) %*% t(blockIndicators))))
            },
            completeLogLik = function(completedNetwork, blockIndicators) {
              network     <- completedNetwork
              network.bar <- 1 - network ; diag(network.bar) <- 0
              return(sum(blockIndicators %*% log(self$mixtureParam)) +
                       sum( network *(blockIndicators %*% log(self$connectParam) %*% t(blockIndicators)) +
                                   network.bar *(blockIndicators %*% log(1-self$connectParam) %*% t(blockIndicators))))
            },
            rSBM = function() {
              blocks <- super$rSBM()$blocks
              adjacencyMatrix <- matrix(rbinom(self$nNodes^2,1, blocks %*% self$connectParam %*% t(blocks)),self$nNodes)
              diag(adjacencyMatrix) <- 0
              return(list(blocks = blocks, adjacencyMatrix = adjacencyMatrix))
            }
          )
  )

#' @export
SBM_BernoulliDirected.fit <-
  R6::R6Class(classname = "SBM_BernoulliDirected.fit",
          inherit = SBM_BernoulliDirected,
          public = list(
            initialize = function(nNodes=NA, mixtureParam=NA, connectParam=NA) {
              super$initialize(nNodes, mixtureParam, connectParam)
            },
            maximization = function(SBM, completedNetwork, blockVarParam) {
              SBM$connectParam <- (t(blockVarParam)%*% completedNetwork %*%blockVarParam) / (t(blockVarParam)%*%((1-diag(self$nNodes)))%*%blockVarParam)
              SBM$connectParam[is.nan(SBM$connectParam)] <- private$zero ; SBM$connectParam[SBM$connectParam > 1-private$zero] <- 1-private$zero ; SBM$connectParam[SBM$connectParam < private$zero] <- private$zero
              SBM$mixtureParam <-  colMeans(blockVarParam)
              return(SBM)
            },
            maximization_MAR = function(SBM, completedNetwork, blockVarParam, samplingMatrix) {
              SBM$connectParam <- (t(blockVarParam)%*% (completedNetwork*samplingMatrix) %*%blockVarParam) / (t(blockVarParam)%*%((1-diag(self$nNodes))*samplingMatrix)%*%blockVarParam)
              SBM$connectParam[is.nan(SBM$connectParam)] <- private$zero ; SBM$connectParam[SBM$connectParam > 1-private$zero] <- 1-private$zero ; SBM$connectParam[SBM$connectParam < private$zero] <- private$zero
              SBM$mixtureParam <- colMeans(blockVarParam)
              return(SBM)
            },
            fixPoint = function(SBM, blockVarParam, completedNetwork) {
              completedNetwork.bar <- 1 - completedNetwork; diag(completedNetwork.bar) <- 0
              blockVarParam.new    <- exp(sweep(completedNetwork %*% blockVarParam %*% t(log(SBM$connectParam)) +
                                                  completedNetwork.bar %*% blockVarParam %*% t(log(1-SBM$connectParam)) +
                                                  t(completedNetwork) %*% blockVarParam %*% t(log(SBM$connectParam)) +
                                                  t(completedNetwork.bar) %*% blockVarParam %*% t(log(1-SBM$connectParam)),2,log(SBM$mixtureParam),"+"))
              num                  <- rowSums(blockVarParam.new)
              blockVarParam.new    <- blockVarParam.new/num
              blockVarParam.new[is.nan(blockVarParam.new)] <- 0.5
              return(blockVarParam.new)
            },
            fixPoint_MAR = function(SBM, blockVarParam, completedNetwork, samplingMatrix) {
              completedNetwork.bar <- (1 - completedNetwork)*samplingMatrix; diag(completedNetwork.bar) <- 0
              blockVarParam.new    <- exp(sweep(completedNetwork %*% blockVarParam %*% t(log(SBM$connectParam)) +
                                                  completedNetwork.bar %*% blockVarParam %*% t(log(1-SBM$connectParam)) +
                                                  t(completedNetwork) %*% blockVarParam %*% t(log(SBM$connectParam)) +
                                                  t(completedNetwork.bar) %*% blockVarParam %*% t(log(1-SBM$connectParam)),2,log(SBM$mixtureParam),"+"))
              num                  <- rowSums(blockVarParam.new)
              blockVarParam.new    <- blockVarParam.new/num
              blockVarParam.new[is.nan(blockVarParam.new)] <- 0.5
              return(blockVarParam.new)
            },
            updateNu = function(SBM, sampling, sampledNetwork, blockVarParam, completedNetwork, taylorVarParam) {
              PI                   <- log(SBM$connectParam) - log(1-SBM$connectParam)
              completedNetwork.new <- completedNetwork
              networkWithZeros     <- completedNetwork
              networkWithZeros[sampledNetwork$missingDyads] <- 0
              pap <- matrix(rowSums(completedNetwork), nrow = self$nNodes, ncol = self$nNodes, byrow = FALSE) - (completedNetwork-networkWithZeros)
              eph <- switch(class(sampling)[1],
                            "sampling_doubleStandard" = log(1-sampling$missingParam[2]) - log(1-sampling$missingParam[1]) + blockVarParam %*% PI %*% t(blockVarParam),
                            "sampling_class"          = blockVarParam %*% PI %*% t(blockVarParam),
                            "sampling_starDegree"     = blockVarParam %*% PI %*% t(blockVarParam)  - sampling$missingParam[2] +
                              2*private$g(taylorVarParam)*(sampling$missingParam[1]*sampling$missingParam[2] + (sampling$missingParam[2]^2)*(1+ pap))  +
                              t(2*private$g(taylorVarParam)*(sampling$missingParam[1]*sampling$missingParam[2] + (sampling$missingParam[2]^2)*(1+ pap))))
              eph <- 1/(1+exp(-eph))
              completedNetwork.new[sampledNetwork$missingDyads] <- eph[sampledNetwork$missingDyads]
              return(completedNetwork.new)
            },
            updateKsi = function(sampling, completedNetwork, sampledNetwork){
              networkWithZeros <- completedNetwork
              networkWithZeros[sampledNetwork$missingDyads] <- 0
              Dtilde   <- rowSums(completedNetwork)
              Dchap    <- rowSums((completedNetwork-networkWithZeros)*(1-(completedNetwork-networkWithZeros))) + Dtilde^2
              ksi      <- sqrt(sampling$missingParam[1]^2 + (sampling$missingParam[2]^2)*Dchap + 2*sampling$missingParam[1]*sampling$missingParam[2]*Dtilde)
            }
          ),
          private = list(
            g = function(x){
              return(-(1/(1+exp(-x)) - 0.5)/(0.5*x))
            }
          )
  )


#' @export
SBM_PoissonDirected <-
  R6::R6Class(classname = "SBM_PoissonDirected",
          inherit = SBM,
          public = list(
            initialize = function(nNodes=NA, mixtureParam=NA, connectParam=NA) {
              super$initialize(nNodes, mixtureParam, connectParam)
            },
            completeLogLik = function(completedNetwork, blockIndicators) {
                loop.bar    <- matrix(1,self$nNodes,self$nNodes) ; diag(loop.bar) <- 0
                return(sum(blockIndicators%*%log(self$mixtureParam)) +
                  sum(completedNetwork*(blockIndicators%*%log(self$connectParam)%*%t(blockIndicators))) -
                  sum(log(factorial(completedNetwork))*(blockIndicators%*%matrix(1,self$nBlocks,self$nBlocks)%*%t(blockIndicators))) -
                  sum(loop.bar*blockIndicators%*%self$connectParam%*%t(blockIndicators)))
            },
            completeLogLik_MAR = function(blockIndicators, sampledNetwork) {
                network     <- sampledNetwork$adjacencyMatrix
                network[is.na(network)] <- 0;
                loop.bar    <- matrix(1,self$nNodes,self$nNodes) ; diag(loop.bar) <- 0
                return(sum(blockIndicators%*%log(self$mixtureParam)) +
                         sum(network*(blockIndicators%*%log(self$connectParam)%*%t(blockIndicators))) -
                         sum(sampledNetwork$samplingMatrix * log(factorial(network))*(blockIndicators%*%matrix(1,self$nBlocks,self$nBlocks)%*%t(blockIndicators))) -
                         sum((sampledNetwork$samplingMatrix * loop.bar)*blockIndicators%*%self$connectParam%*%t(blockIndicators)))
            },
            rSBM = function() {
              blocks <- super$rSBM()$blocks
              adjacencyMatrix <- matrix(rpois(self$nNodes^2, blocks %*% self$connectParam %*% t(blocks)) ,self$nNodes)
              diag(adjacencyMatrix) <- 0
              return(list(blocks = blocks, adjacencyMatrix = adjacencyMatrix))
            }
          )
  )


#' @export
SBM_PoissonDirected.fit <-
  R6::R6Class(classname = "SBM_PoissonDirected.fit",
          inherit = SBM_PoissonDirected,
          public = list(
            initialize = function(nNodes=NA, mixtureParam=NA, connectParam=NA) {
              super$initialize(nNodes, mixtureParam, connectParam)
            },
            maximization = function(SBM, completedNetwork, blockVarParam) {
              SBM$connectParam <- (t(blockVarParam)%*% completedNetwork %*%blockVarParam) / (t(blockVarParam)%*%((1-diag(self$nNodes)))%*%blockVarParam)
              SBM$connectParam[is.nan(SBM$connectParam)] <- private$zero ; SBM$connectParam[SBM$connectParam < private$zero] <- private$zero
              SBM$mixtureParam <-  colMeans(blockVarParam)
              return(SBM)
            },
            maximization_MAR = function(SBM, completedNetwork, blockVarParam, samplingMatrix) {
              SBM$connectParam <- (t(blockVarParam)%*% (completedNetwork*samplingMatrix) %*%blockVarParam) / (t(blockVarParam)%*%((1-diag(self$nNodes))*samplingMatrix)%*%blockVarParam)
              SBM$connectParam[is.nan(SBM$connectParam)] <- private$zero ; SBM$connectParam[SBM$connectParam < private$zero] <- private$zero
              SBM$mixtureParam <- colMeans(blockVarParam)
              return(SBM)
            },
            fixPoint = function(SBM, blockVarParam, completedNetwork) {
              completedNetwork.bar <- 1 - completedNetwork; diag(completedNetwork.bar) <- 0
              blockVarParam.new    <- exp(sweep(completedNetwork %*% blockVarParam %*% t(log(SBM$connectParam)) + (t(completedNetwork) %*% blockVarParam %*% log(SBM$connectParam)) -
                                                  log(factorial(completedNetwork)*t(factorial(completedNetwork))) %*% blockVarParam %*% matrix(1,self$nBlocks,self$nBlocks) -
                                                  (matrix(1,self$nNodes,self$nNodes) - diag(self$nNodes)) %*% blockVarParam %*% t((SBM$connectParam + t(SBM$connectParam))),2,log(SBM$mixtureParam),"+"))
              num                  <- rowSums(blockVarParam.new)
              blockVarParam.new    <- blockVarParam.new/num
              blockVarParam.new[is.nan(blockVarParam.new)] <- 0.5
              return(blockVarParam.new)
            },
            fixPoint_MAR = function(SBM, blockVarParam, completedNetwork, samplingMatrix) {
              blockVarParam.new    <- exp(sweep(completedNetwork %*% blockVarParam %*% t(log(SBM$connectParam)) + (t(completedNetwork) %*% blockVarParam %*% log(SBM$connectParam)) -
                                                  log(factorial(completedNetwork)*t(factorial(completedNetwork))) %*% blockVarParam %*% matrix(1,self$nBlocks,self$nBlocks) -
                                                  (samplingMatrix*(matrix(1,self$nNodes,self$nNodes) - diag(self$nNodes) )) %*% blockVarParam %*% t((SBM$connectParam + t(SBM$connectParam))),2,log(SBM$mixtureParam),"+"))

              # blockVarParam.new    <- exp(log(matrix(self$connectParam,nrow=self$nNodes,ncol=self$nBlocks,byrow=T)) + completedNetwork %*% blockVarParam %*% t(log(SBM$connectParam)) + (t(completedNetwork) %*% blockVarParam %*% log(SBM$connectParam)) -
              #                                     log(factorial(completedNetwork)*t(factorial(completedNetwork))) %*% blockVarParam %*% matrix(1,self$nBlocks,self$nBlocks) -
              #                                     (samplingMatrix*(matrix(1,self$nNodes,self$nNodes) - diag(self$nNodes) )) %*% blockVarParam %*% t(SBM$connectParam + t(SBM$connectParam)))

              num                  <- rowSums(blockVarParam.new)
              blockVarParam.new    <- blockVarParam.new/num
              blockVarParam.new[is.nan(blockVarParam.new)] <- 0.5
              return(blockVarParam.new)
            },
            updateNu = function(SBM, sampling, sampledNetwork, blockVarParam, completedNetwork, taylorVarParam) {
              PI                   <- log(SBM$connectParam) - log(1-SBM$connectParam)
              completedNetwork.new <- completedNetwork
              # networkWithZeros     <- completedNetwork
              # networkWithZeros[sampledNetwork$missingDyads] <- 0
              # pap <- matrix(rowSums(completedNetwork), nrow = self$nNodes, ncol = self$nNodes, byrow = FALSE) - (completedNetwork-networkWithZeros)
              eph <- switch(class(sampling)[1],
                            "sampling_doubleStandard" = warning("Sampling not meaningfull in this context"),
                            "sampling_class"          = blockVarParam %*% PI %*% t(blockVarParam),
                            "sampling_starDegree"     = warning("Sampling not meaningfull in this context"))
              eph <- 1/(1+exp(-eph))
              completedNetwork.new[sampledNetwork$missingDyads] <- eph[sampledNetwork$missingDyads]
              return(completedNetwork.new)
            },
            updateKsi = function(sampling, completedNetwork, sampledNetwork){
              warning("Sampling not meaningfull in this context")
            }
          )
  )


#' @export
SBM_PoissonUndirected <-
  R6::R6Class(classname = "SBM_PoissonUndirected",
          inherit = SBM,
          public = list(
            initialize = function(nNodes=NA, mixtureParam=NA, connectParam=NA) {
              super$initialize(nNodes, mixtureParam, connectParam)
            },
            completeLogLik = function(completedNetwork, blockIndicators) {
              loop.bar    <- matrix(1,self$nNodes,self$nNodes) ; diag(loop.bar) <- 0
              return(sum(blockIndicators%*%log(self$mixtureParam)) +
                       .5 * (sum(completedNetwork*(blockIndicators%*%log(self$connectParam)%*%t(blockIndicators))) -
                       sum(log(factorial(completedNetwork))*(blockIndicators%*%matrix(1,self$nBlocks,self$nBlocks)%*%t(blockIndicators))) -
                       sum(loop.bar*blockIndicators%*%self$connectParam%*%t(blockIndicators))))
            },
            completeLogLik_MAR = function(blockIndicators, sampledNetwork) {
              network     <- sampledNetwork$adjacencyMatrix
              network[is.na(network)] <- 0;
              loop.bar    <- matrix(1,self$nNodes,self$nNodes) ; diag(loop.bar) <- 0
              return(sum(blockIndicators%*%log(self$mixtureParam)) +
                       .5 * (sum(network*(blockIndicators%*%log(self$connectParam)%*%t(blockIndicators))) -
                       sum(sampledNetwork$samplingMatrix * log(factorial(network))*(blockIndicators%*%matrix(1,self$nBlocks,self$nBlocks)%*%t(blockIndicators))) -
                       sum((sampledNetwork$samplingMatrix * loop.bar)*blockIndicators%*%self$connectParam%*%t(blockIndicators))))
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

#' @export
SBM_PoissonUndirected.fit <-
  R6::R6Class(classname = "SBM_PoissonUndirected.fit",
          inherit = SBM_PoissonUndirected,
          public = list(
            initialize = function(nNodes=NA, mixtureParam=NA, connectParam=NA) {
              super$initialize(nNodes, mixtureParam, connectParam)
            },
            maximization = function(SBM, completedNetwork, blockVarParam) {
              SBM$connectParam <- (t(blockVarParam)%*% completedNetwork %*%blockVarParam) / (t(blockVarParam)%*%((1-diag(self$nNodes)))%*%blockVarParam)
              SBM$connectParam[is.nan(SBM$connectParam)] <- private$zero ; SBM$connectParam[SBM$connectParam < private$zero] <- private$zero
              SBM$mixtureParam <-  colMeans(blockVarParam)
              return(SBM)
            },
            maximization_MAR = function(SBM, completedNetwork, blockVarParam, samplingMatrix) {
              SBM$connectParam <- (t(blockVarParam)%*% (completedNetwork*samplingMatrix) %*%blockVarParam) / (t(blockVarParam)%*%((1-diag(self$nNodes))*samplingMatrix)%*%blockVarParam)
              SBM$connectParam[is.nan(SBM$connectParam)] <- private$zero ; SBM$connectParam[SBM$connectParam < private$zero] <- private$zero
              SBM$mixtureParam <- colMeans(blockVarParam)
              return(SBM)
            },
            fixPoint = function(SBM, blockVarParam, completedNetwork) {
              completedNetwork.bar <- 1 - completedNetwork; diag(completedNetwork.bar) <- 0
              blockVarParam.new    <- exp(sweep(completedNetwork %*% blockVarParam %*% t(log(SBM$connectParam)) + (t(completedNetwork) %*% blockVarParam %*% log(SBM$connectParam)) -
                                                  log(factorial(completedNetwork)*t(factorial(completedNetwork))) %*% blockVarParam %*% matrix(1,self$nBlocks,self$nBlocks) -
                                                  (matrix(1,self$nNodes,self$nNodes) - diag(self$nNodes)) %*% blockVarParam %*% t((SBM$connectParam + t(SBM$connectParam))),2,log(SBM$mixtureParam),"+"))
              num                  <- rowSums(blockVarParam.new)
              blockVarParam.new    <- blockVarParam.new/num
              blockVarParam.new[is.nan(blockVarParam.new)] <- 0.5
              return(blockVarParam.new)
            },
            fixPoint_MAR = function(SBM, blockVarParam, completedNetwork, samplingMatrix) {
              blockVarParam.new    <- exp(sweep(completedNetwork %*% blockVarParam %*% t(log(SBM$connectParam)) -
                                                  log(factorial(completedNetwork)) %*% blockVarParam %*% matrix(1,self$nBlocks,self$nBlocks) -
                                                  (samplingMatrix*(matrix(1,self$nNodes,self$nNodes) - diag(self$nNodes) )) %*% blockVarParam %*% t(SBM$connectParam),2,log(SBM$mixtureParam),"+"))
              num                  <- rowSums(blockVarParam.new)
              blockVarParam.new    <- blockVarParam.new/num
              blockVarParam.new[is.nan(blockVarParam.new)] <- 0.5
              return(blockVarParam.new)
            },
            updateNu = function(SBM, sampling, sampledNetwork, blockVarParam, completedNetwork, taylorVarParam) {
              PI                   <- log(SBM$connectParam) - log(1-SBM$connectParam)
              completedNetwork.new <- completedNetwork
              # networkWithZeros     <- completedNetwork
              # networkWithZeros[sampledNetwork$missingDyads] <- 0
              # pap <- matrix(rowSums(completedNetwork), nrow = self$nNodes, ncol = self$nNodes, byrow = FALSE) - (completedNetwork-networkWithZeros)
              eph <- switch(class(sampling)[1],
                            "sampling_doubleStandard" = warning("Sampling not meaningfull in this context"),
                            "sampling_class"          = blockVarParam %*% PI %*% t(blockVarParam),
                            "sampling_starDegree"     = warning("Sampling not meaningfull in this context"))
              eph <- 1/(1+exp(-eph))
              completedNetwork.new[sampledNetwork$missingDyads] <- eph[sampledNetwork$missingDyads]
              return(completedNetwork.new)
            },
            updateKsi = function(sampling, completedNetwork, sampledNetwork){
              warning("Sampling not meaningfull in this context")
            }
          )
  )


