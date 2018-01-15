#' @import R6
#' @export
#' @export
sampling_fit <-
R6Class(classname = "sampling_fit",
  inherit = sampling_model
)

sampling_model_doubleStandard <-
R6Class(classname = "sampling_doubleStandard",
  inherit = sampling_model,
  public = list(
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

sampling_class <-
R6Class(classname = "sampling_class",
  inherit = sampling,
  public = list(
    samplingLogLik = function(sampledNetwork, completedNetwork, blockVarParam) {
      return(sum(t(sampledNetwork$samplingVector) %*% blockVarParam %*% log(self$missingParam) + t(1-sampledNetwork$samplingVector) %*% blockVarParam %*% log(1-self$missingParam)))
    },
    updatePsi = function(completedNetwork, sampledNetwork, blockVarParam, taylorVarParam) {
      return(colSums(blockVarParam*sampledNetwork$samplingVector)/colSums(blockVarParam))
    },
    penality = function(nBlocks) {
      if(self$directed){
        return((nBlocks^2)*log(self$nNodes*(self$nNodes-1)) + 2*(nBlocks-1)*log(self$nNodes))
      } else {
        return((nBlocks*(nBlocks+1)/2)*log(self$nNodes*(self$nNodes-1)/2) + 2*(nBlocks-1)*log(self$nNodes))
      }
    }
  )
)

sampling_starDegree <-
R6Class(classname = "sampling_starDegree",
  inherit = sampling,
  public = list(
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
      if(self$directed){
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

sampling_randomPairMAR <-
R6Class(classname = "sampling_randomPairMAR",
  inherit = sampling,
  public = list(
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

sampling_randomNodesMAR <-
R6Class(classname = "sampling_randomNodesMAR",
  inherit = sampling,
    public = list(
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

sampling_snowball <-
R6Class(classname = "sampling_snowball",
  inherit = sampling,
  public = list(
    samplingLogLik = function(sampledNetwork, completedNetwork, blockVarParam) {
      0
    },
    # samplingLogLik = function(sampledNetwork, completedNetwork) {
    #   return(log((self$missingParam^sampledNetwork$samplingVector)%*%((1-self$missingParam)^(1-sampledNetwork$samplingVector))))
    # },
    penality = function(nBlocks) {
      if(self$directed){
        return((nBlocks^2)*log(self$nNodes*(self$nNodes-1)) + nBlocks*log(self$nNodes))
      }
      else {
        return(nBlocks*(nBlocks+1)/2*log(self$nNodes*(self$nNodes-1)/2) + nBlocks*log(self$nNodes))
      }
    }
  )
)

