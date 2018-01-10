#' @include SBM-R6Class.R

SBM_BernoulliUndirected.fit <-
R6Class(classname = "SBM_BernoulliUndirected.fit",
  inherit = SBM_BernoulliUndirected,
  public = list(
    initialize = function(nNodes=NA, mixtureParam=NA, connectParam=NA) {
      super$initialize(nNodes, mixtureParam, connectParam)
    },
    ## Pourquoi passe-t-on en argument un objet SBM plutôt que de modifier l'objet courant ?
    maximization = function(SBM, completedNetwork, blockVarParam) {
      SBM$connectParam <- (t(blockVarParam) %*% completedNetwork %*% blockVarParam) / (t(blockVarParam) %*% ((1 - diag(private$N))) %*% blockVarParam)
      SBM$connectParam[is.nan(SBM$connectParam)] <- private$zero
      SBM$connectParam[SBM$connectParam > 1 - private$zero] <- 1 - private$zero
      SBM$connectParam[SBM$connectParam <     private$zero] <-     private$zero
      SBM$mixtureParam <-  colMeans(blockVarParam)
      return(SBM)
    },
    ## Même question
    maximization_MAR = function(SBM, completedNetwork, blockVarParam, samplingMatrix) {
      SBM$connectParam <- (t(blockVarParam) %*% (completedNetwork*samplingMatrix) %*% blockVarParam) / (t(blockVarParam) %*% ((1 - diag(self$nNodes))*samplingMatrix) %*% blockVarParam)
      SBM$connectParam[is.nan(SBM$connectParam)] <- private$zero
      SBM$connectParam[SBM$connectParam > 1 - private$zero] <- 1 - private$zero
      SBM$connectParam[SBM$connectParam <     private$zero] <-     private$zero
      SBM$mixtureParam <- colMeans(blockVarParam)
      return(SBM)
    },
    ## ENCORE...
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

SBM_BernoulliDirected.fit <-
R6Class(classname = "SBM_BernoulliDirected.fit",
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
                                          t(completedNetwork) %*% blockVarParam %*% t(log(t(SBM$connectParam))) +
                                          t(completedNetwork.bar) %*% blockVarParam %*% t(log(1-t(SBM$connectParam))),2,log(SBM$mixtureParam),"+"))

      num                  <- rowSums(blockVarParam.new)
      blockVarParam.new    <- blockVarParam.new/num
      blockVarParam.new[is.nan(blockVarParam.new)] <- 0.5
      return(blockVarParam.new)
    },
    fixPoint_MAR = function(SBM, blockVarParam, completedNetwork, samplingMatrix) {
      completedNetwork.bar <- (1 - completedNetwork)*samplingMatrix; diag(completedNetwork.bar) <- 0
      blockVarParam.new    <- exp(sweep(completedNetwork %*% blockVarParam %*% t(log(SBM$connectParam)) +
                                          completedNetwork.bar %*% blockVarParam %*% t(log(1-SBM$connectParam)) +
                                          t(completedNetwork) %*% blockVarParam %*% t(log(t(SBM$connectParam))) +
                                          t(completedNetwork.bar) %*% blockVarParam %*% t(log(1-t(SBM$connectParam))),2,log(SBM$mixtureParam),"+"))

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

SBM_PoissonDirected.fit <-
R6Class(classname = "SBM_PoissonDirected.fit",
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

SBM_PoissonUndirected.fit <-
R6Class(classname = "SBM_PoissonUndirected.fit",
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

