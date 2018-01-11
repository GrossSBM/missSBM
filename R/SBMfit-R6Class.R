#' @include SBM-R6Class.R

SBM_fit <-
R6Class(classname = "SBM_fit",
  inherit = SBM_new,
  private = list(
    tau = NULL # variational parameters for posterior probablility of class belonging
  ),
  public = list(
    maximization = function() {
      pi <- (t(private$tau) %*% private$X %*% private$tau) / (t(private$tau) %*% (1 - diag(private$N)) %*% private$tau)
      pi[is.nan(pi)] <- zero
      pi[pi > 1 - zero] <- 1 - zero
      pi[pi <     zero] <-     zero
      self$connectParam <- pi
      self$mixtureParam <-  colMeans(private$tau)
    },
    cLogLik = function() {
      loglikZ <- sum(private$Z %*% log(private$alpha))
      loglikX <- sum( log( private$d_law(private$X[private$edges], (private$Z %*% private$pi %*% t(private$Z))[private$edges] )  ) )
      loglikZ + loglikX
    },
    lowerBound = function() {
      JZ <- sum(private$tau %*% log(private$alpha))
      JX <- sum( log( private$d_law(private$X[private$edges], (private$tau %*% private$pi %*% t(private$tau))[private$edges])  ) )
      JZ + JX
    }
    #
    # fixPoint = function(SBM, blockVarParam, completedNetwork) {
    #   completedNetwork.bar <- 1 - completedNetwork; diag(completedNetwork.bar) <- 0
    #   blockVarParam.new    <- exp(sweep(completedNetwork %*% blockVarParam %*% t(log(SBM$connectParam)) +
    #                                       completedNetwork.bar %*% blockVarParam %*% t(log(1-SBM$connectParam)),2,log(SBM$mixtureParam),"+"))
    #   num                  <- rowSums(blockVarParam.new)
    #   blockVarParam.new    <- blockVarParam.new/num
    #   blockVarParam.new[is.nan(blockVarParam.new)] <- 0.5
    #   return(blockVarParam.new)
    # },
    # fixPoint_MAR = function(SBM, blockVarParam, completedNetwork, samplingMatrix) {
    #   completedNetwork.bar <- (1 - completedNetwork)*samplingMatrix; diag(completedNetwork.bar) <- 0
    #   blockVarParam.new    <- exp(sweep(completedNetwork %*% blockVarParam %*% t(log(SBM$connectParam)) +
    #                                       completedNetwork.bar %*% blockVarParam %*% t(log(1-SBM$connectParam)),2,log(SBM$mixtureParam),"+"))
    #   num                  <- rowSums(blockVarParam.new)
    #   blockVarParam.new    <- blockVarParam.new/num
    #   blockVarParam.new[is.nan(blockVarParam.new)] <- 0.5
    #   return(blockVarParam.new)
    # }
  ),
  active = list(
    blockVarPar = function(value) {
      if (missing(value)) return(private$tau) else  private$tau <- value
    }
  )
)

SBM_BernoulliUndirected.fit <-
R6Class(classname = "SBM_BernoulliUndirected.fit",
  inherit = SBM_BernoulliUndirected,
  public = list(
    initialize = function(nNodes=NA, mixtureParam=NA, connectParam=NA) {
      super$initialize(nNodes, mixtureParam, connectParam)
    },
    maximization = function(completedNetwork, blockVarParam) {
      pi <- (t(blockVarParam) %*% completedNetwork %*% blockVarParam) / (t(blockVarParam) %*% ((1 - diag(private$N))) %*% blockVarParam)
      pi[is.nan(pi)] <- zero
      pi[pi > 1 - zero] <- 1 - zero
      pi[pi <     zero] <-     zero
      self$connectParam <- pi
      self$mixtureParam <-  colMeans(blockVarParam)
    },
    maximization_MAR = function(completedNetwork, blockVarParam, samplingMatrix) {
      pi <- (t(blockVarParam) %*% (completedNetwork*samplingMatrix) %*% blockVarParam) / (t(blockVarParam) %*% ((1 - diag(self$nNodes))*samplingMatrix) %*% blockVarParam)
      pi[is.nan(pi)] <- zero
      pi[pi > 1 - zero] <- 1 - zero
      pi[pi <     zero] <-     zero
      self$connectParam <- pi
      self$mixtureParam <-  colMeans(blockVarParam)
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
    maximization = function(completedNetwork, blockVarParam) {
      pi <- (t(blockVarParam) %*% completedNetwork %*% blockVarParam) / (t(blockVarParam) %*% ((1 - diag(private$N))) %*% blockVarParam)
      pi[is.nan(pi)] <- zero
      pi[pi > 1 - zero] <- 1 - zero
      pi[pi <     zero] <-     zero
      self$connectParam <- pi
      self$mixtureParam <-  colMeans(blockVarParam)
    },
    maximization_MAR = function(completedNetwork, blockVarParam, samplingMatrix) {
      pi <- (t(blockVarParam) %*% (completedNetwork*samplingMatrix) %*% blockVarParam) / (t(blockVarParam) %*% ((1 - diag(self$nNodes))*samplingMatrix) %*% blockVarParam)
      pi[is.nan(pi)] <- zero
      pi[pi > 1 - zero] <- 1 - zero
      pi[pi <     zero] <-     zero
      self$connectParam <- pi
      self$mixtureParam <-  colMeans(blockVarParam)
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
    maximization = function(completedNetwork, blockVarParam) {
      pi <- (t(blockVarParam) %*% completedNetwork %*% blockVarParam) / (t(blockVarParam) %*% ((1 - diag(private$N))) %*% blockVarParam)
      pi[is.nan(pi)] <- zero
      pi[pi > 1 - zero] <- 1 - zero
      pi[pi <     zero] <-     zero
      self$connectParam <- pi
      self$mixtureParam <-  colMeans(blockVarParam)
    },
    maximization_MAR = function(completedNetwork, blockVarParam, samplingMatrix) {
      pi <- (t(blockVarParam) %*% (completedNetwork*samplingMatrix) %*% blockVarParam) / (t(blockVarParam) %*% ((1 - diag(self$nNodes))*samplingMatrix) %*% blockVarParam)
      pi[is.nan(pi)] <- zero
      pi[pi > 1 - zero] <- 1 - zero
      pi[pi <     zero] <-     zero
      self$connectParam <- pi
      self$mixtureParam <-  colMeans(blockVarParam)
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
    maximization = function(completedNetwork, blockVarParam) {
      pi <- (t(blockVarParam) %*% completedNetwork %*% blockVarParam) / (t(blockVarParam) %*% ((1 - diag(private$N))) %*% blockVarParam)
      pi[is.nan(pi)] <- zero
      pi[pi > 1 - zero] <- 1 - zero
      pi[pi <     zero] <-     zero
      self$connectParam <- pi
      self$mixtureParam <-  colMeans(blockVarParam)
    },
    maximization_MAR = function(completedNetwork, blockVarParam, samplingMatrix) {
      pi <- (t(blockVarParam) %*% (completedNetwork*samplingMatrix) %*% blockVarParam) / (t(blockVarParam) %*% ((1 - diag(self$nNodes))*samplingMatrix) %*% blockVarParam)
      pi[is.nan(pi)] <- zero
      pi[pi > 1 - zero] <- 1 - zero
      pi[pi <     zero] <-     zero
      self$connectParam <- pi
      self$mixtureParam <-  colMeans(blockVarParam)
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

