#' @include SBM-class_prop.R

SBM_fit <-
R6Class(classname = "SBM_fit",
  inherit = SBM_new,
  private = list(
    tau = NULL # variational parameters for posterior probablility of class belonging
  ),
  public = list(
    get_argmin = function() {
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
    }),
  active = list(
    blockVarPar = function(value) {
      if (missing(value)) return(private$tau) else  private$tau <- value
    }
  )
)


## overwrite the SBM initialize function
SBM_fit$set("public", "initialize",
  function(sampledNetwork, nBlocks, clusterInit="spectral") {
    nNodes <- nrow(sampledNetwork)
    if (isSymmetric(sampledNetwork)) directed <- FALSE else directed <- TRUE
    if (length(table(sampledNetwork)) > 2) family <- "Poisson" else family <- "Bernoulli"

    if (is.character(clusterInit)) {
      clusterInit <-
        switch(clusterInit,
               "hierarchical" = init_hierarchical(sampledNetwork, nBlocks),
               "kmeans"       = init_kmeans(      sampledNetwork, nBlocks),
                                init_spectral(    sampledNetwork, nBlocks)
        )
    }

    Z <- matrix(0,nNodes,nBlocks)
    Z[cbind(1:nNodes, clusterInit)] <- 1
    pi0 <- (t(Z) %*% sampledNetwork %*% Z) / (t(Z) %*% (1 - diag(nNodes)) %*% Z)
    alpha0 <-  colMeans(Z)

    super$initialize(family, directed, nNodes = nNodes, mixtureParam = alpha0, connectParam = pi0)
    private$tau <- Z
    private$Z   <- Z
    private$X   <- sampledNetwork
  }
)

SBM_fit$set("public", "fixPoint",
  function() {

    if (private$family == "Bernoulli") {
      ## Bernoulli undirected
      tau <- private$X %*% private$tau %*% t(log(private$pi)) + bar(private$X) %*% private$tau %*% t(log(1 - private$pi))
      if (private$directed) {
        ## Bernoulli directed
        tau <- tau + t(private$X) %*% private$tau %*% t(log(t(private$pi))) + t(bar(private$X)) %*% private$tau %*% t(log(1 - t(private$pi)))
      }
    }

    if (private$family == "Poisson") {
      ## Poisson undirected
      tau <- private$X %*% private$tau %*% t(log(private$pi)) + (t(private$X) %*% private$tau %*% log(private$pi)) -
        log(factorial(private$X)*t(factorial(private$X))) %*% private$tau %*% matrix(1,private$Q, private$Q) -
        (matrix(1,private$N, private$N) - diag(private$N)) %*% private$tau %*% t((private$pi + t(private$pi)))
      # if (private$directed) {
      #   ## Poisson directed
      #   tau <- tau + t(private$X) %*% private$tau %*% t(log(private$pi)) + (t(private$X) %*% private$tau %*% log(private$pi)) -
      #            log(factorial(private$X)*t(factorial(private$X))) %*% private$tau %*% matrix(1,private$Q,private$Q) -
      #            (matrix(1,private$N, private$N) - diag(private$N)) %*% private$tau %*% t((private$pi + t(private$pi)))
      # }
    }

    tau <- exp(sweep(tau, 2, log(private$alpha),"+"))
    tau <- tau/rowSums(tau)
    tau[is.nan(tau)] <- .5
    private$tau <- tau
    tau
  }
)

SBM_fit$set("public", "infer",
  function(threshold, maxIter) {

    conv <- vector("numeric", control$maxIter) ; conv[1] <- NA

    i <- 0; cond <- FALSE
    while(!cond){
      i <- i + 1

      pi_old <- self$pi

      self$Mstep()
      self$VEstep()


      if (i > 1) {
        conv[i] <- sqrt(sum((self$pi-pi_old)^2)) / sqrt(sum((pi_old)^2))
        cond    <- (i > maxIter) |  (conv[i] < threshold)
      }
    }

    if (!(class(self$sampling)[1] == "sampling_snowball")){
      self$vICL <- -2 * (self$cLogLik + self$sampling$samplingLogLik(self$sampledNetwork, self$completedNetwork, self$blockVarParam)) +
        self$sampling$penality(self$SBM$nBlocks)
    } else {
      self$vICL <- -2 * self$cLogLik + self$sampling$penality(self$SBM$nBlocks)

    }
  }
)
