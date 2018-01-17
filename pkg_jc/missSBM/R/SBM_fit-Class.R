#' R6 Class definition of an SBM-fit
#'
#' This class is designed to adjust a Stochastic Block Model on a fully observed network.
#' The doVEM method performs inference via Variational EM
#'
#' @include SBM-Class.R
#' @export
SBM_fit <-
R6Class(classname = "SBM_fit",
  inherit = SBM,
  private = list(
    dyads = NULL, # set of dyads (differ according directed/undirected graph )
    tau   = NULL  # variational parameters for posterior probablility of class belonging
  ),
  public = list(
    update_parameters = function(adjMatrix) {
      pi <- (t(private$tau) %*% adjMatrix %*% private$tau) / (t(private$tau) %*% (1 - diag(private$N)) %*% private$tau)
      pi[is.nan(pi)] <- zero
      pi[pi > 1 - zero] <- 1 - zero
      pi[pi <     zero] <-     zero
      private$pi    <- pi
      private$alpha <- colMeans(private$tau)
    },
    vLogLik = function(adjMatrix) {
      JZ <- sum(private$tau %*% log(private$alpha))
      JX <- sum( log( private$d_law(adjMatrix[private$dyads], (private$tau %*% private$pi %*% t(private$tau))[private$dyads])  ) )
      JZ + JX
    },
    vICL = function(adjMatrix) {
      -2 * self$vLogLik(adjMatrix) + self$penalty
    }),
  active = list(
    blocks = function(value) {
      if (missing(value)) return(private$tau) else  private$tau <- value
    },
    memberships = function(value) {apply(private$tau, 1, which.max)},
    penalty = function(value) {self$df_connectParams * log(sum(private$dyads)) + self$df_mixtureParams * log(private$N)}
  )
)

## overwrite the SBM initialize function
SBM_fit$set("public", "initialize",
  function(adjacencyMatrix, nBlocks, clusterInit="spectral") {
    nNodes <- nrow(adjacencyMatrix)
    if (isSymmetric(adjacencyMatrix)) directed <- FALSE else directed <- TRUE
    if (length(table(adjacencyMatrix)) > 2) family <- "Poisson" else family <- "Bernoulli"

    if (is.character(clusterInit)) {
      clusterInit <-
        switch(clusterInit,
               "hierarchical" = init_hierarchical(adjacencyMatrix, nBlocks),
               "kmeans"       = init_kmeans(      adjacencyMatrix, nBlocks),
                                init_spectral(    adjacencyMatrix, nBlocks)
        )
    }

    Z <- matrix(0,nNodes,nBlocks)
    Z[cbind(1:nNodes, clusterInit)] <- 1
    pi0 <- (t(Z) %*% adjacencyMatrix %*% Z) / (t(Z) %*% (1 - diag(nNodes)) %*% Z)
    alpha0 <-  colMeans(Z)

    super$initialize(family, directed, nNodes = nNodes, mixtureParam = alpha0, connectParam = pi0)

    private$tau <- Z
    dyads <- matrix(TRUE, nNodes, nNodes); diag(dyads) <- FALSE
    if (!private$directed) dyads[lower.tri(dyads)] <- FALSE
    private$dyads <- dyads
  }
)

SBM_fit$set("public", "update_blocks",
  function(adjMatrix, fixPointIter) {
    for (i in 1:fixPointIter) {
      if (private$family == "Bernoulli") {
        ## Bernoulli undirected
        tau <- adjMatrix %*% private$tau %*% t(log(private$pi)) + bar(adjMatrix) %*% private$tau %*% t(log(1 - private$pi))
        if (private$directed) {
          ## Bernoulli directed
          tau <- tau + t(adjMatrix) %*% private$tau %*% t(log(t(private$pi))) + t(bar(adjMatrix)) %*% private$tau %*% t(log(1 - t(private$pi)))
        }
      }

      if (private$family == "Poisson") {
        ## Poisson undirected
        tau <- adjMatrix %*% private$tau %*% t(log(private$pi)) + (t(adjMatrix) %*% private$tau %*% log(private$pi)) -
          log(factorial(adjMatrix)*t(factorial(adjMatrix))) %*% private$tau %*% matrix(1,private$Q, private$Q) -
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
    }
    private$tau <- tau
  }
)

SBM_fit$set("public", "doVEM",
  function(adjMatrix, threshold = 1e-5, maxIter = 100, fixPointIter = 5, trace = FALSE) {

    ## Initialization of quantities that monitor convergence
    delta     <- vector("numeric", maxIter); delta[1] <- NA
    objective <- vector("numeric", maxIter)
    i <- 0; cond <- FALSE

    ## Starting the Variational EM algorithm
    if (trace) cat("\n Adjusting Variational EM for Stochastic Block Model\n")
    while (!cond) {
      i <- i + 1
      if (trace) cat("iteration #:", i, "\r")

      pi_old <- private$pi # save old value of parameters to assess convergence

      self$update_parameters(adjMatrix) # M-step
      self$update_blocks(adjMatrix, fixPointIter)     # Variational E-Step

      ## Check convergence
      if (i > 1) {
        delta[i] <- sqrt(sum((private$pi - pi_old)^2)) / sqrt(sum((pi_old)^2))
        cond     <- (i > maxIter) |  (delta[i] < threshold)
      }
      objective[i] <- self$vLogLik(adjMatrix)
    }
    if (trace) cat("\n")
    res <- data.frame(delta = delta[1:i], objective = objective[1:i])
    res
  }
)

## TODO: add a variationl Bayes procedure ??
