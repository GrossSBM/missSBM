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
    d_law    = NULL, # the probability density distribution of the emission law of the edges
    dyads    = NULL, # set of dyads (differ according to directed/undirected graphs)
    tau      = NULL  # variational parameters for posterior probablility of class belonging
  ),
  public = list(
    init_parameters = function(adjMatrix) { ## NA allowed in adjMatrix
      NAs <- is.na(adjMatrix); adjMatrix[NAs] <- 0
      pi <- (t(private$tau) %*% (adjMatrix * !NAs) %*% private$tau) / (t(private$tau) %*% ((1 - diag(private$N)) * !NAs) %*% private$tau)
      private$pi    <- check_boundaries(pi)
      private$alpha <- colMeans(private$tau)
    },
    update_parameters = function(adjMatrix) { # NA not allowed in adjMatrix (should be imputed)
      pi <- (t(private$tau) %*% adjMatrix %*% private$tau) / (t(private$tau) %*% (1 - diag(private$N)) %*% private$tau)
      private$pi    <- check_boundaries(pi)
      private$alpha <- colMeans(private$tau)
    },
    vBound = function(adjMatrix) {
      JZ <- sum(private$tau %*% log(private$alpha))
      JX <- sum( log( private$d_law(adjMatrix, (private$tau %*% private$pi %*% t(private$tau)))  ) )
      JZ + JX + self$entropy
    },
    vBIC = function(adjMatrix) {
      -2 * self$vBound(adjMatrix) + self$penalty
    },
    vICL = function(adjMatrix) {
      -2 * (self$vBound(adjMatrix) - self$entropy) + self$penalty
    }
  ),
  active = list(
    blocks = function(value) {if (missing(value)) return(private$tau) else  private$tau <- value},
    memberships = function(value) {apply(private$tau, 1, which.max)},
    penalty = function(value) {
      self$df_connectParams * log(sum(private$dyads)) + self$df_mixtureParams * log(nrow(private$tau))
    },
    entropy = function(value) { -sum(xlogx(private$tau))}
  )
)

## overwrite the SBM initialize function
SBM_fit$set("public", "initialize",
  function(adjacencyMatrix, nBlocks, clusterInit = "spectral") {

    try(
      !all.equal(unique(as.numeric(adjacencyMatrix[!is.na(adjacencyMatrix)])), c(0,1)),
      stop("Only binary graphs are supported.")
    )

    # Basic fields intialization and call to super constructor
    nNodes <- nrow(adjacencyMatrix)
    if (isSymmetric(adjacencyMatrix)) directed <- FALSE else directed <- TRUE
    super$initialize(directed, nNodes = nNodes,
                     mixtureParam = rep(NA,nBlocks), connectParam = matrix(NA,nBlocks,nBlocks))
    dyads <- matrix(TRUE, nNodes, nNodes); diag(dyads) <- FALSE
    if (!private$directed) dyads[lower.tri(dyads)] <- FALSE
    private$dyads <- dyads
    private$d_law <-  function(x, prob) {prob^x * (1 - prob)^(1 - x)}

    ## Initial Clustering
    if (is.character(clusterInit)) {
      clusterInit <-
        switch(clusterInit,
               "hierarchical" = init_hierarchical(adjacencyMatrix, nBlocks),
               "kmeans"       = init_kmeans(      adjacencyMatrix, nBlocks),
                                init_spectral(    adjacencyMatrix, nBlocks)
        )
    }
    if (nBlocks > 1) {
      Z <- matrix(0,nNodes,nBlocks)
      Z[cbind(1:nNodes, clusterInit)] <- 1
    } else {
      Z <- matrix(1,nNodes,nBlocks)
    }
    private$tau <- Z

    ## Initialize parameters
    self$init_parameters(adjacencyMatrix)

    invisible(self)
  }
)

SBM_fit$set("public", "update_blocks",
  function(adjMatrix, fixPointIter, log_lambda = 0) {
    NAs <- is.na(adjMatrix)
    adjMatrix[NAs] <- 0
    adjMatrix_bar <- bar(adjMatrix) * !NAs
    for (i in 1:fixPointIter) {
      ## Bernoulli undirected
      tau <- adjMatrix %*% private$tau %*% t(log(private$pi)) + adjMatrix_bar %*% private$tau %*% t(log(1 - private$pi)) + log_lambda
      if (private$directed) {
        ## Bernoulli directed
        tau <- tau + t(adjMatrix) %*% private$tau %*% t(log(t(private$pi))) + t(adjMatrix_bar) %*% private$tau %*% t(log(1 - t(private$pi)))
      }

      tau <- exp(sweep(tau, 2, log(private$alpha),"+"))
      tau <- tau/rowSums(tau)
    }
    private$tau <- tau
  }
)

SBM_fit$set("public", "doVEM",
  function(adjMatrix, threshold = 1e-4, maxIter = 100, fixPointIter = 3, trace = FALSE) {

    ## Initialization of quantities that monitor convergence
    delta     <- vector("numeric", maxIter)
    objective <- vector("numeric", maxIter)
    i <- 0; cond <- FALSE

    ## Starting the Variational EM algorithm
    if (trace) cat("\n Adjusting Variational EM for Stochastic Block Model\n")
    while (!cond) {
      i <- i + 1
      if (trace) cat(" iteration #:", i, "\r")

      pi_old <- private$pi # save old value of parameters to assess convergence

      # Variational E-Step
      self$update_blocks(adjMatrix, fixPointIter)
      # M-step
      self$update_parameters(adjMatrix)

      ## Check convergence
      delta[i] <- sqrt(sum((private$pi - pi_old)^2)) / sqrt(sum((pi_old)^2))
      cond     <- (i > maxIter) |  (delta[i] < threshold)
      objective[i] <- self$vBound(adjMatrix)
    }
    if (trace) cat("\n")
    res <- data.frame(delta = delta[1:i], objective = objective[1:i])
    res
  }
)

