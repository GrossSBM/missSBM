#' R6 Class definition of an SBM-fit
#'
#' This class is designed to adjust a Stochastic Block Model on a fully observed network.
#' The doVEM method performs inference via Variational EM
#'
#' @import R6
#' @include SBM-Class.R
#' @export
SBM_fit <-
R6::R6Class(classname = "SBM_fit",
  inherit = SBM,
  private = list(
    tau  = NULL  # variational parameters for posterior probablility of class belonging
  ),
  public = list(
    init_parameters = function(adjMatrix) { ## NA allowed in adjMatrix
      NAs <- is.na(adjMatrix); adjMatrix[NAs] <- 0
      private$pi    <- check_boundaries((t(private$tau) %*% (adjMatrix * !NAs) %*% private$tau) / (t(private$tau) %*% ((1 - diag(self$nNodes)) * !NAs) %*% private$tau))
      private$alpha <- check_boundaries(colMeans(private$tau))
    },
    update_parameters = function(adjMatrix) { # NA not allowed in adjMatrix (should be imputed)
      private$pi    <- check_boundaries((t(private$tau) %*% adjMatrix %*% private$tau) / (t(private$tau) %*% (1 - diag(self$nNodes)) %*% private$tau))
      private$alpha <- check_boundaries(colMeans(private$tau))
    },
    vExpec = function(adjMatrix) {
      prob   <- private$tau %*% private$pi %*% t(private$tau)
      factor <- ifelse(private$directed, 1, .5)
      adjMatrix_zeroDiag     <- adjMatrix ; diag(adjMatrix_zeroDiag) <- 0           ### Changement ici ###
      adjMatrix_zeroDiag_bar <- 1 - adjMatrix ; diag(adjMatrix_zeroDiag_bar) <- 0   ### Changement ici ###
      sum(private$tau %*% log(private$alpha)) +  factor * sum( adjMatrix_zeroDiag * log(prob) + adjMatrix_zeroDiag_bar *  log(1 - prob))
    },
    vBound = function(adjMatrix) {self$vExpec(adjMatrix) + self$entropy},
    vICL   = function(adjMatrix) {-2 * self$vExpec(adjMatrix) + self$penalty}
  ),
  active = list(
    blocks      = function(value) {if (missing(value)) return(private$tau) else  private$tau <- value},
    memberships = function(value) {apply(private$tau, 1, which.max)},
    penalty     = function(value) {self$df_connectParams * log(self$nDyads) + self$df_mixtureParams * log(self$nNodes)},
    entropy     = function(value) {-sum(xlogx(private$tau))}
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
    super$initialize(
      directed     = ifelse(isSymmetric(adjacencyMatrix), FALSE, TRUE),
      nNodes       = nrow(adjacencyMatrix),
      mixtureParam = rep(NA,nBlocks),
      connectParam = matrix(NA,nBlocks,nBlocks)
    )

    ## Initial Clustering
    if (self$nBlocks > 1) {
      if (is.character(clusterInit)) {
        clusterInit <-
          switch(clusterInit,
            "hierarchical" = init_hierarchical(adjacencyMatrix, self$nBlocks),
            "kmeans"       = init_kmeans(      adjacencyMatrix, self$nBlocks),
                             init_spectral(    adjacencyMatrix, self$nBlocks)
          )
      }
      Z <- matrix(0,self$nNodes,self$nBlocks)
      Z[cbind(1:self$nNodes, clusterInit)] <- 1
    } else {
      Z <- matrix(1, self$nNodes, self$nBlocks)
    }
    private$tau <- Z

    ## Initialize parameters
    self$init_parameters(adjacencyMatrix)

    invisible(self)
  }
)

### !!!TODO!!! Write this in C++ and update each individual coordinate-wise
SBM_fit$set("public", "update_blocks",
  function(adjMatrix, fixPointIter, log_lambda = 0) {
    adjMatrix_bar <- bar(adjMatrix)
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
  function(adjMatrix, threshold = 1e-4, maxIter = 10, fixPointIter = 3, trace = FALSE) {

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

      # Assess convergence
      delta[i] <- sqrt(sum((private$pi - pi_old)^2)) / sqrt(sum((pi_old)^2))
      cond     <- (i > maxIter) |  (delta[i] < threshold)
      cond     <- (i > maxIter)
      objective[i] <- self$vBound(adjMatrix)
    }
    if (trace) cat("\n")
    res <- data.frame(delta = delta[1:i], objective = objective[1:i])
    res
  }
)

