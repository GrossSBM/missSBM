#' R6 Class definition of an SBM-fit
#'
#' This class is designed to adjust a Stochastic Block Model on a fully observed network.
#' The doVEM method performs inference via Variational EM.
#'
#' This class is virtual: inference an be effective only for one of the two child classes (SBM_fit_nocovariate and SBM_fit_covariates)
#'
#' @import R6
#' @include SBM_fit-Class.R
#' @export
SBM_fit <-
R6::R6Class(classname = "SBM_fit",
  inherit = SBM,
  private = list(
    tau  = NULL  # variational parameters for posterior probablility of class belonging
  ),
  public = list(
    vBound = function(adjMatrix) {self$vExpec(adjMatrix) + self$entropy},
    vICL   = function(adjMatrix) {-2 * self$vExpec(adjMatrix) + self$penalty}
  ),
  active = list(
    blocks      = function(value) {if (missing(value)) return(private$tau) else  private$tau <- value},
    memberships = function(value) {apply(private$tau, 1, which.max)},
    penalty     = function(value) {(self$df_connectParams + self$df_covarParams) * log(self$nDyads) + self$df_mixtureParams * log(self$nNodes)},
    entropy     = function(value) {-sum(xlogx(private$tau))}
  )
)

## TODO: have an initialization specific to covariates/nocovariate objects

## overwrite the SBM initialize function
SBM_fit$set("public", "initialize",
  function(adjacencyMatrix, nBlocks, clusterInit = "spectral") {

    ## WHY SUPPRESSING THAT ???
    # try(
    #   !all.equal(unique(as.numeric(adjacencyMatrix[!is.na(adjacencyMatrix)])), c(0,1)),
    #   stop("Only binary graphs are supported.")
    # )

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
        Z <- matrix(0,self$nNodes,self$nBlocks)
        Z[cbind(1:self$nNodes, clusterInit)] <- 1
      } else if (is.numeric(clusterInit)) {
        Z <- matrix(0,self$nNodes,self$nBlocks)
        Z[cbind(1:self$nNodes, clusterInit)] <- 1
      } else {
        stop("unknown type for initial clustering")
      }
    } else {
      Z <- matrix(1, self$nNodes, self$nBlocks)
    }
    private$tau <- Z

    ## Initialize parameters
    self$init_parameters(adjacencyMatrix)

    invisible(self)
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

