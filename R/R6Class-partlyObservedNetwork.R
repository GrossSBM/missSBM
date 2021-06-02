#' An R6 Class used for internal representation of a partially observed network
#'
#' This class is not exported to the user
#'
#' @importFrom R6 R6Class
#' @import Matrix
partlyObservedNetwork <-
  R6::R6Class(classname = "partlyObservedNetwork",
  ## FIELDS : encode network with missing edges
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    Y        = NULL, # adjacency matrix
    X        = NULL, # the covariates matrix
    phi      = NULL, # the covariates array
    directed = NULL, # directed network of not
    R        = NULL, # the sampling matrix (sparse encoding)
    S        = NULL  # the (anti) sampling matrix (sparse encoding)
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## ACTIVE BINDING
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field samplingRate The percentage of observed dyads
    samplingRate = function(value) {sum(private$R)/self$nbDyads},
    #' @field nbNodes The number of nodes
    nbNodes = function(value) {ncol(private$Y)},
    #' @field nbDyads The number of dyads
    nbDyads = function(value) {
      ifelse(private$directed, self$nbNodes*(self$nbNodes - 1), self$nbNodes*(self$nbNodes - 1)/2)
    },
    #' @field is_directed logical indicating if the network is directed or not
    is_directed = function(value) {private$directed},
    #' @field networkData  The adjacency matrix of the network
    networkData = function(value) {private$Y},
    #' @field covarArray the array of covariates
    covarArray = function(value) {private$phi},
    #' @field covarMatrix the matrix of covariates
    covarMatrix = function(value) {if (missing(value)) return(private$X) else  private$X <- value},
    #' @field samplingMatrix matrix of observed and non-observed edges
    samplingMatrix = function(value) {private$R},
    #' @field samplingMatrixBar matrix of observed and non-observed edges
    samplingMatrixBar  = function(value) {private$S},
    #' @field observedNodes a vector of observed and non-observed nodes (observed means at least one non NA value)
    observedNodes   = function(value) {
      if (private$directed) {
        res <- rowSums(private$R) == (self$nbNodes - 1)
      } else {
        res <- (rowSums(private$R | t(private$R))) == (self$nbNodes - 1)
      }
      res
     }
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description constructor
    #' @param adjacencyMatrix The adjacency matrix of the network
    #' @param covariates A list with M entries (the M covariates), each of whom being either a size-N vector or N x N matrix.
    #' @param similarity An R x R -> R function to compute similarities between node covariates. Default is \code{l1_similarity}, that is, -abs(x-y).
    initialize = function(adjacencyMatrix, covariates = list(), similarity = missSBM:::l1_similarity) {

      ## SANITY CHECKS (on data)
      stopifnot(inherits(adjacencyMatrix, "matrix") | inherits(adjacencyMatrix, "dgCMatrix"))
      ## TODO: handle the case when adjacencyMatrix is a sparseMatrix with NA

      ## TODO: later, should also include Poisson/Gaussian models
      stopifnot(all.equal(sort(setdiff(unique(as.numeric(adjacencyMatrix)), NA)), c(0,1)))

      ## type of SBM
      private$directed <- ifelse(isSymmetric(adjacencyMatrix), FALSE, TRUE)

      ## covariates
      ## TODO: for symmetric network, we should only keep the upper triangular part of the covariates
      covar <- format_covariates(covariates, similarity)
      private$X   <- covar$Matrix
      private$phi <- covar$Array

      ## sets of observed / unobserved dyads
      if (private$directed) {
        dyads <- upper.tri(adjacencyMatrix) | lower.tri(adjacencyMatrix)
      } else {
        dyads <- upper.tri(adjacencyMatrix)
      }
      ## where are my observations?
      obs   <- which(!is.na(adjacencyMatrix) & dyads, arr.ind = TRUE )
      miss  <- which( is.na(adjacencyMatrix) & dyads, arr.ind = TRUE )
      ## where are my non-zero entries?
      nzero <- which(!is.na(adjacencyMatrix) & adjacencyMatrix != 0 & dyads, arr.ind = TRUE)

      ## sampling matrix (indicating who is observed)
      private$R <- Matrix::sparseMatrix(obs[,1] , obs[,2] ,x = 1, dims = dim(adjacencyMatrix))
      private$S <- Matrix::sparseMatrix(miss[,1], miss[,2],x = 1, dims = dim(adjacencyMatrix))
      ## network matrix (only none zero, non NA values)
      private$Y   <- Matrix::sparseMatrix(nzero[,1], nzero[,2], x = 1, dims = dim(adjacencyMatrix))

    },
    #' @description method to cluster network data with missing value
    #' @param vBlocks The vector of number of blocks considered in the collection.
    #' @param imputation character indicating the type of imputation among "median", "average"
    #' @importFrom stats binomial glm.fit residuals
    clustering = function(vBlocks,
                          imputation = ifelse(is.null(private$phi), "median", "average")) {

      A <- self$imputation(imputation)
      n <- ncol(A)
      if (self$is_directed) A <- A %*% t(A)
      ## A <- A %*% t(A)
      # A <- as.matrix(1/(1 + exp(-A/sd(A)))) ## caveat: the matrix is dense; pros: lonely node are automatically handled

      ## handling lonely souls
      unconnected <- which(rowSums(abs(A)) == 0)
      connected   <- setdiff(1:n, unconnected)
      A <- A[connected,connected]

      ## normalized Laplacian
      D <- 1/sqrt(rowSums(abs(A)))
      L <- sweep(sweep(A, 1, D, "*"), 2, D, "*")
##      U <- base::svd(L, nu = max(vBlocks), nv = 0)$u
      U <- eigen(L, symmetric = TRUE)$vectors[, 1:max(vBlocks), drop = FALSE]
      res <- future_lapply(vBlocks, function(k) {
        cl <- rep(1L, n)
        if (k != 1) {
          Un <- U[, 1:k, drop = FALSE]
          Un <- sweep(Un, 1, sqrt(rowSums(Un^2)), "/")
          Un[is.nan(Un)] <- 0
          cl_ <- as.integer(
            kmeans_missSBM(Un, k)
          )
         ## handing lonely souls
         cl[connected] <- cl_
         cl[unconnected] <- which.min(rowsum(D, cl_))
        }
        cl
      }, future.seed = TRUE, future.scheduling = structure(TRUE, ordering = "random"))
      res
    },
    #' @description basic imputation from existing clustering
    #' @param type a character, the type of imputation. Either "median" or "average"
    imputation = function(type = c("median", "average", "zero")) {
      adjMat <- private$Y
      if (!is.null(private$phi)) {
        obs <- which(private$R != 0)
        y <- as.vector(adjMat[obs])
        X <- cbind(1, apply(private$phi, 3, function(x) x[obs]))
### TODO: make it work for other model than Bernoulli / family than binomial
        adjMat[obs] <- .logistic(residuals(glm.fit(X, y, family = binomial())))
      }
      miss <- which(private$R == 0)
      suppressMessages(adjMat[miss] <-
        switch(match.arg(type),
               "average"  = mean(adjMat, na.rm = TRUE),
               "median"   = median(adjMat, na.rm = TRUE),
               "zero"     = 0
        ))
      if (!private$directed)
        suppressMessages(adjMat[lower.tri(adjMat)] <- Matrix::t(adjMat)[lower.tri(adjMat)])

      adjMat
    }
  )
)
