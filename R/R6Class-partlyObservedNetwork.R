#' An R6 Class used for internal representation of a partially observed network
#'
#' This class is not exported to the user
#'
#' @importFrom R6 R6Class
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
    Rbar     = NULL  # complementary matrix of R
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
    #' @field missingDyads array indices of missing dyads
    missingDyads    = function(value) {Matrix::which( private$Rbar != 0, arr.ind = TRUE)},
    #' @field observedDyads array indices of observed dyads
    observedDyads   = function(value) {Matrix::which( private$R != 0, arr.ind = TRUE)},
    #' @field samplingMatrix matrix of observed and non-observed edges
    samplingMatrix  = function(value) {private$R},
    #' @field samplingMatrixBar matrix of observed and non-observed edges
    samplingMatrixBar  = function(value) {private$Rbar},
    #' @field observedNodes a vector of observed and non-observed nodes (observed means at least one non NA value)
    observedNodes   = function(value) {
      res <- rep(TRUE, self$nbNodes)
      res[Matrix::rowSums(private$Rbar) > 0] <- FALSE
      res
     },
    #' @field NAs boolean for NA entries in the adjacencyMatrix
    NAs = function(value) {
      if (private$directed)
        private$Rbar
      else
        private$Rbar | Matrix::t(private$Rbar)
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
    initialize = function(adjacencyMatrix, covariates = NULL, similarity = missSBM:::l1_similarity) {

### TODO: handle the case when adjacencyMatrix is a sparseMatrix with NA
      ## adjacency matrix
      stopifnot(is.matrix(adjacencyMatrix))

      ## Only binary graph supported
      stopifnot(all.equal(sort(unique(as.numeric(adjacencyMatrix[!is.na(adjacencyMatrix)]))), c(0,1)))

      private$directed <- ifelse(isSymmetric(adjacencyMatrix), FALSE, TRUE)

      private$Y  <- adjacencyMatrix

      ## covariates
      covar <- format_covariates(covariates, similarity)
      private$X   <- covar$Matrix
      private$phi <- covar$Array

      ## sets of observed / unobserved dyads
      if (private$directed) {
        dyads <- upper.tri(adjacencyMatrix) | lower.tri(adjacencyMatrix)
      } else {
        dyads <- upper.tri(adjacencyMatrix)
      }
      NAs  <- is.na(adjacencyMatrix)
      miss <- which( NAs & dyads, arr.ind = TRUE)
      obs  <- which(!NAs & dyads, arr.ind = TRUE )

      ## sampling matrix (indicating who is observed)
      private$R    <- Matrix::sparseMatrix(obs [,1], obs [,2], x = 1, dims = dim(adjacencyMatrix))
      private$Rbar <- Matrix::sparseMatrix(miss[,1], miss[,2], x = 1, dims = dim(adjacencyMatrix))

      # ## where are my non-zero entries?
      # nzero <- which(!NAs & adjacencyMatrix != 0 & dyads, arr.ind = TRUE)
      # private$Y   <- Matrix::sparseMatrix(nzero[,1], nzero[,2], x = 1, dims = dim(adjacencyMatrix))

    },
    #' @description method to cluster network data with missing value
    #' @param nbBlocks integer, the chosen number of blocks
    #' @param method character with a clustering method among "hierarchical", "spectral", "kmeans".
    #' @importFrom stats binomial glm.fit residuals
    clustering = function(nbBlocks, method = c("spectral", "hierarchical", "kmeans")) {

      if (nbBlocks > 1) {
        adjacencyMatrix <- as.matrix(private$Y)
        if (!is.null(private$phi)) {
          y <- as.vector(adjacencyMatrix)
          X <- cbind(1, apply(private$phi, 3, as.vector))
          NAs <- as.vector(private$nas)
          adjacencyMatrix <- matrix(NA, self$nbNodes, self$nbNodes)
          adjacencyMatrix[!NAs] <- .logistic(residuals(glm.fit(X[!NAs, ], y[!NAs], family = binomial())))
        }
        clustering <-
          switch(match.arg(method),
                 "spectral"     = init_spectral(    adjacencyMatrix, nbBlocks),
                 "kmeans"       = init_kmeans(      adjacencyMatrix, nbBlocks),
                 "hierarchical" = init_hierarchical(adjacencyMatrix, nbBlocks))
      } else {
        clustering <- rep(1L, self$nbNodes)
      }
      clustering
    },
    #' @description basic imputation from existing clustering
    #' @param clustering a vector with size \code{ncol(adjacencyMatrix)}, providing a user-defined clustering with \code{nbBlocks} levels.
    #' @return an adjacency matrix with imputed values
### TODO: include covariates in the imputation!!!
    imputation = function(clustering) {
      adjancency0 <- private$Y
      adjancency0[private$nas] <- 0
      Z <- clustering_indicator(clustering)
      theta0 <- check_boundaries((t(Z) %*% adjancency0 %*% Z) / (t(Z) %*% (1 - diag(self$nbNodes)) %*% Z))
      imputation <- private$Y
      imputation[private$nas] <- (Z %*% theta0 %*% t(Z))[private$nas]
      imputation
    }
  )
)
