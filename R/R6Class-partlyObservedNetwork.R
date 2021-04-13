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
    obs      = NULL, # collection of observed dyads
    miss     = NULL, # collection of missing dyads
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
    #' @field missingDyads array indices of missing dyads
    missingDyads    = function(value) {private$miss},
    #' @field observedDyads array indices of observed dyads
    observedDyads   = function(value) {private$obs},
    #' @field samplingMatrix matrix of observed and non-observed edges
    samplingMatrix  = function(value) {private$R},
    #' @field samplingMatrixBar matrix of observed and non-observed edges
    samplingMatrixBar  = function(value) {private$S},
    #' @field observedNodes a vector of observed and non-observed nodes (observed means at least one non NA value)
    observedNodes   = function(value) {
      res <- rep(TRUE, self$nbNodes)
      res[Matrix::rowSums(private$S) > 0] <- FALSE
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
      stopifnot(all.equal(sort(unique(as.numeric(adjacencyMatrix[!is.na(adjacencyMatrix)]))), c(0,1)))

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
      private$obs   <- which(!is.na(adjacencyMatrix) & dyads, arr.ind = TRUE )
      private$miss  <- which( is.na(adjacencyMatrix) & dyads, arr.ind = TRUE )
      ## where are my non-zero entries?
      nzero <- which(!is.na(adjacencyMatrix) & adjacencyMatrix != 0 & dyads, arr.ind = TRUE)

      ## sampling matrix (indicating who is observed)
      private$R <- Matrix::sparseMatrix(private$obs[,1] , private$obs[,2] ,x = 1, dims = dim(adjacencyMatrix))
      private$S <- Matrix::sparseMatrix(private$miss[,1], private$miss[,2],x = 1, dims = dim(adjacencyMatrix))
      ## network matrix (only none zero, non NA values)
      private$Y   <- Matrix::sparseMatrix(nzero[,1], nzero[,2], x = 1, dims = dim(adjacencyMatrix))

    },
    #' @description method to cluster network data with missing value
    #' @param vBlocks The vector of number of blocks considered in the collection.
    #' @param method character with a clustering method among "hierarchical", "spectral", "kmeans".
    #' @importFrom stats binomial glm.fit residuals
    #' @importFrom ClusterR KMeans_rcpp
    clustering = function(vBlocks, imputation = c("median", "average") ) {
      A <- self$imputation(imputation)
      ## normalized  Laplacian with Gaussian kernel
      D <- rowSums(A)
      A <- 1/(1 + exp(-A/sd(A)))
      L <- diag(D^(-1/2)) %*% A %*% diag(D^(-1/2))
      U <- eigen(L, symmetric = TRUE)$vectors[,1:max(vBlocks)]
      lapply(vBlocks, function(k)
        as.integer(
          ClusterR::KMeans_rcpp(U[, 1:k, drop = FALSE], k, num_init = 20)$clusters
        )
      )
    },
    #' @description basic imputation from existing clustering
    #' @param type a character, the type of imputation. Either "median" or "average"
    imputation = function(type = c("median", "average")) {
      adjacencyMatrix <- private$Y
      if (!is.null(private$phi)) {
        y <- as.vector(adjacencyMatrix[self$observedDyads])
        X <- cbind(1, apply(private$phi, 3, function(x) x[self$observedDyads]))
        adjacencyMatrix <- matrix(NA, self$nbNodes, self$nbNodes)
### TODO: make it work for other model than Bernoulli / family than binomial
        adjacencyMatrix[self$observedDyads] <- .logistic(residuals(glm.fit(X, y, family = binomial())))
      }
      suppressMessages(adjacencyMatrix[self$missingDyads] <-
        switch(match.arg(type),
               "average"  = mean(adjacencyMatrix, na.rm = TRUE),
               "median"   = median(adjacencyMatrix, na.rm = TRUE)
        ))
      if (!private$directed)
        suppressMessages(adjacencyMatrix[lower.tri(adjacencyMatrix)] <- Matrix::t(adjacencyMatrix)[lower.tri(adjacencyMatrix)])

      adjacencyMatrix
    }
  )
)
