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
    Y        = NULL, # adjacency matrix (sparse encoding: NA and 0 are not encoded)
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
        dyads <- .row(dim(adjacencyMatrix)) != .col(dim(adjacencyMatrix))
      } else {
        dyads <- .row(dim(adjacencyMatrix)) < .col(dim(adjacencyMatrix))
      }
      ## where are my observations?
      NAs   <- is.na(adjacencyMatrix)
      obs   <- which(!NAs & dyads, arr.ind = TRUE )
      miss  <- which( NAs & dyads, arr.ind = TRUE )
      ## where are my non-zero entries?
      nzero <- which(!NAs & adjacencyMatrix != 0 & dyads, arr.ind = TRUE)

      ## sampling matrix (indicating who is observed)
      private$R <- Matrix::sparseMatrix(obs[,1] , obs[,2] ,x = 1, dims = dim(adjacencyMatrix))
      private$S <- Matrix::sparseMatrix(miss[,1], miss[,2],x = 1, dims = dim(adjacencyMatrix))
      ## network matrix (only none zero, non NA values)
      private$Y <- Matrix::sparseMatrix(nzero[,1], nzero[,2], x = 1, dims = dim(adjacencyMatrix))

    },
    #' @description method to cluster network data with missing value
    #' @param vBlocks The vector of number of blocks considered in the collection.
    #' @param imputation character indicating the type of imputation among "median", "average"
    clustering = function(vBlocks,
                            imputation = ifelse(is.null(private$phi), "median", "average")) {
      A   <- self$imputation(imputation)
      res <- spectral_clustering(A, vBlocks)
      res
    },
    #' @description basic imputation from existing clustering
    #' @param type a character, the type of imputation. Either "median" or "average"
    #' @importFrom stats binomial glm.fit residuals
    #' @importFrom Matrix Diagonal
    imputation = function(type = c("median", "average", "zero")) {
      adjMat <- private$Y
      type   <- match.arg(type)
      ## When there are covariartes
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
      # if (type == "median" & median(adjMat, na.rm = TRUE) != 0) {
      #   miss <- which(private$R == 0)
      #   suppressMessages(adjMat[miss]) <- median(adjMat, na.rm = TRUE)
      # }
      # if (type == "average")
      #   suppressMessages(adjMat[miss]) <- mean(adjMat, na.rm = TRUE)

      if (!private$directed)
        suppressMessages(adjMat[lower.tri(adjMat)] <- Matrix::t(adjMat)[lower.tri(adjMat)])

      adjMat
    }
  )
)
