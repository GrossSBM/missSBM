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
    D        = NULL, # list of potential dyads in the network
    nas      = NULL, # all NA in Y
    D_obs    = NULL, # array indices of missing dyads
    D_miss   = NULL, # array indices of observed dyads
    R        = NULL, # matrix of observed and non-observed edges
    S        = NULL  # vector of observed and non-observed nodes
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## ACTIVE BINDING
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field samplingRate The percentage of observed dyads
    samplingRate = function(value) {length(private$D_obs)/self$nbDyads},
    #' @field nbNodes The number of nodes
    nbNodes = function(value) {ncol(private$Y)},
    #' @field nbDyads The number of dyads
    nbDyads = function(value) {
      ifelse(private$directed, self$nbNodes*(self$nbNodes - 1), self$nbNodes*(self$nbNodes - 1)/2)
    },
    #' @field is_directed logical indicating if the network is directed or not
    is_directed = function(value) {private$directed},
    #' @field netMatrix  The adjacency matrix of the network
    netMatrix = function(value) {private$Y},
    #' @field covarArray the array of covariates
    covarArray = function(value) {private$phi},
    #' @field covarMatrix the matrix of covariates
    covarMatrix = function(value) {if (missing(value)) return(private$X) else  private$X <- value},
    #' @field dyads a list of potential dyads in the network
    dyads           = function(value) {private$D},
    #' @field missingDyads array indices of missing dyads
    missingDyads    = function(value) {private$D_miss},
    #' @field observedDyads array indices of observed dyads
    observedDyads   = function(value) {private$D_obs},
    #' @field samplingMatrix matrix of observed and non-observed edges
    samplingMatrix  = function(value) {private$R},
    #' @field observedNodes a vector of observed and non-observed nodes
    observedNodes   = function(value) {private$S},
    #' @field NAs boolean for NA entries in the adjacencyMatrix
    NAs             = function(value) {private$nas}
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

      ## adjacency matrix
      stopifnot(is.matrix(adjacencyMatrix))

      ## Only binary graph supported
      stopifnot(all.equal(sort(unique(as.numeric(adjacencyMatrix[!is.na(adjacencyMatrix)]))), c(0,1)))

      if (isSymmetric(adjacencyMatrix)) private$directed <- FALSE else private$directed <- TRUE
      private$Y  <- adjacencyMatrix

      ## covariates
      covar <- format_covariates(covariates, similarity)
      private$X   <- covar$Matrix
      private$phi <- covar$Array

      ## sets of observed / unobserved dyads
      private$nas <- is.na(adjacencyMatrix)
      if (private$directed) {
        ## remove diagonal (no loops)
        private$D <- which(upper.tri(adjacencyMatrix) | lower.tri(adjacencyMatrix))
      } else {
        private$D <- which(upper.tri(adjacencyMatrix))
      }
      private$D_miss <- intersect(which( private$nas), private$D )
      private$D_obs  <- intersect(which(!private$nas), private$D )

      ## sets of observed / unobserved nodes
      S <- rep(FALSE, self$nbNodes)
      S[!is.na(rowSums(adjacencyMatrix))] <- TRUE
      private$S <- S

      ## sampling matrix (indicating who is observed) : USELESS ??
      R <- matrix(0, self$nbNodes, self$nbNodes)
      R[private$D_obs] <- 1
      if (!private$directed)  R <- t(R) | R
      private$R <- R
    },
    #' @description method to cluster network data with missing value
    #' @param nbBlocks integer, the chosen number of blocks
    #' @param method character with a clustering method among "hierarchical", "spectral", "kmeans".
    #' @importFrom stats binomial glm.fit residuals
    clustering = function(nbBlocks, method = c("hierarchical", "spectral", "kmeans")) {

      if (nbBlocks > 1) {
        adjacencyMatrix <- private$Y
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
