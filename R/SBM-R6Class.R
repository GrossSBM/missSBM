zero <- .Machine$double.eps

bar <- function(X) {
  X.bar <- 1 - X ; diag(X.bar) <- 0
  X.bar
}

#' @import R6
SBM <- # this 'virtual' class is the mother of all subtypes of SBM
R6Class(classname = "SBM",
  ## fields for internal use (refering to mathematical notations)
  private = list(
    N        = NULL, # number of nodes
    Q        = NULL, # number of blocks
    alpha    = NULL, # vector of block parameters (a.k.a. alpha)
    pi       = NULL, # matrix of connectivity (a.k.a. pi)
    Z        = NULL, # a sampled indicator of blocks
    X        = NULL  # a sampled adjacency matrix
  ),
  public = list(
    ## constructor
    initialize = function(nNodes=NA, mixtureParam=NA, connectParam=NA) {
      private$N        <- nNodes
      private$alpha    <- mixtureParam
      private$pi       <- connectParam
      private$Q        <- length(mixtureParam)
    },
    ## a method to generate a vector of clusters indicators
    rBlocks = function() {
      private$Z <- t(rmultinom(private$N, size = 1, prob = private$alpha))
    },
    ## a "virtual" method to sample an adjacency matrix for the current SBM
    rAdjMatrix = function() {NA}
  ),
  active = list(
    ## active binding to access fields
    nNodes       = function(value) {private$N}        , # number of nodes
    nBlocks      = function(value) {private$Q}        , # number of blocks
    isDirected   = function(value) {private$directed} , # directed network or not
    mixtureParam = function(value) {                # vector of block parameters (a.k.a. alpha)
      if (missing(value)) return(private$alpha) else private$alpha <- value
    },
    connectParam = function(value) {                # matrix of connectivity (a.k.a. pi)
      if (missing(value)) return(private$pi) else private$pi <- value
    }   ,
    blocks       = function(value) {                # indicator of blocks
      if (missing(value)) return(private$Z) else private$Z <- value
      }    ,
    clusters     = function(value) {                # vector of clusters
      if (!is.null(private$Z)) apply(private$Z, 1, which.max) else NA
    }    ,
    adjacencyMatrix = function(value) {
      if (missing(value)) return(private$X) else private$X <- value
    }
  )
)

#' @import R6
SBM_new <- # this 'virtual' class is the mother of all subtypes of SBM
R6Class(classname = "SBM_new",
  ## fields for internal use (refering to mathematical notations)
  private = list(
    N        = NULL, # number of nodes
    Q        = NULL, # number of blocks
    alpha    = NULL, # vector of block parameters (a.k.a. alpha)
    pi       = NULL, # matrix of connectivity (a.k.a. pi)
    directed = NULL, # directed or undirected network
    edges    = NULL, # matrix of booleans indicating the observed edges
    family   = NULL, # emission law
    d_law    = NULL, # the density of the emission law of the edges
    r_law    = NULL, # random generation for the emission law of the edges
    Z        = NULL, # a sampled indicator of blocks
    X        = NULL  # a sampled adjacency matrix
  ),
  public = list(
    ## constructor
    initialize = function(family = "Bernoulli", directed = FALSE, nNodes=NA, mixtureParam=NA, connectParam=NA) {
      private$N        <- nNodes
      private$alpha    <- mixtureParam
      private$pi       <- connectParam
      private$Q        <- length(mixtureParam)
      private$family   <- family
      private$directed <- directed
      private$r_law <- switch(family,
                              "Bernoulli" = function(n, prob) {rbinom(n, 1, prob)},
                              "Poisson"   = function(n, prob) {rpois( n,    prob)})
      private$d_law <- switch(family,
                              "Bernoulli" = function(x, prob) {dbinom(x, 1, prob)},
                              "Poisson"   = function(x, prob) {dpois( x,    prob)})
      edges <- matrix(TRUE, nNodes, nNodes); diag(edges) <- FALSE
      if (!directed) edges[lower.tri(edges)] <- FALSE
      private$edges <- edges

    },
    ## a method to generate a vector of clusters indicators
    rBlocks = function() {
      private$Z <- t(rmultinom(private$N, size = 1, prob = private$alpha))
    },
    ## a "virtual" method to sample an adjacency matrix for the current SBM
    rAdjMatrix = function() {
      X <- matrix(private$r_law(private$N^2, private$Z %*% private$pi %*% t(private$Z)), private$N)
      if (!private$directed) X <- X * lower.tri(X) + t(X * lower.tri(X))
      diag(X) <- 0
      private$X <- X
    }
  ),
  active = list(
    ## active binding to access fields
    nNodes       = function(value) {private$N}        , # number of nodes
    nBlocks      = function(value) {private$Q}        , # number of blocks
    isDirected   = function(value) {private$directed} , # directed network or not
    mixtureParam = function(value) {                # vector of block parameters (a.k.a. alpha)
      if (missing(value)) return(private$alpha) else private$alpha <- value
    },
    connectParam = function(value) {                # matrix of connectivity (a.k.a. pi)
      if (missing(value)) return(private$pi) else private$pi <- value
    }   ,
    blocks       = function(value) {                # indicator of blocks
      if (missing(value)) return(private$Z) else private$Z <- value
      }    ,
    clusters     = function(value) {                # vector of clusters
      if (!is.null(private$Z)) apply(private$Z, 1, which.max) else NA
    }    ,
    adjacencyMatrix = function(value) {
      if (missing(value)) return(private$X) else private$X <- value
    }
  )
)

SBM_BernoulliUndirected <-
R6Class(classname = "SBM_BernoulliUndirected",
  inherit = SBM,
  public = list(
    ## effective method to sample an adjacency matrix for the current SBM
    rAdjMatrix = function() {
      X <- matrix(rbinom(private$N^2, 1, private$Z %*% private$pi %*% t(private$Z)), private$N)
      X <- X * lower.tri(X) + t(X * lower.tri(X))
      diag(X) <- 0
      private$X <- X
    },
    completeLogLik_MAR = function(blockIndicators, sampledNetwork) {
      network                 <- sampledNetwork$adjacencyMatrix
      network[is.na(network)] <- 0
      network.bar             <- (1 - network); diag(network.bar) <- 0
      return(sum(blockIndicators %*% log(private$alpha)) +
              .5 * sum( network * (blockIndicators %*% log(private$pi) %*% t(blockIndicators)) +
                network.bar * (blockIndicators %*% log(1 - private$pi) %*% t(blockIndicators))))
    },
    completeLogLik = function(completedNetwork, blockIndicators) {
      network     <- completedNetwork
      network.bar <- 1 - network ; diag(network.bar) <- 0
      return(sum(blockIndicators %*% log(private$alpha)) +
               .5 * sum( network * (blockIndicators %*% log(private$pi) %*% t(blockIndicators)) +
                 network.bar * (blockIndicators %*% log(1 - private$pi) %*% t(blockIndicators))))
    }
  ),
  active = list(
    cLogLik = function(value) {
      X     <- private$X
      X.bar <- 1 - X ; diag(X.bar) <- 0
      return(sum(private$Z %*% log(private$alpha)) +
               .5 * sum( X * (private$Z %*% log(private$pi) %*% t(private$Z)) +
                           X.bar * (private$Z %*% log(1 - private$pi) %*% t(private$Z))))
    }
  )
)

SBM_BernoulliDirected <-
R6Class(classname = "SBM_BernoulliDirected",
  inherit = SBM,
  public = list(
    ## a function to sample an adjacency matrix for the current SBM
    rAdjMatrix = function() {
      X <- matrix(rbinom(private$N^2,1, private$Z %*% private$pi %*% t(private$Z)), private$N)
      diag(X) <- 0
      private$X <- X
    },
    completeLogLik_MAR = function(blockIndicators, sampledNetwork) {
      network                 <- sampledNetwork$adjacencyMatrix
      network[is.na(sampledNetwork$adjacencyMatrix)] <- (blockIndicators %*% self$connectParam %*% t(blockIndicators))[is.na(sampledNetwork$adjacencyMatrix)]
      network.bar             <- (1 - network); diag(network.bar) <- 0
      return(sum(blockIndicators %*% log(self$mixtureParam)) +
               sum( network *(blockIndicators %*% log(self$connectParam) %*% t(blockIndicators)) +
                      network.bar * (blockIndicators %*% log(1-self$connectParam) %*% t(blockIndicators))))
    },
    completeLogLik = function(completedNetwork, blockIndicators) {
      network     <- completedNetwork
      network.bar <- 1 - network ; diag(network.bar) <- 0
      return(sum(blockIndicators %*% log(self$mixtureParam)) +
               sum( network *(blockIndicators %*% log(self$connectParam) %*% t(blockIndicators)) +
                      network.bar *(blockIndicators %*% log(1-self$connectParam) %*% t(blockIndicators))))
    }
  ),
  active = list(
      cLogLik = function(value) {
        X     <- private$X
        X.bar <- 1 - X ; diag(X.bar) <- 0
        return(sum(private$Z %*% log(private$alpha)) +
                 sum( X * (private$Z %*% log(private$pi) %*% t(private$Z)) +
                             X.bar * (private$Z %*% log(1 - private$pi) %*% t(private$Z))))
    }
  )
)

SBM_PoissonDirected <-
R6Class(classname = "SBM_PoissonDirected",
  inherit = SBM,
    public = list(
      ## a function to sample an adjacency matrix for the current SBM
      rAdjMatrix = function() {
        X <- matrix(rpois(private$N^2, private$Z %*% private$pi %*% t(private$Z)), private$N)
        diag(X) <- 0
        private$X <- X
      },
      completeLogLik = function(completedNetwork, blockIndicators) {
        loop.bar    <- matrix(1,self$nNodes,self$nNodes) ; diag(loop.bar) <- 0
        return(sum(blockIndicators%*%log(self$mixtureParam)) +
                 sum(completedNetwork*(blockIndicators%*%log(self$connectParam)%*%t(blockIndicators))) -
                 sum(log(factorial(completedNetwork))*(blockIndicators%*%matrix(1,self$nBlocks,self$nBlocks)%*%t(blockIndicators))) -
                 sum(loop.bar*blockIndicators%*%self$connectParam%*%t(blockIndicators)))
      },
      completeLogLik_MAR = function(blockIndicators, sampledNetwork) {
        network     <- sampledNetwork$adjacencyMatrix
        network[is.na(network)] <- 0;
        loop.bar    <- matrix(1,self$nNodes,self$nNodes) ; diag(loop.bar) <- 0
        return(sum(blockIndicators%*%log(self$mixtureParam)) +
                 sum(network*(blockIndicators%*%log(self$connectParam)%*%t(blockIndicators))) -
                 sum(sampledNetwork$samplingMatrix * log(factorial(network))*(blockIndicators%*%matrix(1,self$nBlocks,self$nBlocks)%*%t(blockIndicators))) -
                 sum((sampledNetwork$samplingMatrix * loop.bar)*blockIndicators%*%self$connectParam%*%t(blockIndicators)))
      }
    ),
  active = list(
    cLogLik = function(value) {
      loop.bar    <- matrix(1,private$N,private$N) ; diag(loop.bar) <- 0
      return(sum(private$Z %*% log(private$alpha)) +
               sum(private$X * (private$Z %*% log(private$pi) %*% t(private$Z))) -
               sum(log(factorial(private$X)) * (private$Z %*% matrix(1,private$Q, private$Q) %*% t(private$Z))) -
               sum(loop.bar * private$Z %*% private$pi %*% t(private$Z)))
    }
  )
)

SBM_PoissonUndirected <-
R6Class(classname = "SBM_PoissonUndirected",
  inherit = SBM,
  public = list(
    ## a function to sample an adjacency matrix for the current SBM
    rAdjMatrix = function() {
      X <- matrix(rpois(private$N^2, private$Z %*% private$pi %*% t(private$Z)), private$N)
      X <- X * lower.tri(X) + t(X * lower.tri(X))
      diag(X) <- 0
      private$X <- X
    },
    completeLogLik = function(completedNetwork, blockIndicators) {
      loop.bar    <- matrix(1,self$nNodes,self$nNodes) ; diag(loop.bar) <- 0
      return(sum(blockIndicators%*%log(self$mixtureParam)) +
               .5 * (sum(completedNetwork*(blockIndicators%*%log(self$connectParam)%*%t(blockIndicators))) -
                       sum(log(factorial(completedNetwork))*(blockIndicators%*%matrix(1,self$nBlocks,self$nBlocks)%*%t(blockIndicators))) -
                       sum(loop.bar*blockIndicators%*%self$connectParam%*%t(blockIndicators))))
    },
    completeLogLik_MAR = function(blockIndicators, sampledNetwork) {
      network     <- sampledNetwork$adjacencyMatrix
      network[is.na(network)] <- 0;
      loop.bar    <- matrix(1,self$nNodes,self$nNodes) ; diag(loop.bar) <- 0
      return(sum(blockIndicators%*%log(self$mixtureParam)) +
               .5 * (sum(network*(blockIndicators%*%log(self$connectParam)%*%t(blockIndicators))) -
                       sum(sampledNetwork$samplingMatrix * log(factorial(network))*(blockIndicators%*%matrix(1,self$nBlocks,self$nBlocks)%*%t(blockIndicators))) -
                       sum((sampledNetwork$samplingMatrix * loop.bar)*blockIndicators%*%self$connectParam%*%t(blockIndicators))))
    }
  ),
  active = list(
      cLogLik = function(value) {
        loop.bar    <- matrix(1,private$N,private$N) ; diag(loop.bar) <- 0
        return(sum(private$Z %*% log(private$alpha)) +
                 .5 * (sum(private$X * (private$Z %*% log(private$pi) %*% t(private$Z))) -
                       sum(log(factorial(private$X)) * (private$Z %*% matrix(1,private$Q, private$Q) %*% t(private$Z))) -
                       sum(loop.bar * private$Z %*% private$pi %*% t(private$Z))))
      }
    )
)

