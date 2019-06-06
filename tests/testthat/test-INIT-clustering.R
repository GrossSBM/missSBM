context("test clustering function used in initialization")

library(aricode)
### A SBM model used for all tests
set.seed(178303)
N <- 50
Q <- 3
alpha <- rep(1, Q)/Q       # mixture parameter
pi <- diag(.45, Q, Q) + .05   # connectivity matrix
directed <- FALSE         # if the network is directed or not

### Draw a SBM model
sbm <- missSBM::simulate(N, alpha, pi, directed) # simulation of a Bernoulli non-directed SBM

A_full <- sbm$adjacencyMatrix             # the adjacency matrix

## Draw random missing entries: MAR case (dyads)
psi <- 0.8
sampledNet <- missSBM::sample(A_full, "dyad", psi)
A_dyad <- sampledNet$adjacencyMatrix

psi <- 0.8
sampledNet <- missSBM::sample(A_full, "node", psi)
A_node <- sampledNet$adjacencyMatrix

test_that("Spectral clustering is consistent", {

  for (A in list(A_full, A_dyad, A_node)) {

    ## internal function
    cl_spectral_internal <- missSBM:::init_spectral(A, Q)
    expect_is(cl_spectral_internal, "integer")
    expect_equal(length(cl_spectral_internal), N)

    ## top level function
    cl_spectral <-
      missSBM:::init_clustering(
        adjacencyMatrix = A,
        nBlocks = Q,
        clusterInit = "spectral"
      )
    expect_is(cl_spectral, "integer")
    expect_equal(length(cl_spectral), N)

    ## must be equivalent (up to label switching)
    expect_equal(ARI(cl_spectral, cl_spectral_internal), 1.0)
  }

})

test_that("Kmeans clustering is consistent", {

  for (A in list(A_full, A_dyad, A_node)) {

    ## internal function
    cl_kmeans_internal <- missSBM:::init_kmeans(A, Q)
    expect_is(cl_kmeans_internal, "integer")
    expect_equal(length(cl_kmeans_internal), N)

    ## top level function
    cl_kmeans <-
      missSBM:::init_clustering(
        adjacencyMatrix = A,
        nBlocks = Q,
        clusterInit = "kmeans"
      )
    expect_is(cl_kmeans, "integer")
    expect_equal(length(cl_kmeans), N)

    ## must be equivalent (up to label switching)
    expect_equal(ARI(cl_kmeans, cl_kmeans_internal), 1.0)

  }
})

test_that("Hierarchical clustering is consistent", {

  for (A in list(A_full, A_dyad, A_node)) {
    ## internal function
    cl_hierarchical_internal <- missSBM:::init_hierarchical(A, Q)
    expect_is(cl_hierarchical_internal, "integer")
    expect_equal(length(cl_hierarchical_internal), N)

    ## top level function
    cl_hierarchical <-
      missSBM:::init_clustering(
        adjacencyMatrix = A,
        nBlocks = Q,
        clusterInit = "hierarchical"
      )
    expect_is(cl_hierarchical, "integer")
    expect_equal(length(cl_hierarchical), N)

    ## must be equivalent (up to label switching)
    expect_equal(ARI(cl_hierarchical, cl_hierarchical_internal), 1.0)
  }
})


test_that("Clustering initializations are relevant", {

  for (A in list(A_full, A_dyad, A_node)) {

    for (method in c("spectral", "kmeans", "hierarchical")) {

      ## top level function
      cl <-
        missSBM:::init_clustering(
          adjacencyMatrix = A,
          nBlocks = Q,
          clusterInit = method
        )

      relevance <- .5
      expect_gt(ARI(cl, sbm$memberships), relevance)

    }
  }
})


## ========================================================================
## A SBM model with covariates

set.seed(178303)
N <- 40
Q <- 2
alpha <- rep(1,Q)/Q                     # mixture parameter
pi <- diag(.45, Q, Q) + .05                 # connectivity matrix
gamma <- missSBM:::logit(pi)
directed <- FALSE

### Draw a SBM model (Bernoulli, undirected) with covariates
M <- 1
covariates_node <- replicate(M, rnorm(N,mean = 0, sd = 1), simplify = FALSE)
covariates_dyad <- replicate(M, matrix(rnorm(N * N ,mean = 0, sd = 1), N, N), simplify = FALSE)
covarParam  <- rnorm(M, -1, 1)

sbm <- missSBM::simulate(N, alpha, gamma, directed, covariates_dyad, covarParam)

test_that("Init clustering with covariate is consistent", {

  A_full <- sbm$adjacencyMatrix
  psi <- runif(M, -5, 5)
  A_dyad <- missSBM::sample(A_full, "covar-dyad", psi, covariates = covariates_dyad)$adjacencyMatrix
  A_node <- missSBM::sample(A_full, "covar-node", psi, covariates = covariates_node)$adjacencyMatrix

  for (A in list(A_full, A_dyad, A_node)) {
    for (method in c("hierarchical", "spectral", "kmeans")) {
    cl <-
      missSBM:::init_clustering(
        adjacencyMatrix = A,
        nBlocks = Q,
        covarArray = sbm$covarArray,
        clusterInit = method
      )
      relevance <- .4
      expect_is(cl, "integer")
      expect_equal(length(cl), N)
      expect_gt(ARI(cl, sbm$memberships), relevance)
    }
  }
})

