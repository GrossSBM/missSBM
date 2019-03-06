context("test-clustering_initialization")

library(aricode)

### A SBM model used for all tests
set.seed(178303)
N <- 400
Q <- 5
alpha <- rep(1,Q)/Q       # mixture parameter
pi <- diag(.45,Q) + .05   # connectivity matrix
directed <- FALSE         # if the network is directed or not

### Draw a SBM model
mySBM <- simulateSBM(N, alpha, pi, directed) # simulation of a Bernoulli non-directed SBM
A_full <- mySBM$adjMatrix             # the adjacency matrix

## Draw random missing entries: MAR case (dyads)
psi <- 0.3
sampledNet <- samplingSBM(A_full, "dyad", psi)
A_dyad <- sampledNet$adjMatrix

psi <- 0.3
sampledNet <- samplingSBM(A_full, "node", psi)
A_node <- sampledNet$adjMatrix

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

      relevance <- .6
      expect_gt(ARI(cl, mySBM$memberships), relevance)

    }
  }
})


