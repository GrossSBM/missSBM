context("test-clustering-initialization-with-covariates")

library(aricode)
library(blockmodels)

set.seed(178303)
### A SBM model : ###
N <- 300
Q <- 3
alpha <- rep(1,Q)/Q                     # mixture parameter
pi <- diag(.45,Q) + .05                 # connectivity matrix
directed <- FALSE

### Draw a SBM model (Bernoulli, undirected) with covariates
M <- 10
covarMatrix <- matrix(rnorm(N*M,mean = 0, sd = 1), N, M)
covarParam  <- rnorm(M,0,1)

mySBM <- simulateSBM(N, alpha, pi, directed, covarMatrix, covarParam)
A_full <- mySBM$adjMatrix

psi <- runif(ncol(covariates), -5, 5)
A_dyad <- samplingSBM(A_full, "dyad", psi, covariates)$adjMatrix
A_node <- samplingSBM(A_full, "node", psi, covariates)$adjMatrix

test_that("Init clustering with covariate is consistent", {

  for (A in list(A_full, A_dyad, A_node)) {
    for (method in c("hierarchical", "spectral", "kmeans")) {
    cl <-
      missSBM:::init_clustering(
        adjacencyMatrix = A,
        nBlocks = Q,
        covariates = mySBM$covarArray,
        clusterInit = method
      )
      relevance <- .6
      expect_is(cl, "integer")
      expect_equal(length(cl), N)
      expect_gt(ARI(cl, mySBM$memberships), relevance)
    }
  }
})

