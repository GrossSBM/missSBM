context("test clustering function used in initialization")

library(aricode)
### A SBM model used for all tests
set.seed(178303)
N <- 50
Q <- 3
pi <- rep(1, Q)/Q           # block proportion
theta <- list(mean = diag(.45, Q, Q) + .05) # connectivity matrix
directed <- FALSE              # if the network is directed or not

### Draw a SBM model
sbm <- sbm::sampleSimpleSBM(N, pi, theta) # simulation of a Bernoulli non-directed SBM

A_full <- sbm$netMatrix             # the adjacency matrix

## Draw random missing entries: MAR case (dyads)
psi <- 0.8
A_dyad <- missSBM::observeNetwork(A_full, "dyad", psi)

psi <- 0.8
A_node <- missSBM::observeNetwork(A_full, "node", psi)

test_that("Spectral clustering is consistent", {

  for (A in list(A_full, A_dyad, A_node)) {

    ## internal function
    cl_spectral_internal <- missSBM:::init_spectral(A, Q)
    expect_is(cl_spectral_internal, "integer")
    expect_equal(length(cl_spectral_internal), N)

    ## top level function
    cl_spectral <- missSBM:::init_spectral(A, Q)
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
    cl_kmeans <- missSBM:::init_kmeans(A, Q)
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
    cl_hierarchical <- missSBM:::init_hierarchical(A,Q)
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
      cl <- switch(method,
                   "spectral"     = missSBM:::init_spectral(A, Q),
                   "kmeans"       = missSBM:::init_kmeans(A, Q),
                   "hierarchical" = missSBM:::init_hierarchical(A, Q)
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
pi <- rep(1,Q)/Q               # block proportion
theta <- list(mean = diag(.45, Q, Q) + .05)    # connectivity matrix
directed <- FALSE

### Draw a SBM model (Bernoulli, undirected) with covariates
M <- 1
covariates_node <- replicate(M, rnorm(N,mean = 0, sd = 1), simplify = FALSE)
covariates_dyad <- replicate(M, matrix(rnorm(N * N ,mean = 0, sd = 1), N, N), simplify = FALSE)
covarParam  <- rnorm(M, -1, 1)

sbm <- sbm::sampleSimpleSBM(N, pi, theta, covariates = covariates_dyad, covariatesParam = covarParam)

test_that("Init clustering with covariate is consistent", {

  A_full <- sbm$netMatrix
  psi <- runif(M, -5, 5)
  A_dyad <- missSBM::observeNetwork(A_full, "covar-dyad", psi, covariates = covariates_dyad)
  A_node <- missSBM::observeNetwork(A_full, "covar-node", psi, covariates = covariates_node)

  X <- cbind(1, apply(sbm$covarArray, 3, as.vector))

  for (A in list(A_full, A_dyad, A_node)) {
    for (method in c("hierarchical", "spectral", "kmeans")) {

    y <- as.vector(A)
    NAs <- as.vector(is.na(A))
    A_ <- matrix(NA, nrow(A), ncol(A))
    A_[!NAs] <- .logistic(residuals(glm.fit(X[!NAs, ], A[!NAs], family = binomial())))

     cl <- switch(method,
                 "spectral"     = missSBM:::init_spectral(A_, Q),
                 "kmeans"       = missSBM:::init_kmeans(A_, Q),
                 "hierarchical" = missSBM:::init_hierarchical(A_, Q)
     )

      relevance <- .4
      expect_is(cl, "integer")
      expect_equal(length(cl), N)
      expect_gt(ARI(cl, sbm$memberships), relevance)
    }
  }
})

