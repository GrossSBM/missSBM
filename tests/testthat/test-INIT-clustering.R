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

A_full <- sbm$networkData             # the adjacency matrix

## Draw random missing entries: MAR case (dyads)
psi <- 0.7
A_dyad <- missSBM::observeNetwork(A_full, "dyad", psi)

psi <- 0.7
A_node <- missSBM::observeNetwork(A_full, "node", psi)

test_that("Spectral clustering is consistent", {

  for (A in list(A_full, A_dyad, A_node)) {

    myNet <- missSBM:::partlyObservedNetwork$new(A)
    cl <- myNet$clustering(1:Q)

    expect_is(cl, "list")
    expect_is(cl[[Q]], "integer")

    expect_equal(unique(c(sapply(cl, length))), N)
    ## must be equivalent (up to label switching)
    expect_gt(ARI(cl[[Q]], sbm$memberships), 0.8)
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

  A_full <- sbm$networkData
  psi <- runif(M, -5, 5)
  A_dyad <- missSBM::observeNetwork(A_full, "covar-dyad", psi, covariates = covariates_dyad)
  A_node <- missSBM::observeNetwork(A_full, "covar-node", psi, covariates = covariates_node)

  for (A in list(A_full, A_dyad, A_node)) {

    myNet <- missSBM:::partlyObservedNetwork$new(A, covariates = covariates_dyad)
    clusterings <- myNet$clustering(1:Q)
    expect_is(clusterings, "list")

    cl <- clusterings[[Q]]
    relevance <- .6
    expect_is(cl, "integer")
    expect_equal(length(cl), N)
    expect_gt(ARI(cl, sbm$memberships), relevance)

  }
})

