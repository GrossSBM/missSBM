context("test-network-data-samplings-with-covariates")

set.seed(178303)
### A SBM model : ###
N <- 300
Q <- 3
alpha <- rep(1,Q)/Q                     # mixture parameter
pi <- diag(.45,Q) + .05                 # connectivity matrix
directed <- FALSE

### Draw a SBM model (Bernoulli, undirected) with covariates
M <- 10
covariates <- matrix(rnorm(N*M,mean = 0, sd = 1), N, M)
covarParam  <- rnorm(M,0,1)

mySBM <- simulateSBM(N, alpha, pi, directed, covariates, covarParam)
A <- mySBM$adjMatrix

test_that("Consistency of dyad-centered network sampling in the presence of covariates", {

  psi <- runif(ncol(covariates), -5, 5)
  sampling_prob_dyad <- missSBM:::logistic(missSBM:::roundProduct(mySBM$covarArray, psi))

  dyad <- samplingSBM(A, "dyad", sampling_prob_dyad)

  expect_is(dyad, "sampledNetwork", "R6")
  expect_equal(dim(dyad$adjMatrix), dim(A))

})

test_that("Consistency of node-centered network sampling", {

  psi <- runif(ncol(covariates), -5, 5)
  sampling_prob_node <- missSBM:::logistic(mySBM$covarMatrix %*% psi)

  node <- samplingSBM(A, "node", sampling_prob_node)
  expect_is(node, "sampledNetwork", "R6")
  expect_equal(dim(node$adjMatrix), dim(A))

})
