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
  dyad <- samplingSBM(A, "dyad", psi, covariates)
  expect_is(dyad, "sampledNetwork", "R6")
  expect_equal(dim(dyad$adjMatrix), dim(A))
  ## expect error if psi and covariates do not have confortable sizes
  expect_error(samplingSBM(A, "dyad", 0.1, covariates))

})

test_that("Consistency of node-centered network sampling", {

  psi <- runif(ncol(covariates), -5, 5)
  node <- samplingSBM(A, "node", psi, covariates)
  expect_is(node, "sampledNetwork", "R6")
  expect_equal(dim(node$adjMatrix), dim(A))
  ## expect error if psi and covariates do not have confortable sizes
  expect_error(samplingSBM(A, "node", 0.1, covariates))

})
