context("test network samplers (Class networkSampler and chidren)")

set.seed(178303)
### A SBM model : ###
N <- 500
Q <- 3
alpha <- rep(1, Q)/Q        # mixture parameter
pi <- diag(.45, Q, Q) + .05 # connectivity matrix
gamma <- log(pi/(1 - pi))   # logit transform of pi for the model with covariates
directed <- FALSE

### Draw a SBM model (Bernoulli, undirected)
mySBM <- missSBM::simulate(N, alpha, pi, directed)
A <- mySBM$adjMatrix

### Draw a SBM model (Bernoulli, undirected) with covariates
M <- 10
covariates  <- replicate(M, rnorm(N,mean = 0, sd = 1), simplify = FALSE)
covarMatrix <- simplify2array(covariates)
covarArray  <- missSBM:::getCovarArray(covarMatrix, missSBM:::l1_similarity)
covarParam  <- rnorm(M, 0, 1)
sbm <- missSBM::simulate(N, alpha, gamma, directed, covariates, covarParam)
A_cov <- sbm$adjMatrix

## tolerance for tests
tol <- 1e-2

test_that("Consistency of dyad-centered sampler", {

  ## without covariates
  psi <- 0.1
  mySampler <- missSBM:::simpleDyadSampler$new(psi, N, directed)
  expect_is(mySampler, "simpleDyadSampler")
  expect_equal(mySampler$type, "dyad")
  expect_equal(mySampler$df, 1)
  expect_equal(mySampler$parameters, psi)
  mySampler$rSamplingMatrix()
  expect_equal(dim(mySampler$samplingMatrix), c(N,N))
  psi_hat <- sum(mySampler$samplingMatrix)/N^2
  expect_lt(abs(psi_hat - psi), tol)

  ## with covariates
  psi <- runif(M, -5, 5)
  mySampler <- missSBM:::simpleDyadSampler$new(psi, N, directed, covarArray)
  expect_is(mySampler, "simpleDyadSampler")
  expect_equal(mySampler$type, "dyad")
  expect_equal(mySampler$df, M)
  expect_equal(mySampler$parameters, psi)
  mySampler$rSamplingMatrix()
  expect_equal(dim(mySampler$samplingMatrix), c(N,N))
})

test_that("Consistency of simple node-centered sampling", {

  ## without covariates
  psi <- 0.1
  mySampler <- missSBM:::simpleNodeSampler$new(psi, N, directed)
  expect_is(mySampler, "simpleNodeSampler")
  expect_equal(mySampler$type, "node")
  expect_equal(mySampler$df, 1)
  expect_equal(mySampler$parameters, psi)
  expect_equal(length(mySampler$prob), N)
  mySampler$rSamplingMatrix()
  expect_equal(dim(mySampler$samplingMatrix), c(N,N))
  samplingRate <- sum(mySampler$samplingMatrix)/N^2
  expect_lt(abs(samplingRate - psi*(2 - psi)), .1)

  ## with covariates
  psi <- runif(M, -5, 5)
  mySampler <- missSBM:::simpleNodeSampler$new(psi, N, directed, covarMatrix)
  expect_is(mySampler, "simpleNodeSampler")
  expect_equal(mySampler$type, "node")
  expect_equal(mySampler$df, M)
  expect_equal(mySampler$parameters, psi)
  expect_equal(length(mySampler$prob), N)
  mySampler$rSamplingMatrix()
  expect_equal(dim(mySampler$samplingMatrix), c(N,N))
})

test_that("Consistency of double-standard sampling", {

  psi <- c(0.1, 0.5)
  mySampler <- missSBM:::doubleStandardSampler$new(psi, A, directed)
  expect_is(mySampler, "doubleStandardSampler")
  expect_equal(mySampler$type, "double-standard")
  expect_equal(mySampler$df, 2)
  expect_equal(mySampler$parameters, psi)
  mySampler$rSamplingMatrix()
  expect_equal(dim(mySampler$samplingMatrix), c(N,N))

})

test_that("Consistency of block-node sampling", {

  psi <- c(.1, .2, .7)
  mySampler <- missSBM:::blockNodeSampler$new(psi, N, directed, mySBM$memberships)
  expect_is(mySampler, "blockNodeSampler")
  expect_equal(mySampler$type, "block-node")
  expect_equal(mySampler$df, mySBM$nBlocks)
  expect_equal(mySampler$parameters, psi)
  mySampler$rSamplingMatrix()
  expect_equal(dim(mySampler$samplingMatrix), c(N,N))

})


test_that("Consistency of degree network sampling", {

  psi <- c(.01, .01)
  mySampler <- missSBM:::degreeSampler$new(psi, rowSums(A), directed)
  expect_is(mySampler, "degreeSampler")
  expect_equal(mySampler$type, "degree")
  expect_equal(mySampler$df, 2)
  expect_equal(mySampler$parameters, psi)
  mySampler$rSamplingMatrix()
  expect_equal(dim(mySampler$samplingMatrix), c(N,N))

})

test_that("Consistency of block-dyad sampling", {

  psi <- diag(.45, Q, Q) + .05
  mySampler <- missSBM:::blockDyadSampler$new(psi, N, directed, mySBM$memberships)
  expect_is(mySampler, "blockDyadSampler")
  expect_equal(mySampler$type, "block-dyad")
  expect_equal(mySampler$df, Q * (Q + 1) / 2)
  expect_equal(mySampler$parameters, psi)
  mySampler$rSamplingMatrix()
  expect_equal(dim(mySampler$samplingMatrix), c(N,N))

})

