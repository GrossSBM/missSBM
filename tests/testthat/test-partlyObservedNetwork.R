set.seed(1234)
library(sbm)
library(aricode)

context("test partlyObservedNetwork class")

N_cov   <- 300
N_nocov <- 300
Q <- 3
M <- 4
source("utils_test.R", local = TRUE)

psi <- 0.75 # missingness

test_that("partlyObservedNetwork: 'Bernoulli' model, undirected, no covariate, no missing", {

  sampler_undirected_nocov$rNetwork(store = TRUE)

  adjMatrix <- sampler_undirected_nocov$networkData
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix)

  ## check object net
  expect_true(inherits(net, "partlyObservedNetwork"))
  expect_true(inherits(net$networkData   , "dgCMatrix"))
  expect_true(inherits(net$samplingMatrix, "dgCMatrix"))
  expect_true(inherits(net$imputation()  , "dgCMatrix"))
  expect_true(Matrix::isSymmetric(net$imputation()))
  expect_equal(net$nbNodes, N_nocov)
  expect_equal(net$is_directed, FALSE)
  expect_equal(net$covarArray, NULL)
  expect_equal(net$covarMatrix, NULL)
  expect_equal(sum(net$observedNodes), N_nocov)
  expect_equal(net$samplingRate, 1)
  expect_equal(dim(net$networkData), c(N_nocov, N_nocov))
  expect_equal(dim(net$samplingMatrix), c(N_nocov, N_nocov))
  expect_equal(dim(net$samplingMatrix), dim(net$samplingMatrixBar))
  expect_equal(dim(net$imputation()), dim(net$networkData))

  ## check that clustering works
  expect_gt(max(sapply(net$clustering(1:(2*Q)), ARI, sampler_undirected_nocov$memberships)), 0.95)

})

test_that("partlyObservedNetwork: 'Bernoulli' model, directed, no covariate, no missing", {

  sampler_directed_nocov$rNetwork(store = TRUE)

  adjMatrix <- sampler_directed_nocov$networkData
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix)

  ## check object net
  expect_true(inherits(net, "partlyObservedNetwork"))
  expect_true(inherits(net$networkData   , "dgCMatrix"))
  expect_true(inherits(net$samplingMatrix, "dgCMatrix"))
  expect_true(inherits(net$imputation()  , "dgCMatrix"))
  expect_true(!Matrix::isSymmetric(net$imputation()))
  expect_equal(net$nbNodes, N_nocov)
  expect_equal(net$is_directed, TRUE)
  expect_equal(net$covarArray, NULL)
  expect_equal(net$covarMatrix, NULL)
  expect_equal(sum(net$observedNodes), N_nocov)
  expect_equal(net$samplingRate, 1)
  expect_equal(dim(net$networkData), c(N_nocov, N_nocov))
  expect_equal(dim(net$samplingMatrix), c(N_nocov, N_nocov))
  expect_equal(dim(net$samplingMatrix), dim(net$samplingMatrixBar))
  expect_equal(dim(net$imputation()), dim(net$networkData))

  ## check that clustering works
  expect_gt(max(sapply(net$clustering(1:(2*Q)), ARI, sampler_directed_nocov$memberships)), 0.95)

})

test_that("partlyObservedNetwork: 'Bernoulli' model, undirected, no covariate, missing", {

  sampler_undirected_nocov$rNetwork(store = TRUE)

  adjMatrix <- missSBM::observeNetwork(sampler_undirected_nocov$networkData, "dyad", psi)
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix)

  ## check object net
  expect_true(inherits(net, "partlyObservedNetwork"))
  expect_true(inherits(net$networkData   , "dgCMatrix"))
  expect_true(inherits(net$samplingMatrix, "dgCMatrix"))
  expect_true(inherits(net$imputation()  , "dgCMatrix"))
  expect_true(Matrix::isSymmetric(net$imputation()))
  expect_equal(net$nbNodes, N_nocov)
  expect_equal(net$is_directed, FALSE)
  expect_equal(net$covarArray, NULL)
  expect_equal(net$covarMatrix, NULL)
  expect_equal(length(net$observedNodes), N_nocov)
  expect_lt(net$samplingRate, psi + 0.1)
  expect_equal(dim(net$networkData), c(N_nocov, N_nocov))
  expect_equal(dim(net$samplingMatrix), c(N_nocov, N_nocov))
  expect_equal(dim(net$samplingMatrix), dim(net$samplingMatrixBar))
  expect_equal(dim(net$imputation()), dim(net$networkData))

  ## check that clustering works
  expect_gt(max(sapply(net$clustering(1:(2*Q)), ARI, sampler_undirected_nocov$memberships)), 0.95)

  adjMatrix <- missSBM::observeNetwork(sampler_undirected_nocov$networkData, "node", psi)
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix)

  ## check that clustering works
  expect_gt(max(sapply(net$clustering(1:(2*Q)), ARI, sampler_undirected_nocov$memberships)), 0.95)

})

test_that("partlyObservedNetwork: 'Bernoulli' model, directed, no covariate, missing", {

  sampler_directed_nocov$rNetwork(store = TRUE)

  adjMatrix <- missSBM::observeNetwork(sampler_directed_nocov$networkData, "dyad", psi)
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix)

  ## check object net
  expect_true(inherits(net, "partlyObservedNetwork"))
  expect_true(inherits(net$networkData   , "dgCMatrix"))
  expect_true(inherits(net$samplingMatrix, "dgCMatrix"))
  expect_true(inherits(net$imputation()  , "dgCMatrix"))
  expect_true(!Matrix::isSymmetric(net$imputation()))
  expect_equal(net$nbNodes, N_nocov)
  expect_equal(net$is_directed, TRUE)
  expect_equal(net$covarArray, NULL)
  expect_equal(net$covarMatrix, NULL)
  expect_equal(length(net$observedNodes), N_nocov)
  expect_lt(net$samplingRate, psi + 0.1)
  expect_equal(dim(net$networkData), c(N_nocov, N_nocov))
  expect_equal(dim(net$samplingMatrix), c(N_nocov, N_nocov))
  expect_equal(dim(net$samplingMatrix), dim(net$samplingMatrixBar))
  expect_equal(dim(net$imputation()), dim(net$networkData))

  ## check that clustering works
  expect_gt(max(sapply(net$clustering(1:(2*Q)), ARI, sampler_directed_nocov$memberships)), 0.95)

  ## node sampling
  adjMatrix <- missSBM::observeNetwork(sampler_directed_nocov$networkData, "node", psi)
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix)

  ## check that clustering works
  expect_gt(max(sapply(net$clustering(1:(2*Q)), ARI, sampler_directed_nocov$memberships)), 0.95)


})

test_that("partlyObservedNetwork: 'Bernoulli' model, undirected, covariates, no missing", {

  sampler_undirected_cov$rNetwork(store = TRUE)

  adjMatrix <- sampler_undirected_cov$networkData
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix, covariates = covarList_undirected)

  ## check object net
  expect_true(inherits(net, "partlyObservedNetwork"))
  expect_true(inherits(net$networkData   , "dgCMatrix"))
  expect_true(inherits(net$samplingMatrix, "dgCMatrix"))
  expect_true(inherits(net$imputation()  , "dgCMatrix"))
  expect_true(Matrix::isSymmetric(net$imputation()))
  expect_equal(net$nbNodes, N_nocov)
  expect_equal(net$is_directed, FALSE)
  expect_equal(dim(net$covarArray), c(N_cov, N_cov, M))
  expect_equal(sum(net$observedNodes), N_nocov)
  expect_equal(net$samplingRate, 1)
  expect_equal(dim(net$networkData), c(N_nocov, N_nocov))
  expect_equal(dim(net$samplingMatrix), c(N_nocov, N_nocov))
  expect_equal(dim(net$samplingMatrix), dim(net$samplingMatrixBar))
  expect_equal(dim(net$imputation()), dim(net$networkData))

  ## check that clustering works
  expect_gt(max(sapply(net$clustering(1:(2*Q)), ARI, sampler_undirected_cov$memberships)), 0.7)

})

test_that("partlyObservedNetwork: 'Bernoulli' model, directed, covariates, no missing", {

  sampler_directed_cov$rNetwork(store = TRUE)

  adjMatrix <- sampler_directed_cov$networkData
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix, covariates = covarList_directed)

  ## check object net
  expect_true(inherits(net, "partlyObservedNetwork"))
  expect_true(inherits(net$networkData   , "dgCMatrix"))
  expect_true(inherits(net$samplingMatrix, "dgCMatrix"))
  expect_true(inherits(net$imputation()  , "dgCMatrix"))
  expect_true(!Matrix::isSymmetric(net$imputation()))
  expect_equal(net$nbNodes, N_nocov)
  expect_equal(net$is_directed, TRUE)
  expect_equal(dim(net$covarArray), c(N_cov, N_cov, M))
  expect_equal(sum(net$observedNodes), N_nocov)
  expect_equal(net$samplingRate, 1)
  expect_equal(dim(net$networkData), c(N_nocov, N_nocov))
  expect_equal(dim(net$samplingMatrix), c(N_nocov, N_nocov))
  expect_equal(dim(net$samplingMatrix), dim(net$samplingMatrixBar))
  expect_equal(dim(net$imputation()), dim(net$networkData))

  ## check that clustering works
  expect_gt(max(sapply(net$clustering(1:(2*Q)), ARI, sampler_directed_cov$memberships)), 0.5)

})

test_that("partlyObservedNetwork: 'Bernoulli' model, undirected, covariates, missing", {

  sampler_undirected_cov$rNetwork(store = TRUE)

  adjMatrix <- observeNetwork(sampler_undirected_cov$networkData, "dyad", psi)
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix, covariates = covarList_undirected)

  ## check object net
  expect_true(inherits(net, "partlyObservedNetwork"))
  expect_true(inherits(net$networkData   , "dgCMatrix"))
  expect_true(inherits(net$samplingMatrix, "dgCMatrix"))
  expect_true(inherits(net$imputation()  , "dgCMatrix"))
  expect_true(Matrix::isSymmetric(net$imputation()))
  expect_equal(net$nbNodes, N_nocov)
  expect_equal(net$is_directed, FALSE)
  expect_equal(dim(net$covarArray), c(N_cov, N_cov, M))
  expect_equal(length(net$observedNodes), N_nocov)
  expect_lt(net$samplingRate, psi + 0.1)
  expect_equal(dim(net$networkData), c(N_nocov, N_nocov))
  expect_equal(dim(net$samplingMatrix), c(N_nocov, N_nocov))
  expect_equal(dim(net$samplingMatrix), dim(net$samplingMatrixBar))
  expect_equal(dim(net$imputation()), dim(net$networkData))

  ## check that clustering works
  expect_gt(max(sapply(net$clustering(1:(2*Q)), ARI, sampler_undirected_cov$memberships)), 0.5)

  adjMatrix <- observeNetwork(sampler_undirected_cov$networkData, "covar-dyad", runif(M, 0, 2), covariates = covarList_undirected)
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix, covariates = covarList_undirected)

  expect_gt(max(sapply(net$clustering(1:(2*Q)), ARI, sampler_undirected_cov$memberships)), 0.5)

  adjMatrix <- missSBM::observeNetwork(sampler_undirected_cov$networkData, "node", psi)
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix, covariates = covarList_undirected)

  ## check that clustering works
  expect_gt(max(sapply(net$clustering(1:(2*Q)), ARI, sampler_undirected_cov$memberships)), 0.5)

  # adjMatrix <- missSBM::observeNetwork(sampler_undirected_cov$networkData, "covar-node",  runif(M, 0, 2), covariates = covarList_node)
  # net <- missSBM:::partlyObservedNetwork$new(adjMatrix, covariates = covarList_node)
  #
  # ## check that clustering works
  # expect_gt(max(sapply(net$clustering(1:(2*Q)), ARI, sampler_undirected_cov$memberships)), 0.3)

})

test_that("partlyObservedNetwork: 'Bernoulli' model, directed, covariates, missing", {

  sampler_directed_cov$rNetwork(store = TRUE)

  adjMatrix <- observeNetwork(sampler_directed_cov$networkData, "dyad", psi)
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix, covariates = covarList_directed)

  ## check object net
  expect_true(inherits(net, "partlyObservedNetwork"))
  expect_true(inherits(net$networkData   , "dgCMatrix"))
  expect_true(inherits(net$samplingMatrix, "dgCMatrix"))
  expect_true(inherits(net$imputation()  , "dgCMatrix"))
  expect_true(!Matrix::isSymmetric(net$imputation()))
  expect_equal(net$nbNodes, N_nocov)
  expect_equal(net$is_directed, TRUE)
  expect_equal(dim(net$covarArray), c(N_cov, N_cov, M))
  expect_equal(length(net$observedNodes), N_nocov)
  expect_lt(net$samplingRate, psi + 0.1)
  expect_equal(dim(net$networkData), c(N_nocov, N_nocov))
  expect_equal(dim(net$samplingMatrix), c(N_nocov, N_nocov))
  expect_equal(dim(net$samplingMatrix), dim(net$samplingMatrixBar))
  expect_equal(dim(net$imputation()), dim(net$networkData))

  ## check that clustering works
  expect_gt(max(sapply(net$clustering(1:(2*Q)), ARI, sampler_directed_cov$memberships)), 0.5)

  adjMatrix <- observeNetwork(sampler_directed_cov$networkData, "covar-dyad", runif(M, 0, 2), covariates = covarList_directed)
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix, covariates = covarList_directed)

  expect_gt(max(sapply(net$clustering(1:(2*Q)), ARI, sampler_directed_cov$memberships)), 0.5)

  adjMatrix <- missSBM::observeNetwork(sampler_directed_cov$networkData, "node", psi)
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix, covariates = covarList_directed)

  expect_gt(max(sapply(net$clustering(1:(2*Q)), ARI, sampler_directed_cov$memberships)), 0.5)

  adjMatrix <- missSBM::observeNetwork(sampler_directed_cov$networkData, "covar-node", runif(M, 0, 2), covariates = covarList_node)
  net <- missSBM:::partlyObservedNetwork$new(adjMatrix, covariates = covarList_node)

  expect_gt(max(sapply(net$clustering(1:(2*Q)), ARI, sampler_directed_cov$memberships)), 0.5)

})

