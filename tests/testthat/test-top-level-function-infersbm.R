# context("test-test-top-level-function-infersbm")
#
# library(aricode)
#
# set.seed(1890718)
# ### A SBM model : ###
# N <- 400
# Q <- 5
# alpha <- rep(1,Q)/Q       # mixture parameter
# pi <- diag(.45,Q) + .05   # connectivity matrix
# directed <- FALSE         # if the network is directed or not
#
# ### Draw a SBM model
# mySBM <- simulateSBM(N, alpha, pi, directed) # simulation of ad Bernoulli non-directed SBM
# A <- mySBM$adjacencyMatrix             # the adjacency matrix
#
# test_that("infer SBM with dyad sampling works", {
#
#   psi <- 0.3
#   sampledNet <- samplingSBM(adjacencyMatrix, "dyad", psi)
#
#   vBlocks <- 1:8
#   out <- inferSBM(
#     adjacencyMatrix = A,
#     vBlocks         = 1:8,
#     sampling        = "dyad",
#     smoothing       = "forward",
#     control_VEM     = list(threshold = 1e-3, maxIter = 200, fixPointIter = 3)
#   )
#
# })
#
# test_that("infer SBM with node sampling works", {
#
# })
#
# test_that("infer SBM with double standard sampling works", {
#
# })
#
# test_that("infer SBM with block sampling works", {
#
# })
#
# test_that("infer SBM with degree sampling works", {
#
# })
