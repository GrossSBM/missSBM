context("test missSBM_collection's warm q-chain initialization (estimate_chain())")

source("utils_test.R", local = TRUE)
sampler_undirected_nocov$rNetwork(store = TRUE)

build_collection <- function(vB) {
  adjMatrix <- missSBM::observeNetwork(sampler_undirected_nocov$networkData, "dyad", .9,
                                        clusters = sampler_undirected_nocov$memberships)
  partlyObservedNet <- missSBM:::partlyObservedNetwork$new(adjMatrix)
  cl <- partlyObservedNet$clustering(vB)
  missSBM_collection$new(
    partlyObservedNet = partlyObservedNet, sampling = "dyad", clusterInit = cl,
    control = list(useCov = FALSE, trace = FALSE))
}

control <- list(threshold = 1e-2, maxIter = 50, fixPointIter = 3, trace = FALSE, iterates = 0)

test_that("estimate_chain() fits a contiguous vBlocks sequence (gap == 1)", {
  collection <- build_collection(1:6)
  out <- collection$estimate_chain(control)
  expect_identical(out, collection) # mutates in place, invisibly returns self

  expect_equal(collection$vBlocks, 1:6)
  expect_true(all(is.finite(collection$ICL)))
  expect_false(any(collection$degenerate))
})

test_that("estimate_chain() handles a non-contiguous vBlocks sequence (gap > 1)", {
  collection <- build_collection(c(1, 3, 6))
  collection$estimate_chain(control)

  expect_equal(collection$vBlocks, c(1, 3, 6))
  expect_true(all(is.finite(collection$ICL)))
})

test_that("estimate_chain() falls back to the cold-started clustering when nothing is splittable", {
  ## candidates_split() requires classes with >= 4 members (see R6Class-missSBM_fit.R); a
  ## 3-node network can never satisfy that, so the chain step from Q=1 to Q=2 must fall back to
  ## this slot's own cold-started placeholder instead of erroring
  adj <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3, 3)
  net <- missSBM:::partlyObservedNetwork$new(adj)
  cl <- list(rep(1, 3), c(1, 1, 2))
  collection <- missSBM_collection$new(
    partlyObservedNet = net, sampling = "dyad", clusterInit = cl,
    control = list(useCov = FALSE, trace = FALSE))

  expect_no_error(collection$estimate_chain(control))
  expect_equal(collection$vBlocks, c(1, 2))
  expect_true(all(is.finite(collection$ICL)))
})

test_that("estimateMissSBM()'s warmChain control routes through estimate_chain()", {
  adj <- missSBM::observeNetwork(sampler_undirected_nocov$networkData, "dyad", .9,
                                  clusters = sampler_undirected_nocov$memberships)
  collection <- estimateMissSBM(adj, vBlocks = 1:5, sampling = "dyad",
                                 control = missSBM_param(trace = FALSE, iterates = 0,
                                                                  polish = FALSE, warmChain = TRUE))
  expect_s3_class(collection, "missSBM_collection")
  expect_equal(collection$vBlocks, 1:5)
  expect_true(all(is.finite(collection$ICL)))
})

test_that("warmChain keeps VEM collapse rare on a network where cold-started exploration degenerates heavily", {
  skip_on_cran() # network-scale reproduction, too slow for CRAN's routine checks

  ## NOTE: a head-to-head warmChain vs. cold comparison at a fixed seed was tried here and
  ## dropped -- future_lapply(future.seed = TRUE) switches the session's RNG kind to
  ## L'Ecuyer-CMRG on first use, and that switch persists for the rest of the R session
  ## (including this test file), so set.seed() before each call does not reproducibly control
  ## the comparison once some earlier test has already triggered the switch. Manually, on a
  ## fresh session, warmChain = TRUE brought collapsed-class models down from 12/14 to 3/14 on
  ## this same scenario (control = missSBM_param(iterates = 0), isolating initialization
  ## quality)
  frenchblog <- igraph::delete_vertices(missSBM::frenchblog2007, which(igraph::degree(missSBM::frenchblog2007) == 0))
  frenchblog <- igraph::delete_vertices(frenchblog, 61:igraph::vcount(frenchblog))
  blog <- igraph::as_adjacency_matrix(frenchblog, sparse = FALSE)

  set.seed(3052008)
  sbm_full <- estimateMissSBM(blog, 1:6, "node", control = missSBM_param(trace = FALSE, iterates = 1, polish = FALSE))
  samplingParameters <- ifelse(sbm_full$bestModel$fittedSBM$blockProp < 0.1, 0.2, 0.8)
  blog_obs <- observeNetwork(blog, sampling = "block-node", parameters = samplingParameters,
                              clusters = sbm_full$bestModel$fittedSBM$memberships)

  blocks <- 1:14
  set.seed(42)
  res_chain <- suppressWarnings(estimateMissSBM(blog_obs, blocks, "block-node",
    control = missSBM_param(trace = FALSE, iterates = 0, polish = FALSE, warmChain = TRUE)))

  expect_true(all(is.finite(res_chain$ICL)))
  expect_lte(sum(res_chain$degenerate), 8) # well below cold-started's ~12/14 on this scenario
})

test_that("estimateMissSBM() sorts and de-duplicates a misordered vBlocks", {
  adj <- missSBM::observeNetwork(sampler_undirected_nocov$networkData, "dyad", .9,
                                  clusters = sampler_undirected_nocov$memberships)
  expect_warning(
    collection <- estimateMissSBM(adj, vBlocks = c(3, 1, 2, 2), sampling = "dyad",
                                   control = missSBM_param(trace = FALSE, iterates = 0)),
    "not strictly increasing"
  )
  expect_equal(collection$vBlocks, 1:3)
})
