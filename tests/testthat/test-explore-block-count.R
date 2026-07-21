context("test that split/merge exploration never corrupts the collection's block count")

test_that("estimateMissSBM() with 'block-node' sampling and many requested blocks keeps vBlocks == requested blocks", {
  skip_on_cran() # network-scale reproduction of a diagnosed bug, too slow for CRAN's routine checks

  frenchblog <- igraph::delete_vertices(missSBM::frenchblog2007, which(igraph::degree(missSBM::frenchblog2007) == 0))
  frenchblog <- igraph::delete_vertices(frenchblog, 61:igraph::vcount(frenchblog))
  blog <- igraph::as_adjacency_matrix(frenchblog, sparse = FALSE)

  set.seed(3052008)
  sbm_full <- estimateMissSBM(blog, 1:6, "node", control = missSBM_param(trace = FALSE, iterates = 1))
  samplingParameters <- ifelse(sbm_full$bestModel$fittedSBM$blockProp < 0.1, 0.2, 0.8)
  blog_obs <- observeNetwork(blog, sampling = "block-node", parameters = samplingParameters,
                              clusters = sbm_full$bestModel$fittedSBM$memberships)

  blocks <- 1:12
  sbm_block <- estimateMissSBM(blog_obs, blocks, "block-node", control = missSBM_param(trace = FALSE, iterates = 2))

  ## before the fix, VEM component collapse during split()/merge() exploration could silently
  ## corrupt the collection's block-count bookkeeping (e.g. vBlocks == c(1,2,3,4,5,4,5,6,6,7,9,10)
  ## instead of 1:12): duplicated/missing slots that make the ICL/ELBO plots nonsensical
  expect_equal(sbm_block$vBlocks, blocks)
  expect_equal(length(unique(sbm_block$vBlocks)), length(sbm_block$vBlocks))
})
