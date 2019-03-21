context("test-misssbm_collection")

set.seed(1890718)
### A SBM model : ###
N <- 300
Q <- 3
alpha <- rep(1, Q)/Q       # mixture parameter
pi <- diag(.45, Q, Q) + .05   # connectivity matrix
directed <- FALSE         # if the network is directed or not

### Draw a SBM model
mySBM <- missSBM::simulate(N, alpha, pi, directed) # simulation of ad Bernoulli non-directed SBM
A <- mySBM$adjMatrix             # the adjacency matrix

mc.cores <- 5


test_that("missSBMcollection works", {

#  l_psi <- list(
#    "dyad" = c(.3),
#    "node" = c(.3),
#    "double-standard" = c(0.4, 0.8),
#    "block-node" = c(.3, .8, .5),
#    "block-dyad" = mySBM$connectParam,
#    "degree" = c(.01, .01)
#  )

#  for (k in seq_along(l_psi)) {

    sampling <- "dyad"

    sampledNet <- missSBM::sample(A, sampling, .5, clusters = mySBM$memberships)

    ## Instatntiate the collection of missingSBM_fit
    collection <- missSBM_collection$new(
      adjMatrix = sampledNet$adjMatrix,
      vBlocks = 1:5,
      sampling = sampling,
      clusterInit = 'hierarchical', NULL, NULL, "none", mc.cores, TRUE)

    ## control parameter for the VEM
    control <- list(threshold = 1e-4, maxIter = 200, fixPointIter = 5, trace = FALSE)

    ## VEM Estimation on each element of the collection
    collection$estimate(control, mc.cores, TRUE)

    collection$smooth_ICL("none", control, mc.cores, NA, TRUE)

    collection$smooth_ICL("forward", control, mc.cores, NA, TRUE)

    collection$smooth_ICL("backward", control, mc.cores, NA, TRUE)

    collection$smooth_ICL("both", control, mc.cores, 2, TRUE)

    expect_is(collection, "missSBM_collection")

    expect_equal(which.min(collection$ICL), 3)

#  }
})
