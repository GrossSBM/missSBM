context("test missSBM-fit without covariate")

Q <- 3
N_nocov <- 100
N <- N_nocov
source("utils_test.R", local = TRUE)

sampler_undirected_nocov$rNetwork(store = TRUE)

samplings <- list(
  list(name = "dyad", psi = 0.5, class = "dyadSampling_fit"),
  list(name = "node", psi = 0.5, class = "nodeSampling_fit"),
  list(name = "double-standard", psi =  c(.3, .6), class = "doubleStandardSampling_fit"),
  list(name = "block-node", psi = c(.3, .5, .7), class = "blockNodeSampling_fit"),
  list(name = "block-dyad", psi = psi <- matrix(.5,3,3) + diag(3)*.3, class = "blockDyadSampling_fit")
)

## control parameter for the VEM
control <- list(threshold = 1e-3, maxIter = 100, fixPointIter = 5, trace = FALSE)

test_that("missSBM-fit works and is consistent for all samplings", {

  ## Consistency
  tol_truth <- .4
  tol_ARI   <- .8

  cat("Tested sampling:")
  for (sampling in samplings) {
    cat("\n -", sampling$name)

    ## sampled the network
    adjMatrix  <- missSBM::observeNetwork(sampler_undirected_nocov$networkData, sampling$name, sampling$psi, sampler_undirected_nocov$memberships)
    myNet <- missSBM:::partlyObservedNetwork$new(adjMatrix)
    cl <- myNet$clustering(1: (2*Q) )

    ## Perform inference
    missSBM <- missSBM:::missSBM_fit$new(myNet, sampling$name, cl[[Q]], FALSE)
    out <- missSBM$doVEM(control)

    ## Sanity check
    expect_true(inherits(missSBM, "missSBM_fit"))
    expect_true(inherits(missSBM$fittedSBM, "SimpleSBM_fit"))
    expect_is(missSBM$fittedSampling, sampling$class)
    expect_equal(out, missSBM$monitoring)

    ## Optimization success
    expect_gte(diff(range(out$elbo, na.rm = TRUE)), 0)

    if (ARI(missSBM$fittedSBM$memberships, sampler_undirected_nocov$memberships) > tol_ARI) {
      ## SBM: parameters estimation
      expect_lt(error(missSBM$fittedSBM$blockProp, sampler_undirected_nocov$blockProp, sort = TRUE), tol_truth)

      ## sampling design: parameters estimation
      expect_lt(error(missSBM$fittedSampling$parameters, sampling$psi, sort = TRUE), tol_truth)

      ## clustering
      expect_gt(ARI(missSBM$fittedSBM$memberships, sampler_undirected_nocov$memberships), tol_ARI)
    }
  }

})

