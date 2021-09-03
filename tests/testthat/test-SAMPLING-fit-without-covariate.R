context("test network sampling fit with covariates (Classes nodeSampling_fit and dyadSampling_fit)")

set.seed(321)

N_nocov <- 300
source("utils_test.R", local = TRUE)
N <- N_nocov

samplings <- list(
  list(name = "dyad", psi = 0.5, class = "dyadSampling_fit", k = log(N * (N-1)/2)),
  list(name = "node", psi = 0.5, class = "nodeSampling_fit", k = log(N)),
  list(name = "double-standard", psi =  c(.3, .6), class = "doubleStandardSampling_fit", k = log(N * (N-1)/2)),
  list(name = "block-node", psi = c(.3, .5, .7), class = "blockNodeSampling_fit", k = log(N)),
  list(name = "block-dyad", psi = psi <- matrix(.5,3,3) + diag(3)*.3, class = "blockDyadSampling_fit", k = log(N * (N-1)/2))
#  list(name = "degree", psi = c(-.05, .01), class = "degreeSampling_fit", k = log(N))
)

test_that("Consistency of sampling fit for undirected bernoulli withou covariate", {

  sampler_undirected_nocov$rNetwork(store = TRUE)
  Z0  <- sampler_undirected_nocov$indMemberships

  tol_truth <- 1e-1
  cat("Tested sampling:")
  for (sampling in samplings) {
    cat("\n -", sampling$name)

    ## sampled the network
    adjMatrix  <- missSBM::observeNetwork(sampler_undirected_nocov$networkData, sampling$name, sampling$psi, sampler_undirected_nocov$memberships)
    partlyObservedNet <- missSBM:::partlyObservedNetwork$new(adjMatrix)
    fittedSampling <- switch(
      sampling$name,
      "dyad"            = missSBM:::dyadSampling_fit$new(partlyObservedNet),
      "node"            = missSBM:::nodeSampling_fit$new(partlyObservedNet),
      "double-standard" = missSBM:::doubleStandardSampling_fit$new(partlyObservedNet),
      "block-node"      = missSBM:::blockNodeSampling_fit$new(partlyObservedNet, Z0),
      "block-dyad"      = missSBM:::blockDyadSampling_fit$new(partlyObservedNet, Z0),
      "degree"          = missSBM:::degreeSampling_fit$new(partlyObservedNet, Z0, sbm$connectParam$mean)
    )

    expect_is(fittedSampling, sampling$class)
    expect_equal(fittedSampling$df, length(sampling$psi))
    expect_equal(fittedSampling$penalty, sampling$k * length(sampling$psi))
    expect_lte(fittedSampling$vExpec, 0)

    if (sampling$name %in% c("dyad", "node")) {
      expect_lt(error(fittedSampling$parameters, sampling$psi), tol_truth)
    }

  }
})

samplings <- list(
  list(name = "dyad", psi = 0.5, class = "dyadSampling_fit", k = log(N * (N-1))),
  list(name = "node", psi = 0.5, class = "nodeSampling_fit", k = log(N)),
  list(name = "double-standard", psi =  c(.3, .6), class = "doubleStandardSampling_fit", k = log(N * (N-1))),
  list(name = "block-node", psi = c(.3, .5, .7), class = "blockNodeSampling_fit", k = log(N)),
  list(name = "block-dyad", psi = psi <- matrix(seq(.9, .1, -.1),3,3), class = "blockDyadSampling_fit", k = log(N * (N-1)))
#  list(name = "degree", psi = c(-.05, .01), class = "degreeSampling_fit", k = log(N))
)

test_that("Consistency of sampling fit for directed network, no covariates", {

  sampler_directed_nocov$rNetwork(store = TRUE)
  Z0  <- sampler_directed_nocov$indMemberships

  tol_truth <- 1e-1
  cat("Tested sampling:")
  for (sampling in samplings) {
    cat("\n -", sampling$name)

    ## sampled the network
    adjMatrix  <- missSBM::observeNetwork(sampler_directed_nocov$networkData, sampling$name, sampling$psi, sampler_directed_nocov$memberships)
    partlyObservedNet <- missSBM:::partlyObservedNetwork$new(adjMatrix)
    fittedSampling <- switch(
      sampling$name,
      "dyad"            = missSBM:::dyadSampling_fit$new(partlyObservedNet),
      "node"            = missSBM:::nodeSampling_fit$new(partlyObservedNet),
      "double-standard" = missSBM:::doubleStandardSampling_fit$new(partlyObservedNet),
      "block-node"      = missSBM:::blockNodeSampling_fit$new(partlyObservedNet, Z0),
      "block-dyad"      = missSBM:::blockDyadSampling_fit$new(partlyObservedNet, Z0),
      "degree"          = missSBM:::degreeSampling_fit$new(partlyObservedNet, Z0, sbm$connectParam$mean)
    )

    expect_is(fittedSampling, sampling$class)
    expect_equal(fittedSampling$df, length(sampling$psi))
    expect_equal(fittedSampling$penalty, sampling$k * length(sampling$psi))
    expect_lte(fittedSampling$vExpec, 0)

    if (sampling$name %in% c("dyad", "node")) {
      expect_lt(error(fittedSampling$parameters, sampling$psi), tol_truth)
    }

  }
})
