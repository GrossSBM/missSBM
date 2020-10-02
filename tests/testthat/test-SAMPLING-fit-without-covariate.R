context("test network sampling fit with covariates (Classes nodeSampling_fit and dyadSampling_fit)")

source("utils_test.R")

set.seed(178303)
### SBM model
N <- 100
Q <- 3
pi <- rep(1, Q)/Q           # block proportion
theta <- diag(.45, Q, Q) + .05 # connectivity matrix
directed <- FALSE              # if the network is directed or not

sbm <- missSBM::simulate(N, pi, theta, directed) # simulation of ad Bernoulli non-directed SBM
Z0  <- missSBM:::clustering_indicator(sbm$memberships)

samplings <- list(
  list(name = "dyad", psi = 0.5, class = "dyadSampling_fit", k = log(N * (N-1)/2)),
  list(name = "node", psi = 0.5, class = "nodeSampling_fit", k = log(N)),
  list(name = "double-standard", psi =  c(.3, .6), class = "doubleStandardSampling_fit", k = log(N * (N-1)/2)),
  list(name = "block-node", psi = c(.3, .5, .7), class = "blockSampling_fit", k = log(N)),
  list(name = "degree", psi = c(-.05, .01), class = "degreeSampling_fit", k = log(N))
)

test_that("Consistency of sampling fit", {

  tol_truth <- 1e-2
  cat("Tested sampling:")
  for (sampling in samplings) {
    cat("\n -", sampling$name)

    ## sampled the network
    adjMatrix  <- missSBM::sample(sbm$netMatrix, sampling$name, sampling$psi, sbm$memberships)
    sampledNet <- missSBM:::sampledNetwork$new(adjMatrix)
    fittedSampling <- switch(
      sampling$name,
      "dyad"            = missSBM:::dyadSampling_fit$new(sampledNet),
      "node"            = missSBM:::nodeSampling_fit$new(sampledNet),
      "double-standard" = missSBM:::doubleStandardSampling_fit$new(sampledNet),
      "block-node"      = missSBM:::blockSampling_fit$new(sampledNet, Z0),
      "degree"          = missSBM:::degreeSampling_fit$new(sampledNet, Z0, sbm$connectParam)
    )

    expect_is(fittedSampling, sampling$class)
    expect_equal(fittedSampling$df, length(sampling$psi))
    expect_equal(fittedSampling$penalty, sampling$k * length(sampling$psi))
    expect_lt(fittedSampling$vExpec, 0)

    if (sampling$name %in% c("dyad", "node"))
      expect_lt(error(fittedSampling$parameters, sampling$psi), tol_truth)

  }
})
