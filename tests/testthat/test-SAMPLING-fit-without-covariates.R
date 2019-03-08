context("test network sampling fit with covariates (Classes nodeSampling_fit and dyadSampling_fit)")

error <- function(beta1, beta2, sort = FALSE) {
  if (sort)
    err <- sum((sort(beta1) - sort(beta2))^2)/length(beta2)
  else
    err <- sum((beta1 - beta2)^2)/length(beta2)
  err
}

set.seed(178303)
### SBM model
N <- 400
Q <- 5
alpha <- rep(1,Q)/Q       # mixture parameter
pi <- diag(.45,Q) + .05   # connectivity matrix
directed <- FALSE         # if the network is directed or not

sbm <- simulateSBM(N, alpha, pi, directed) # simulation of ad Bernoulli non-directed SBM
Z0  <- missSBM:::clustering_indicator(sbm$memberships)

samplings <- list(
  list(name = "dyad", psi = 0.3, class = "dyadSampling_fit", k = log(N * (N-1)/2)),
  list(name = "node", psi = 0.2, class = "nodeSampling_fit", k = log(N)),
  list(name = "double-standard", psi =  c(.3, .6), class = "doubleStandardSampling_fit", k = log(N * (N-1)/2)),
  list(name = "block-node", psi = c(.1, .3, .2, .5, .7), class = "blockSampling_fit", k = log(N)),
  list(name = "degree", psi = c(-.05, .01), class = "degreeSampling_fit", k = log(N))
)

test_that("Consistency of sampling fit", {

  tol_truth <- 1e-2
  cat("Tested sampling:")
  for (sampling in samplings) {
    cat("\n -", sampling$name)

    ## sampled the network
    sampledNet <- sampleNetwork(sbm$adjMatrix, sampling$name, sampling$psi, sbm$memberships)

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
