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
N <- 100
Q <- 3
pi <- rep(1, Q)/Q           # block proportion
theta <- list(mean = diag(.45, Q, Q) + .05) # connectivity matrix
directed <- FALSE              # if the network is directed or not

sbm <- sbm::sampleSimpleSBM(N, pi, theta) # simulation of ad Bernoulli non-directed SBM
Z0  <- sbm$indMemberships

samplings <- list(
  list(name = "dyad", psi = 0.5, class = "dyadSampling_fit", k = log(N * (N-1)/2)),
  list(name = "node", psi = 0.5, class = "nodeSampling_fit", k = log(N)),
  list(name = "double-standard", psi =  c(.3, .6), class = "doubleStandardSampling_fit", k = log(N * (N-1)/2)),
  list(name = "block-node", psi = c(.3, .5, .7), class = "blockNodeSampling_fit", k = log(N)),
  list(name = "block-dyad", psi = psi <- matrix(.5,3,3) + diag(3)*.3, class = "blockDyadSampling_fit", k = log(N * (N-1)/2))
#  list(name = "degree", psi = c(-.05, .01), class = "degreeSampling_fit", k = log(N))
)

test_that("Consistency of sampling fit", {

  tol_truth <- 1e-2
  cat("Tested sampling:")
  for (sampling in samplings) {
    cat("\n -", sampling$name)

    ## sampled the network
    adjMatrix  <- missSBM::observeNetwork(sbm$networkData, sampling$name, sampling$psi, sbm$memberships)
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
    } else {
      expect_lt(error(fittedSampling$parameters, sampling$psi), tol_truth * 10 )
    }

  }
})
