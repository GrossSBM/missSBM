context("degree sampling: estimation of the sampling parameters and imputation")

## The sampling model of 'degree' depends on the degrees of the full network: its parameters can
## only be recovered if the degrees are computed on the observed part plus the imputations, and
## the imputations must be probabilities. Both were broken (see NEWS), while the other tests only
## checked the coherence of the fits.

simulate_degree_sampling <- function(seed, directed = FALSE) {
  set.seed(seed)
  n <- 200; Q <- 3
  pi <- matrix(.05, Q, Q); diag(pi) <- c(.6, .3, .15); pi[1, 2] <- pi[2, 1] <- .15
  sbm <- sbm::sampleSimpleSBM(n, c(.2, .3, .5), list(mean = pi), directed = directed)
  A <- sbm$networkData
  list(
    A     = A,
    P     = sbm$indMemberships %*% pi %*% t(sbm$indMemberships),
    A_obs = missSBM::observeNetwork(A, "degree", c(-3.5, .12))
  )
}

test_that("degree sampling: the sampling parameters are recovered", {
  data <- simulate_degree_sampling(1)
  fit <- estimateMissSBM(data$A_obs, vBlocks = 3, sampling = "degree",
                         control = missSBM_param(trace = FALSE, iterates = 0))$models[[1]]
  psi <- fit$fittedSampling$parameters
  ## truth: (-3.5, 0.12), with high-degree nodes more often observed
  expect_gt(psi[2], 0.05)
  expect_lt(psi[2], 0.2)
  expect_gt(psi[1], -4.5)
  expect_lt(psi[1], -2)
})

for (directed in c(FALSE, TRUE)) {
  test_that(paste("degree sampling: the imputations are probabilities close to the truth,", if (directed) "directed" else "undirected"), {
    data <- simulate_degree_sampling(2, directed)
    fit <- estimateMissSBM(data$A_obs, vBlocks = 3, sampling = "degree",
                           control = missSBM_param(trace = FALSE, iterates = 0))$models[[1]]
    imputed <- as.matrix(fit$imputedNetwork)
    missing <- is.na(data$A_obs) & row(imputed) != col(imputed)
    expect_true(all(imputed[missing] > 0 & imputed[missing] < 1))
    expect_lt(sqrt(mean((imputed[missing] - data$P[missing])^2)), 0.08)
    ## observed dyads are left unchanged
    observed <- !is.na(data$A_obs) & row(imputed) != col(imputed)
    expect_equal(unname(imputed[observed]), unname(data$A_obs[observed]))
  })
}
