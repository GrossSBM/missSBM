context("ICL built on the observed data")

set.seed(1890)
N <- 120
sbm_sim <- sbm::sampleSimpleSBM(N, c(.4, .6), list(mean = matrix(c(.5, .1, .1, .3), 2, 2)))
A_node <- observeNetwork(sbm_sim$networkData, "node", .6)
A_dyad <- observeNetwork(sbm_sim$networkData, "dyad", .6)
ctrl <- function(...) missSBM_param(trace = FALSE, iterates = 0, polish = FALSE, ...)

## the MAP memberships of a fit, and the log-likelihood of the observed dyads and of the blocks
## at given memberships
map_Z <- function(fit) {
  tau <- fit$fittedSBM$probMemberships
  Z <- matrix(0, nrow(tau), ncol(tau)); Z[cbind(seq_len(nrow(tau)), max.col(tau))] <- 1; Z
}
vLL_obs <- function(fit, Z) {
  priv <- environment(fit$fittedSBM$get_state)$private
  missSBM:::vLL_complete_sparse_bernoulli_nocovariate(priv$Y, priv$R, Z, fit$fittedSBM$connectParam$mean, fit$fittedSBM$blockProp)
}

test_that("the ICL is -2 times the bound at the MAP memberships, plus the penalty", {
  for (design in c("node", "block-node", "dyad", "block-dyad")) {
    A <- if (design %in% c("node", "block-node")) A_node else A_dyad
    fit <- estimateMissSBM(A, 2, design, control = ctrl())$bestModel
    Z <- map_Z(fit)
    expect_equal(fit$ICL, -2 * (fit$fittedSampling$vExpec_at(Z) + vLL_obs(fit, Z)) + fit$penalty)
    ## near-hard memberships: close to the bound without the entropy of the blocks
    expect_equal(fit$ICL, -2 * (fit$loglik - fit$fittedSBM$entropy) + fit$penalty, tolerance = 1e-3)
    ## no variational distribution on the missing dyads
    expect_equal(fit$entropyImputed, 0)
  }
  ## block-node: the sampling term at the MAP memberships, node by node
  fit <- estimateMissSBM(A_node, 2, "block-node", control = ctrl())$bestModel
  rho <- fit$fittedSampling$parameters; z <- max.col(fit$fittedSBM$probMemberships)
  obs <- missSBM:::partlyObservedNetwork$new(A_node)$observedNodes
  expect_equal(fit$fittedSampling$vExpec_at(map_Z(fit)), sum(log(rho[z[obs]])) + sum(log(1 - rho[z[!obs]])))
})

test_that("with a single block, the ICL is equal for node and block-node, dyad and block-dyad", {
  ICL <- function(A, design) estimateMissSBM(A, 1, design, control = ctrl())$ICL
  expect_equal(ICL(A_node, "block-node"), ICL(A_node, "node"), tolerance = 1e-4)
  expect_equal(ICL(A_dyad, "block-dyad"), ICL(A_dyad, "dyad"), tolerance = 1e-4)
})

test_that("degree: the ICL includes the expected log-likelihood and the entropy of the imputations", {
  set.seed(2)
  A_degree <- observeNetwork(sbm_sim$networkData, "degree", c(-1, .05))
  fit <- suppressWarnings(estimateMissSBM(A_degree, 2, "degree", control = ctrl()))$bestModel
  expect_true(is.finite(fit$ICL))
  expect_gt(fit$entropyImputed, 0)
  nu <- fit$imputations # the imputations of the sampling
  Z  <- map_Z(fit)
  expect_equal(fit$ICL, -2 * (fit$fittedSampling$vExpec + fit$fittedSBM$vExpec_at(Z) +
                                fit$fittedSBM$vExpec_imputed(nu, Z) + fit$entropyImputed) + fit$penalty)
})

test_that("double-standard: the missing dyads are integrated out given the blocks", {
  set.seed(2)
  A_ds <- observeNetwork(sbm_sim$networkData, "double-standard", c(.4, .8))
  fit <- estimateMissSBM(A_ds, 2, "double-standard", control = ctrl())$bestModel
  ## as for the designs ignorable given the blocks: no entropy of the imputations
  Z <- map_Z(fit)
  expect_equal(fit$ICL, -2 * (fit$fittedSampling$vExpec_at(Z) + vLL_obs(fit, Z)) + fit$penalty)
})

test_that("the penalty counts the observed dyads for the connectivity", {
  for (directed in c(FALSE, TRUE)) {
    set.seed(7)
    sim <- sbm::sampleSimpleSBM(N, c(.4, .6), list(mean = matrix(c(.5, .1, .1, .3), 2, 2)), directed = directed)
    A   <- observeNetwork(sim$networkData, "dyad", .3)
    fit <- estimateMissSBM(A, 2, "dyad", control = ctrl())$bestModel
    nb_obs <- if (directed) sum(!is.na(A)) - N else sum(!is.na(A[upper.tri(A)]))
    nb_all <- if (directed) N * (N - 1) else N * (N - 1) / 2
    nb_connect <- if (directed) 4 else 3
    expect_equal(fit$fittedSBM$nbObservedDyads, nb_obs)
    ## the connectivity counts the observed dyads, the block proportions the nodes, the
    ## sampling parameter (dyad: 1) all the dyads
    expect_equal(fit$penalty, nb_connect * log(nb_obs) + log(N) + log(nb_all))
  }
  ## a fully observed network: all the dyads are observed
  fit <- estimateMissSBM(sbm_sim$networkData, 2, "dyad", control = ctrl())$bestModel
  expect_equal(fit$fittedSBM$nbObservedDyads, N * (N - 1) / 2)
})
