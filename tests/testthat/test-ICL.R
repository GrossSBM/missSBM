context("ICL: observed data or completed network")

set.seed(1890)
N <- 120
sbm_sim <- sbm::sampleSimpleSBM(N, c(.4, .6), list(mean = matrix(c(.5, .1, .1, .3), 2, 2)))
A_node <- observeNetwork(sbm_sim$networkData, "node", .6)
A_dyad <- observeNetwork(sbm_sim$networkData, "dyad", .6)
ctrl <- function(...) missSBM_param(trace = FALSE, iterates = 0, polish = FALSE, ...)

test_that("the ICL built on the observed data integrates the missing dyads out", {
  for (design in c("node", "block-node", "dyad", "block-dyad")) {
    A <- if (design %in% c("node", "block-node")) A_node else A_dyad
    fit <- estimateMissSBM(A, 2, design, control = ctrl())$bestModel
    expect_equal(fit$criterion, "observed")
    expect_equal(fit$ICL, fit$ICLobserved)
    expect_equal(fit$ICLobserved,
                 -2 * (fit$fittedSampling$vExpec + fit$fittedSBM$vExpec) + fit$penalty)
    ## the completed network adds -2 times the expected log-likelihood of the imputed dyads,
    ## a positive term
    expect_gt(fit$ICLcompleted, fit$ICLobserved)
  }
})

test_that("the ICL can be the one of the completed network, as in the previous versions", {
  fit <- estimateMissSBM(A_dyad, 2, "dyad", control = ctrl(ICL = "completed"))$bestModel
  expect_equal(fit$criterion, "completed")
  expect_equal(fit$ICL, fit$ICLcompleted)
  expect_equal(fit$ICLcompleted,
               -2 * (fit$fittedSampling$vExpec + fit$fittedSBM$vExpec_corrected) + fit$penalty)
  expect_error(missSBM_param(ICL = "full"))
})

test_that("with a single block, both ICL are equal for node and block-node, dyad and block-dyad", {
  for (crit in c("observed", "completed")) {
    ICL <- function(A, design) estimateMissSBM(A, 1, design, control = ctrl(ICL = crit))$ICL
    expect_equal(ICL(A_node, "block-node"), ICL(A_node, "node"), tolerance = 1e-4)
    expect_equal(ICL(A_dyad, "block-dyad"), ICL(A_dyad, "dyad"), tolerance = 1e-4)
  }
})

test_that("degree and double-standard: the ICL built on the observed data includes the entropy of the imputations", {
  set.seed(2)
  A_degree <- observeNetwork(sbm_sim$networkData, "degree", c(-1, .05))
  A_ds     <- observeNetwork(sbm_sim$networkData, "double-standard", c(.4, .8))
  fit <- suppressWarnings(estimateMissSBM(A_degree, 2, "degree", control = ctrl()))$bestModel
  expect_true(is.finite(fit$ICLobserved))
  ## the expected log-likelihood of the missing dyads, at the imputations of the SBM, is the
  ## term that the completed network adds
  expect_equal(fit$fittedSBM$vExpec_imputed(fit$fittedSBM$imputation),
               fit$fittedSBM$vExpec_corrected - fit$fittedSBM$vExpec)
  expect_equal(fit$ICL, fit$ICLobserved)

  fit <- estimateMissSBM(A_ds, 2, "double-standard", control = ctrl())$bestModel
  expect_equal(fit$ICLobserved, -2 * (fit$vExpec + fit$entropyImputed) + fit$penalty)
  expect_equal(fit$ICLobserved, -2 * (fit$loglik - fit$fittedSBM$entropy) + fit$penalty)
})
