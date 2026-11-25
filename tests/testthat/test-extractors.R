context("test the extractors of the fits and collections, for every design")

## every user-facing extractor runs on every design, and predict() (the imputed network) keeps
## the observed dyads, stays in [0, 1], and is symmetric for undirected networks
n <- 40; Q <- 2
connect     <- matrix(c(.6, .1, .1, .4), Q, Q)  # undirected
connect_dir <- matrix(c(.6, .05, .2, .4), Q, Q) # directed
params <- list("dyad" = .6, "node" = .7, "block-node" = c(.8, .5),
               "block-dyad" = matrix(c(.8, .4, .4, .6), Q, Q), "double-standard" = c(.4, .7),
               "degree" = c(-1, .1), "snowball" = c(2, .2))
ctrl <- missSBM_param(trace = FALSE, iterates = 0, polish = FALSE)

check_extractors <- function(A, sampling, directed, covariates = list(), useCov = TRUE) {
  control <- ctrl; control$useCov <- useCov
  collection <- suppressWarnings(estimateMissSBM(A, 1:Q, sampling, covariates = covariates, control = control))
  expect_output(print(collection))
  expect_output(summary(collection))
  for (type in c("icl", "elbo", "monitoring")) expect_s3_class(plot(collection, type), "ggplot")

  fit <- collection$models[[Q]]
  expect_output(print(fit), "missSBM-fit")
  expect_output(summary(fit), "missSBM-fit")
  expect_length(coef(fit, "mixture"), Q)
  expect_true(is.list(coef(fit, "connectivity")))
  expect_equal(coef(fit, "sampling"), fit$fittedSampling$parameters)
  if (useCov && length(covariates)) expect_length(coef(fit, "covariates"), length(covariates))
  for (type in c("imputed", "expected", "monitoring")) expect_s3_class(plot(fit, type), "ggplot")

  expect_equal(dim(fitted(fit)), dim(A))
  P <- as.matrix(predict(fit))
  expect_equal(P, as.matrix(fit$imputedNetwork))
  obs <- !is.na(A); diag(obs) <- FALSE
  expect_equal(P[obs], A[obs])
  off_diagonal <- P[row(P) != col(P)]
  expect_false(anyNA(off_diagonal))
  expect_true(all(off_diagonal >= 0 & off_diagonal <= 1))
  if (!directed) expect_true(isSymmetric(unname(P)))
  expect_true(is.finite(fit$ICL))
}

test_that("extractors are consistent for every design, undirected and directed", {
  for (directed in c(FALSE, TRUE)) {
    set.seed(1)
    net <- sbm::sampleSimpleSBM(n, c(.5, .5), list(mean = if (directed) connect_dir else connect), directed = directed)
    for (sampling in names(params)) {
      set.seed(2)
      A <- suppressWarnings(observeNetwork(net$networkData, sampling, params[[sampling]], clusters = net$memberships))
      check_extractors(A, sampling, directed)
    }
  }
})

test_that("extractors are consistent with covariates, in the SBM and in the sampling", {
  set.seed(3)
  net <- sbm::sampleSimpleSBM(n, c(.5, .5), list(mean = connect))
  x_node <- rnorm(n)
  x_dyad <- matrix(rnorm(n * n), n, n); x_dyad <- x_dyad + t(x_dyad)
  A <- observeNetwork(net$networkData, "covar-node", 1, covariates = list(x_node))
  for (useCov in c(TRUE, FALSE)) check_extractors(A, "covar-node", FALSE, list(x_node), useCov)
  A <- observeNetwork(net$networkData, "covar-dyad", .5, covariates = list(x_dyad))
  for (useCov in c(TRUE, FALSE)) check_extractors(A, "covar-dyad", FALSE, list(x_dyad), useCov)
})
