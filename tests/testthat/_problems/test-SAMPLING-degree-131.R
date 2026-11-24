# Extracted from test-SAMPLING-degree.R:131

# prequel ----------------------------------------------------------------------
context("degree sampling: estimation of the sampling parameters and imputation")
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
for (directed in c(FALSE, TRUE)) {
  test_that(paste("degree sampling: sparse updates match the dense formulas,", if (directed) "directed" else "undirected"), {
    set.seed(3)
    n <- 60
    A <- matrix(rbinom(n * n, 1, .15), n, n)
    if (!directed) A[lower.tri(A)] <- t(A)[lower.tri(A)]
    diag(A) <- 0
    A[sample(n, 25), ] <- NA
    if (!directed) A[lower.tri(A)] <- t(A)[lower.tri(A)]
    pon <- missSBM:::partlyObservedNetwork$new(A)
    Z <- missSBM:::clustering_indicator(sample(1:2, n, replace = TRUE))
    PI0 <- matrix(c(.3, .05, .05, .2), 2, 2)
    ## the initial glm of psi separates perfectly here (unobserved nodes have no observed edge)
    fit <- suppressWarnings(missSBM:::degreeSampling_fit$new(pon, Z, PI0))
    S <- pon$samplingMatrixBar
    nu <- S; nu@x <- runif(length(S@x), .05, .95)
    PI <- S; PI@x <- runif(length(S@x), .05, .95)

    ## dense reference, as in missSBM <= 1.1.1 with the degree fixes
    NAs <- as.matrix(S) != 0
    Y   <- as.matrix(pon$networkData)
    nuD <- as.matrix(nu)
    if (!directed) { NAs <- NAs | t(NAs); Y <- Y + t(Y); nuD <- nuD + t(nuD) }
    psi0 <- fit$parameters
    N_obs <- pon$observedNodes
    D  <- rowSums(Y + nuD)
    V  <- nuD; V[!NAs] <- NA
    D2 <- rowSums(V * (1 - V), na.rm = TRUE) + D^2
    ksi <- pmax(sqrt(psi0[1]^2 + psi0[2]^2 * D2 + 2 * psi0[1] * psi0[2] * D), 1e-6)
    hk <- missSBM:::h(ksi)
    C0 <- .5 * n - sum(!N_obs); s1 <- .5 * sum(D) - sum(D[!N_obs])
    denom <- 2 * (sum(hk) * sum(hk * D2) - sum(hk * D)^2)
    psi <- c((s1 * sum(hk * D) - C0 * sum(hk * D2)) / denom, (C0 * sum(hk * D) - s1 * sum(hk)) / denom)
    V[!NAs] <- 0
    Dij <- matrix(D, n, n) - V
    C  <- 2 * hk * (psi[1] * psi[2] + psi[2]^2 * (1 + Dij))
    eta <- C + t(C) - psi[2]
    ref <- missSBM:::check_boundaries(missSBM:::.logistic(missSBM:::.logit(as.matrix(PI)) + eta))[as.matrix(S) != 0]

    fit$update_parameters(nu)
    expect_equal(fit$parameters, psi, tolerance = 1e-12)
    expect_equal(fit$update_imputation(PI)@x, ref, tolerance = 1e-12)
  })
}
small_degree_fit <- function(seed, n = 60, maxIter = 100, threshold = 1e-6) {
  set.seed(seed)
  sim <- sbm::sampleSimpleSBM(n, c(.5, .5), list(mean = matrix(c(.4, .05, .05, .2), 2, 2)))
  A   <- missSBM::observeNetwork(sim$networkData, "degree", c(-1.5, .2))
  pon <- missSBM:::partlyObservedNetwork$new(A)
  ## the initial glm of (a, b) can separate perfectly on a small network
  fit <- suppressWarnings(missSBM:::missSBM_fit$new(pon, "degree", sim$memberships))
  fit$doVEM(list(threshold = threshold, maxIter = maxIter, fixPointIter = 3, trace = FALSE))
  list(fit = fit, pon = pon)
}

# test -------------------------------------------------------------------------
x <- small_degree_fit(11)
samp <- x$fit$fittedSampling
priv <- environment(samp$update_parameters)$private
nu <- environment(x$fit$doVEM)$private$nu
a <- samp$parameters[1]
b <- samp$parameters[2]
nuD <- as.matrix(nu)
nuD <- nuD + t(nuD)
exact <- sum(sapply(seq_along(priv$D_obs), function(i) {
    p <- 1
    for (v in nuD[i, nuD[i, ] > 0]) p <- c(p * (1 - v), 0) + c(0, p * v)
    D <- priv$D_obs[i] + seq_along(p) - 1
    x <- a + b * D
    sum(p * if (priv$N_obs[i]) plogis(x, log.p = TRUE) else plogis(-x, log.p = TRUE))
  }))
plug_in <- sum(ifelse(priv$N_obs, plogis(a + b * priv$D, log.p = TRUE), plogis(-(a + b * priv$D), log.p = TRUE)))
expect_lte(samp$vExpec, exact)
expect_gte(plug_in, exact)
expect_gt(max(priv$ksi), 1)
