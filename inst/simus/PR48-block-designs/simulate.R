## Simulations for PR #48: block-node and block-dyad samplings, SBM part fitted on the observed
## dyads (the missing dyads are integrated out) instead of the observed dyads plus the
## imputations of the missing ones.
##
## Two experiments, run with a given version of missSBM (source tree), the same data for all
## versions (the seeds only depend on the setting and the replicate):
##   - "grid": data simulated and fitted with block-node or block-dyad, 3 observation rates x
##     2 connectivity structures, Q = 1:6 explored: selected Q, ARI, RMSE of the connectivity
##     and of the sampling parameters, degenerate fits.
##   - "design": data simulated under node / block-node / degree (resp. dyad / block-dyad /
##     double-standard), fitted with the three designs of the family: does the ICL pick the
##     design that generated the data?
##
## Usage (from the root of the package):
##   [SIMU_ICL=observed|completed] Rscript inst/simus/PR48-block-designs/simulate.R <path to missSBM source> <label> [nb cores]
## writes inst/simus/PR48-block-designs/results/<label>.rds. SIMU_ICL sets the ICL used by the
## exploration and the selection (missSBM_param(ICL = ...)), for the versions that have it.

args  <- commandArgs(TRUE)
pkg   <- args[1]
label <- args[2]
cores <- if (length(args) > 2) as.integer(args[3]) else 1L
out_dir <- file.path("inst", "simus", "PR48-block-designs", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

suppressMessages({
  devtools::load_all(pkg, quiet = TRUE)
  library(sbm)
})

## ---------------------------------------------------------------------------------------
## settings
n      <- 200
Q      <- 3
alpha  <- c(.25, .35, .40)
vBlocks <- 1:6
nb_rep <- as.integer(Sys.getenv("SIMU_NB_REP", 20)) # fewer replicates for a quick try
ICL    <- Sys.getenv("SIMU_ICL", "")
control <- if (nzchar(ICL)) missSBM_param(trace = FALSE, ICL = ICL) else missSBM_param(trace = FALSE)

connectivity <- list(
  assortative    = matrix(c(.60, .05, .05,
                            .05, .30, .05,
                            .05, .05, .15), Q, Q),
  disassortative = matrix(c(.05, .40, .25,
                            .40, .10, .15,
                            .25, .15, .05), Q, Q)
)

## block-node: probability that a node of each block is observed
rho <- list(
  low    = c(.50, .30, .10),
  medium = c(.80, .50, .20),
  high   = c(.95, .85, .70)
)
## block-dyad: probability that a dyad between two blocks is observed
make_psi <- function(diag, off) { psi <- matrix(off, Q, Q); diag(psi) <- diag; psi }
psi <- list(
  low    = make_psi(c(.40, .30, .20), .10),
  medium = make_psi(c(.70, .50, .30), .30),
  high   = make_psi(c(.95, .90, .85), .80)
)
## experiment 2 only (medium rate): degree, P(node i observed) = logistic(a + b * degree of i),
## and double-standard, P(dyad observed | 0) and P(dyad observed | 1)
degree_param <- c(-2.4, .08)
double_standard_param <- c(.3, .7)

## ---------------------------------------------------------------------------------------
## helpers
sample_net <- function(conn, seed) {
  set.seed(seed)
  sbm::sampleSimpleSBM(n, alpha, list(mean = conn))
}

observe <- function(net, design, rate) {
  switch(design,
    "block-node" = observeNetwork(net$networkData, "block-node", rho[[rate]], clusters = net$memberships),
    "block-dyad" = observeNetwork(net$networkData, "block-dyad", psi[[rate]], clusters = net$memberships),
    "node"       = observeNetwork(net$networkData, "node", mean(rho[[rate]][net$memberships])),
    "dyad"       = observeNetwork(net$networkData, "dyad", {
      z <- net$memberships; P <- psi[[rate]][z, z]; mean(P[upper.tri(P)])
    }),
    "degree"          = observeNetwork(net$networkData, "degree", degree_param),
    "double-standard" = observeNetwork(net$networkData, "double-standard", double_standard_param)
  )
}

upper <- function(M) M[upper.tri(M)]

## metrics of the fit with Q blocks (the true number), permutation-free: the fitted
## probabilities of each dyad (resp. node) are compared to the true ones
metrics_at_Q <- function(collection, net, design, rate) {
  nbBlocks <- sapply(collection$models, function(m) m$fittedSBM$nbBlocks)
  fit <- collection$models[[which(nbBlocks == Q)[1]]]
  tau <- fit$fittedSBM$probMemberships
  z   <- net$memberships
  conn <- net$connectParam$mean
  P_hat <- tau %*% fit$fittedSBM$connectParam$mean %*% t(tau)
  samp <- fit$fittedSampling$parameters
  rmse_samp <- switch(design,
    "block-node" = sqrt(mean((tau %*% samp - rho[[rate]][z])^2)),
    "block-dyad" = sqrt(mean(upper(tau %*% samp %*% t(tau) - psi[[rate]][z, z])^2)),
    NA_real_
  )
  data.frame(
    ARI_Q     = aricode::ARI(fit$fittedSBM$memberships, z),
    rmse_conn = sqrt(mean(upper(P_hat - conn[z, z])^2)),
    rmse_samp = rmse_samp,
    degenerate_Q = length(unique(fit$fittedSBM$memberships)) < Q
  )
}

fit_design <- function(A, design, seed) {
  set.seed(seed)
  time <- system.time(
    collection <- suppressWarnings(estimateMissSBM(A, vBlocks, design, control = control))
  )[["elapsed"]]
  list(collection = collection, time = time)
}

## ---------------------------------------------------------------------------------------
## experiment 1: grid
grid <- expand.grid(rep = seq_len(nb_rep), rate = names(rho), conn = names(connectivity),
                    design = c("block-node", "block-dyad"), stringsAsFactors = FALSE)
grid$seed <- 1000 * rep(seq_len(nrow(grid) / nb_rep), each = nb_rep) + grid$rep # one block of seeds per setting

run_grid <- function(k) {
  g <- grid[k, ]
  res <- tryCatch({
    net <- sample_net(connectivity[[g$conn]], g$seed)
    set.seed(g$seed + 1)
    A <- observe(net, g$design, g$rate)
    f <- fit_design(A, g$design, g$seed + 2)
    best <- f$collection$bestModel
    cbind(
      data.frame(obs_rate = mean(!is.na(upper(A))),
                 Q_sel    = best$fittedSBM$nbBlocks,
                 ARI_sel  = aricode::ARI(best$fittedSBM$memberships, net$memberships),
                 nb_degenerate = sum(sapply(f$collection$models, function(m)
                   length(unique(m$fittedSBM$memberships)) < m$fittedSBM$nbBlocks)),
                 time = f$time, error = NA_character_),
      metrics_at_Q(f$collection, net, g$design, g$rate)
    )
  }, error = function(e) data.frame(error = conditionMessage(e)))
  cbind(g, res)
}

## ---------------------------------------------------------------------------------------
## experiment 2: choice of the design by the ICL (assortative network, medium rate)
families <- list(node = c("node", "block-node", "degree"), dyad = c("dyad", "block-dyad", "double-standard"))
design_cases <- expand.grid(rep = seq_len(nb_rep), truth = unlist(families), stringsAsFactors = FALSE)
design_cases$family <- ifelse(design_cases$truth %in% families$node, "node", "dyad")
design_cases$seed   <- 100000 + 1000 * match(design_cases$truth, unique(design_cases$truth)) + design_cases$rep

run_design <- function(k) {
  d <- design_cases[k, ]
  candidates <- families[[d$family]]
  res <- tryCatch({
    net <- sample_net(connectivity$assortative, d$seed)
    set.seed(d$seed + 1)
    A <- observe(net, d$truth, "medium")
    fits <- lapply(candidates, function(design) fit_design(A, design, d$seed + 2)$collection$bestModel)
    ICL <- sapply(fits, function(m) m$ICL)
    ## one column per candidate, named by its rank in the family (1: MAR, 2: block, 3: other MNAR)
    per_candidate <- function(prefix, values) setNames(as.list(values), paste0(prefix, seq_along(values)))
    data.frame(
      obs_rate = mean(!is.na(upper(A))),
      chosen   = candidates[which.min(ICL)],
      per_candidate("Q_", sapply(fits, function(m) m$fittedSBM$nbBlocks)),
      per_candidate("ARI_", sapply(fits, function(m) aricode::ARI(m$fittedSBM$memberships, net$memberships))),
      per_candidate("ICL_", ICL),
      error = NA_character_
    )
  }, error = function(e) data.frame(error = conditionMessage(e)))
  cbind(d, res)
}

## ---------------------------------------------------------------------------------------
## run
t0 <- Sys.time()
res_grid   <- parallel::mclapply(seq_len(nrow(grid)), run_grid, mc.cores = cores, mc.preschedule = FALSE)
res_design <- parallel::mclapply(seq_len(nrow(design_cases)), run_design, mc.cores = cores, mc.preschedule = FALSE)
bind <- function(l) {
  cols <- unique(unlist(lapply(l, names)))
  do.call(rbind, lapply(l, function(x) { x[setdiff(cols, names(x))] <- NA; x[cols] }))
}
results <- list(
  label   = label,
  ICL     = if (nzchar(ICL)) ICL else "default",
  version = as.character(utils::packageVersion("missSBM", lib.loc = NULL)),
  commit  = tryCatch(system2("git", c("-C", pkg, "rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) NA),
  date    = Sys.time(),
  elapsed = as.numeric(difftime(Sys.time(), t0, units = "mins")),
  grid    = bind(res_grid),
  design  = bind(res_design)
)
saveRDS(results, file.path(out_dir, paste0(label, ".rds")))
cat(sprintf("%s (%s): %d grid fits, %d design cases, %.1f min, %d errors\n", label, results$commit,
            nrow(results$grid), nrow(results$design), results$elapsed,
            sum(!is.na(results$grid$error)) + sum(!is.na(results$design$error))))
