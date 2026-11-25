## Issue #49: time and peak memory of the phases of a fit, to measure the remaining dense n x n
## computations before and after their removal.
##
## Usage, from the root of the package:
##   OPENBLAS_NUM_THREADS=1 Rscript inst/simus/issue49-memory/baseline.R <label> [cores]
## writes inst/simus/issue49-memory/results/<label>.rds
args  <- commandArgs(TRUE)
label <- args[1]
cores <- if (length(args) > 1) as.integer(args[2]) else 8
## the package: installed in SIMU_LIB if set (compiled with optimization, as for a user), from the
## sources otherwise (devtools::load_all() compiles the C++ code without optimization: timings of
## the C++ parts are then much slower than for an installed package)
lib <- Sys.getenv("SIMU_LIB", "")
if (nzchar(lib)) {
  suppressMessages(library(missSBM, lib.loc = lib))
  attach(asNamespace("missSBM"), name = "missSBM_internals", warn.conflicts = FALSE)
} else suppressMessages(devtools::load_all(quiet = TRUE))
future::plan("sequential")
library(parallel)

Q <- 3; alpha <- c(.25, .35, .40)
## connectivity scaled by 200 / n: a mean degree of about 26 whatever n
connect <- function(n) {pi <- matrix(.05, Q, Q); diag(pi) <- c(.6, .3, .15); pi * 200 / n}
sampling_param <- function(design, n) switch(design,
  "node" = .5, "dyad" = .5,
  "block-dyad" = {p <- matrix(.3, Q, Q); diag(p) <- c(.7, .5, .3); p},
  "double-standard" = c(.3, .7),
  "degree" = c(-2.4, .08))

## elapsed time of an expression, and the increase of the peak memory of the R heap (Mb) over
## what was in use before it (gc() reset before)
measure <- function(expr) {
  g0 <- gc(reset = TRUE)
  t <- system.time(res <- force(expr))[["elapsed"]]
  g <- gc()
  list(value = res, time = t, peak_mb = sum(g[, ncol(g)]) - sum(g0[, 2])) # "max used" - "used", Mb
}

grid <- expand.grid(n = c(2000, 4000), design = c("node", "block-dyad", "double-standard", "degree"),
                    stringsAsFactors = FALSE)

run <- function(k) {
  g <- grid[k, ]; set.seed(k)
  s0 <- sbm::sampleSimpleSBM(g$n, alpha, list(mean = connect(g$n)))
  A  <- suppressWarnings(observeNetwork(s0$networkData, g$design, sampling_param(g$design, g$n),
                                        clusters = s0$memberships))
  input_mb <- as.numeric(object.size(A)) / 2^20
  m_net  <- measure(partlyObservedNetwork$new(A))                       # point 5
  pon    <- m_net$value
  m_clus <- measure(pon$clustering(Q))                                  # points 2 and 3
  cl     <- m_clus$value[[1]]
  prof   <- tempfile()
  Rprof(prof, interval = .02)
  m_fit  <- measure({fit <- suppressWarnings(missSBM_fit$new(pon, g$design, cl))
                     fit$doVEM(list(threshold = 1e-2, maxIter = 50, fixPointIter = 3, trace = FALSE))
                     fit})
  Rprof(NULL)
  m_icl  <- measure(m_fit$value$ICL)
  by_total <- summaryRprof(prof)$by.total
  share <- function(pattern) {
    i <- grepl(pattern, rownames(by_total))
    if (any(i)) max(by_total$total.pct[i]) else 0
  }
  data.frame(g, input_mb = input_mb, observed = mean(!is.na(A)),
             net_time = m_net$time, net_peak = m_net$peak_mb,
             clustering_time = m_clus$time, clustering_peak = m_clus$peak_mb,
             fit_time = m_fit$time, fit_peak = m_fit$peak_mb, fit_iterations = nrow(m_fit$value$monitoring),
             icl_time = m_icl$time,
             pct_mask_dense = share("mask_dense_at_pattern"),
             pct_E_step = share("E_step"),
             pct_vLL = share("vLL_complete"), pct_M_step = share("update_parameters"))
}

res <- do.call(rbind, mclapply(seq_len(nrow(grid)), run, mc.cores = cores, mc.preschedule = FALSE))
attr(res, "commit") <- system("git rev-parse --short HEAD", intern = TRUE)
saveRDS(res, sprintf("inst/simus/issue49-memory/results/%s.rds", label))
print(res, digits = 3)
