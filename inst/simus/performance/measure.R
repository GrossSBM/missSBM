## Time and peak memory of the phases of a fit (the network object, the initial clustering, the
## VEM at Q = 3, the ICL) on large sparse networks, for each family of designs.
##
## Usage (from the root of the package; see load-missSBM.R for SIMU_LIB, needed for timings):
##   SIMU_LIB=<lib> OPENBLAS_NUM_THREADS=1 Rscript inst/simus/performance/measure.R [cores] [label]
## writes inst/simus/performance/results/<label>.rds (label: the version by default). Fewer jobs
## than cores leave the measures undisturbed.
args  <- commandArgs(TRUE)
cores <- if (length(args) > 0) as.integer(args[1]) else 8
source(file.path("inst", "simus", "load-missSBM.R"))
label <- if (length(args) > 1) args[2] else simu_version
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
  m_net  <- measure(partlyObservedNetwork$new(A))
  pon    <- m_net$value
  m_clus <- measure(pon$clustering(Q))
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
             pct_E_step = share("E_step"), pct_M_step = share("update_parameters"))
}

res <- do.call(rbind, mclapply(seq_len(nrow(grid)), run, mc.cores = cores, mc.preschedule = FALSE))
attr(res, "version") <- paste(simu_version, simu_source)
dir.create("inst/simus/performance/results", showWarnings = FALSE)
saveRDS(res, sprintf("inst/simus/performance/results/%s.rds", label))
print(res, digits = 3)
