## The penalty of the ICL counts the observed dyads for the connectivity, on which it is
## estimated, rather than all the possible dyads. Both penalties are computed on the same fits:
## the one counting all the dyads adds (nbConnectParam + nbCovariates) log(nbDyads / nbObservedDyads).
## Networks near the threshold of detectability, where the two can select different models.
##
## Usage (from the root of the package; see load-missSBM.R for SIMU_LIB):
##   [SIMU_LIB=<lib>] OPENBLAS_NUM_THREADS=1 Rscript inst/simus/icl-penalty/simulate.R [cores]
args  <- commandArgs(TRUE)
cores <- if (length(args) > 0) as.integer(args[1]) else 16
source(file.path("inst", "simus", "load-missSBM.R"))
library(parallel)

nets <- list(
  ER    = list(Q = 1, alpha = c(.5, .5), pi = matrix(.1, 2, 2)), # two identical blocks: Erdos-Renyi
  sbm2a = list(Q = 2, alpha = c(.5, .5), pi = matrix(c(.20, .05, .05, .20), 2)),
  sbm2b = list(Q = 2, alpha = c(.5, .5), pi = matrix(c(.25, .05, .05, .25), 2)),
  sbm3a = list(Q = 3, alpha = rep(1/3, 3), pi = matrix(.05, 3, 3) + diag(.2, 3)),
  sbm3b = list(Q = 3, alpha = rep(1/3, 3), pi = matrix(.05, 3, 3) + diag(.3, 3))
)
rates <- list(dyad = c(low = .15, medium = .35), node = c(low = .35, medium = .6))
grid <- expand.grid(net = names(nets), n = c(100, 200), design = c("dyad", "node"),
                    rate = c("low", "medium"), rep = 1:20, stringsAsFactors = FALSE)

run <- function(k) {
  g <- grid[k, ]; set.seed(k)
  m <- nets[[g$net]]
  s0 <- sbm::sampleSimpleSBM(g$n, m$alpha, list(mean = m$pi))
  A  <- observeNetwork(s0$networkData, g$design, rates[[g$design]][[g$rate]])
  col <- estimateMissSBM(A, 1:5, g$design, control = missSBM_param(trace = FALSE))
  icl_observed <- sapply(col$models, function(f) f$ICL)
  icl_all      <- sapply(col$models, function(f) {
    sbm <- f$fittedSBM
    f$ICL + (sbm$nbConnectParam + sbm$nbCovariates) * log(sbm$nbDyads / sbm$nbObservedDyads)
  })
  data.frame(g, Qtrue = m$Q, observed_dyads = mean(!is.na(A[upper.tri(A)])),
             Q_all = which.min(icl_all), Q_observed = which.min(icl_observed))
}

t0  <- Sys.time()
res <- do.call(rbind, mclapply(seq_len(nrow(grid)), run, mc.cores = cores))
attr(res, "version") <- paste(simu_version, simu_source)
saveRDS(res, "inst/simus/icl-penalty/results/results.rds")
cat(nrow(res), "networks,", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "\n")
