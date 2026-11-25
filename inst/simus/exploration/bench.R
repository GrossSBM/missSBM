## Fits each instance (see setup.R) with Q = 1 to 14 and several settings of the exploration:
##   none:    no exploration, no polish (VEM from the spectral initialization only)
##   polish:  node-swap polish only
##   default: 1 forward/backward pass + polish (the default of missSBM_param())
##   iter2, iter3: 2 or 3 passes + polish
##
## Usage (from the root of the package; see load-missSBM.R for SIMU_LIB):
##   [SIMU_LIB=<lib>] OPENBLAS_NUM_THREADS=1 Rscript inst/simus/exploration/bench.R [nb cores] [label] [settings]
## settings: a comma-separated subset of the above (all by default). Writes results/<label>.rds
## (label: the version by default).

args  <- commandArgs(TRUE)
cores <- if (length(args) > 0) as.integer(args[1]) else 1L
source(file.path("inst", "simus", "load-missSBM.R"))
label <- if (length(args) > 1) args[2] else simu_version
dir   <- file.path("inst", "simus", "exploration", "results")

I <- readRDS(file.path(dir, "instances.rds"))
settings <- list(none    = list(iterates = 0, polish = FALSE),
                 polish  = list(iterates = 0, polish = TRUE),
                 default = list(iterates = 1, polish = TRUE),
                 iter2   = list(iterates = 2, polish = TRUE),
                 iter3   = list(iterates = 3, polish = TRUE))
if (length(args) > 2) settings <- settings[strsplit(args[3], ",")[[1]]]

jobs <- do.call(rbind, lapply(names(I$inst), function(i)
  expand.grid(instance = i, design = I$inst[[i]]$designs, setting = names(settings), stringsAsFactors = FALSE)))
run <- function(k) {
  j <- jobs[k, ]
  ctrl <- c(settings[[j$setting]], trace = FALSE)
  set.seed(42)
  t <- system.time(col <- suppressWarnings(
    estimateMissSBM(I$inst[[j$instance]]$A, 1:14, j$design, control = do.call(missSBM_param, ctrl))))[["elapsed"]]
  list(job = j, time = t, ICL = unname(col$ICL), degenerate = unname(col$degenerate),
       cl = lapply(col$models, function(m) m$fittedSBM$memberships), best = col$bestModel$fittedSBM$memberships)
}
res <- parallel::mclapply(seq_len(nrow(jobs)), run, mc.cores = cores, mc.preschedule = FALSE)
saveRDS(list(label = label, version = simu_version, source = simu_source, cores = cores, jobs = jobs, res = res), file.path(dir, paste0(label, ".rds")))
cat(label, "done,", sum(sapply(res, inherits, "try-error")), "errors\n")
