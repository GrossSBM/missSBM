## Fits each instance (see setup.R) with Q = 1 to 14 and several settings of the exploration:
##   none:    no exploration, no polish (VEM from the spectral initialization only)
##   polish:  node-swap polish only
##   default: 1 forward/backward pass + polish (the default of missSBM_param())
##   iter2, iter3: 2 or 3 passes + polish
##
## Usage (from the root of the package):
##   [OPENBLAS_NUM_THREADS=1] Rscript inst/simus/exploration-blog/bench.R <path to missSBM source> <label> [ICL] [nb cores] [settings]
## ICL: "observed", "completed", or "" for the versions without this option; settings: a
## comma-separated subset of the above (all by default). Writes results/<label>.rds.

args  <- commandArgs(TRUE)
pkg   <- args[1]; label <- args[2]
crit  <- if (length(args) > 2) args[3] else ""
cores <- if (length(args) > 3) as.integer(args[4]) else 1L
dir   <- file.path("inst", "simus", "exploration-blog", "results")
suppressMessages(devtools::load_all(pkg, quiet = TRUE))
future::plan("sequential") # the jobs run in parallel instead

I <- readRDS(file.path(dir, "instances.rds"))
settings <- list(none    = list(iterates = 0, polish = FALSE),
                 polish  = list(iterates = 0, polish = TRUE),
                 default = list(iterates = 1, polish = TRUE),
                 iter2   = list(iterates = 2, polish = TRUE),
                 iter3   = list(iterates = 3, polish = TRUE))
if (length(args) > 4) settings <- settings[strsplit(args[5], ",")[[1]]]

jobs <- do.call(rbind, lapply(names(I$inst), function(i)
  expand.grid(instance = i, design = I$inst[[i]]$designs, setting = names(settings), stringsAsFactors = FALSE)))
run <- function(k) {
  j <- jobs[k, ]
  ctrl <- c(settings[[j$setting]], trace = FALSE)
  if (nzchar(crit)) ctrl$ICL <- crit
  set.seed(42)
  t <- system.time(col <- suppressWarnings(
    estimateMissSBM(I$inst[[j$instance]]$A, 1:14, j$design, control = do.call(missSBM_param, ctrl))))[["elapsed"]]
  list(job = j, time = t, ICL = unname(col$ICL), degenerate = unname(col$degenerate),
       cl = lapply(col$models, function(m) m$fittedSBM$memberships), best = col$bestModel$fittedSBM$memberships)
}
res <- parallel::mclapply(seq_len(nrow(jobs)), run, mc.cores = cores, mc.preschedule = FALSE)
saveRDS(list(label = label, crit = crit, jobs = jobs, res = res), file.path(dir, paste0(label, ".rds")))
cat(label, "done,", sum(sapply(res, inherits, "try-error")), "errors\n")
