## Table of inst/simus/performance/README.md, written to results/summary.md: two runs of
## measure.R (by default v1.1.0 and v1.2.0).
##
## Usage (from the root of the package):
##   Rscript inst/simus/performance/summarize.R [label before] [label after]
args <- commandArgs(TRUE)
labels <- if (length(args) == 2) args else c("v1.1.0", "v1.2.0")
a <- readRDS(sprintf("inst/simus/performance/results/%s.rds", labels[1]))
b <- readRDS(sprintf("inst/simus/performance/results/%s.rds", labels[2]))
f <- function(t0, m0, t1, m1) sprintf("%.2f s, +%.0f Mb → %.2f s, +%.0f Mb", t0, m0, t1, m1)
tab <- data.frame(n = a$n, Design = a$design,
  `partlyObservedNetwork$new()` = f(a$net_time, a$net_peak, b$net_time, b$net_peak),
  `clustering()` = f(a$clustering_time, a$clustering_peak, b$clustering_time, b$clustering_peak),
  `Fit (Q = 3)` = f(a$fit_time, a$fit_peak, b$fit_time, b$fit_peak),
  `Total time` = sprintf("%.1f s → %.1f s", a$net_time + a$clustering_time + a$fit_time,
                         b$net_time + b$clustering_time + b$fit_time),
  check.names = FALSE)
sink("inst/simus/performance/results/summary.md")
cat(sprintf("%s (%s) → %s (%s). Peaks: increase of the R heap during the phase.\n\n",
            labels[1], attr(a, "version"), labels[2], attr(b, "version")))
print(knitr::kable(tab[order(tab$Design, tab$n), ], row.names = FALSE))
sink()
