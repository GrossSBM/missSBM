## Table of inst/simus/issue49-memory/README.md: before (baseline.rds) and after (final.rds)
a <- readRDS("inst/simus/issue49-memory/results/baseline.rds")
b <- readRDS("inst/simus/issue49-memory/results/final.rds")
f <- function(t0, m0, t1, m1) sprintf("%.2f s, +%.0f Mb → %.2f s, +%.0f Mb", t0, m0, t1, m1)
tab <- data.frame(n = a$n, Design = a$design,
  `partlyObservedNetwork$new()` = f(a$net_time, a$net_peak, b$net_time, b$net_peak),
  `clustering()` = f(a$clustering_time, a$clustering_peak, b$clustering_time, b$clustering_peak),
  `Fit (Q = 3)` = f(a$fit_time, a$fit_peak, b$fit_time, b$fit_peak),
  `Total time` = sprintf("%.1f s → %.1f s", a$net_time + a$clustering_time + a$fit_time,
                         b$net_time + b$clustering_time + b$fit_time),
  check.names = FALSE)
sink("inst/simus/issue49-memory/results/summary.md")
cat(sprintf("Before (commit %s) → after (commit %s). Peaks: increase of the R heap during the phase.\n\n",
            attr(a, "commit"), attr(b, "commit")))
print(knitr::kable(tab[order(tab$Design, tab$n), ], row.names = FALSE))
sink()
