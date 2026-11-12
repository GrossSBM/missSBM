## Tables comparing two runs of simulate.R (by default master and PR48), written to
## inst/simus/PR48-block-designs/results/summary.md
##
## Usage (from the root of the package):
##   Rscript inst/simus/PR48-block-designs/summarize.R [label before] [label after]

args   <- commandArgs(TRUE)
labels <- if (length(args) >= 2) args[1:2] else c("master", "PR48")
dir    <- file.path("inst", "simus", "PR48-block-designs", "results")
runs   <- lapply(labels, function(l) readRDS(file.path(dir, paste0(l, ".rds"))))
names(runs) <- labels

fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)
pair <- function(f, d = 2) paste(fmt(f(runs[[1]]), d), "→", fmt(f(runs[[2]]), d))

## ---------------------------------------------------------------------------------------
## experiment 1: grid
settings <- unique(runs[[1]]$grid[, c("design", "conn", "rate")])
settings$rate <- factor(settings$rate, levels = c("low", "medium", "high"))
settings <- settings[order(settings$design, settings$conn, settings$rate), ]

grid_row <- function(s) {
  sub <- function(run) with(run$grid, run$grid[design == s$design & conn == s$conn & rate == as.character(s$rate) & is.na(error), ])
  stat <- function(f) function(run) f(sub(run))
  data.frame(
    Design = s$design, Connectivity = s$conn, Rate = as.character(s$rate),
    `Observed dyads` = fmt(mean(sub(runs[[2]])$obs_rate)),
    `Q = 3 selected` = pair(stat(function(x) mean(x$Q_sel == 3))),
    `Mean selected Q` = pair(stat(function(x) mean(x$Q_sel)), 1),
    `ARI, selected Q` = pair(stat(function(x) mean(x$ARI_sel))),
    `ARI, Q = 3` = pair(stat(function(x) mean(x$ARI_Q))),
    `RMSE connectivity, Q = 3` = pair(stat(function(x) mean(x$rmse_conn)), 3),
    `RMSE sampling, Q = 3` = pair(stat(function(x) mean(x$rmse_samp)), 3),
    `Degenerate fits per collection` = pair(stat(function(x) mean(x$nb_degenerate)), 1),
    `Q = 3 selected, observed-data ICL` = fmt(mean(sub(runs[[2]])$Q_sel_obs == 3)),
    `ARI, observed-data ICL` = fmt(mean(sub(runs[[2]])$ARI_sel_obs)),
    check.names = FALSE
  )
}
tab_grid <- do.call(rbind, lapply(seq_len(nrow(settings)), function(k) grid_row(settings[k, ])))

## ---------------------------------------------------------------------------------------
## experiment 2: choice of the design by the ICL
truths <- c("node", "block-node", "dyad", "block-dyad")
design_row <- function(t) {
  sub <- function(run) run$design[run$design$truth == t & is.na(run$design$error), ]
  data.frame(
    `Simulated with` = t,
    `Observed dyads` = fmt(mean(sub(runs[[2]])$obs_rate)),
    `True design chosen` = pair(function(run) mean(sub(run)$chosen == t)),
    `Selected Q, MAR design` = pair(function(run) mean(sub(run)$Q_MAR), 1),
    `Selected Q, block design` = pair(function(run) mean(sub(run)$Q_block), 1),
    `ARI, MAR design` = pair(function(run) mean(sub(run)$ARI_MAR)),
    `ARI, block design` = pair(function(run) mean(sub(run)$ARI_block)),
    `True design chosen, observed-data ICL` = fmt(mean(sub(runs[[2]])$chosen_obs == t)),
    `Selected Q (MAR / block), observed-data ICL` = paste(fmt(mean(sub(runs[[2]])$Q_MAR_obs), 1), "/", fmt(mean(sub(runs[[2]])$Q_block_obs), 1)),
    check.names = FALSE
  )
}
tab_design <- do.call(rbind, lapply(truths, design_row))

## ---------------------------------------------------------------------------------------
errors <- sapply(runs, function(run) sum(!is.na(run$grid$error)) + sum(!is.na(run$design$error)))
md <- c(
  sprintf("Each cell: %s → %s, mean over %d replicates. Runs: %s.", labels[1], labels[2],
          max(runs[[1]]$grid$rep),
          paste(sprintf("%s at %s (%.0f min)", labels, sapply(runs, `[[`, "commit"), sapply(runs, `[[`, "elapsed")), collapse = ", ")),
  if (any(errors > 0)) sprintf("Fits that failed: %s.", paste(labels, errors, sep = " ", collapse = ", ")),
  sprintf("The last columns: %s only, models selected with the ICL built on the observed data (see README).", labels[2]),
  "",
  "### Experiment 1: block-node and block-dyad, Q = 1 to 6 explored, true Q = 3",
  "",
  knitr::kable(tab_grid, format = "pipe", align = "lllrrrrrrrrrr", row.names = FALSE),
  "",
  "### Experiment 2: choice of the design by the ICL (assortative network, medium rate)",
  "",
  knitr::kable(tab_design, format = "pipe", align = "lrrrrrrrr", row.names = FALSE)
)
writeLines(md, file.path(dir, "summary.md"))
cat(md, sep = "\n")
