## Tables comparing runs of simulate.R (by default master, PR48-completed and PR48-observed),
## written to inst/simus/PR48-block-designs/results/summary.md
##
## Usage (from the root of the package):
##   Rscript inst/simus/PR48-block-designs/summarize.R [labels of the runs, in order]

args   <- commandArgs(TRUE)
labels <- if (length(args) >= 2) args else c("master", "PR48-completed", "PR48-observed")
dir    <- file.path("inst", "simus", "PR48-block-designs", "results")
runs   <- lapply(labels, function(l) readRDS(file.path(dir, paste0(l, ".rds"))))
names(runs) <- labels

fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)
## one cell: the statistic f for each run, in order
chain <- function(f, d = 2) paste(sapply(runs, function(run) fmt(f(run), d)), collapse = " → ")

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
    `Observed dyads` = fmt(mean(sub(runs[[1]])$obs_rate)),
    `Q = 3 selected` = chain(stat(function(x) mean(x$Q_sel == 3))),
    `Mean selected Q` = chain(stat(function(x) mean(x$Q_sel)), 1),
    `ARI, selected Q` = chain(stat(function(x) mean(x$ARI_sel))),
    `ARI, Q = 3` = chain(stat(function(x) mean(x$ARI_Q))),
    `RMSE connectivity, Q = 3` = chain(stat(function(x) mean(x$rmse_conn)), 3),
    `RMSE sampling, Q = 3` = chain(stat(function(x) mean(x$rmse_samp)), 3),
    `Degenerate fits per collection` = chain(stat(function(x) mean(x$nb_degenerate)), 1),
    check.names = FALSE
  )
}
tab_grid <- do.call(rbind, lapply(seq_len(nrow(settings)), function(k) grid_row(settings[k, ])))

## ---------------------------------------------------------------------------------------
## experiment 2: choice of the design by the ICL, among the three designs of the family
families <- list(node = c("node", "block-node", "degree"), dyad = c("dyad", "block-dyad", "double-standard"))
design_row <- function(t) {
  candidates <- families[[if (t %in% families$node) "node" else "dyad"]]
  k <- match(t, candidates)
  sub <- function(run) run$design[run$design$truth == t & is.na(run$design$error), ]
  data.frame(
    `Simulated with` = t,
    `Observed dyads` = fmt(mean(sub(runs[[1]])$obs_rate)),
    `True design chosen` = chain(function(run) mean(sub(run)$chosen == t)),
    `Chosen instead (last run)` = {
      wrong <- sub(runs[[length(runs)]])$chosen
      wrong <- table(factor(wrong[wrong != t], levels = setdiff(candidates, t)))
      paste(names(wrong), wrong, sep = ": ", collapse = ", ")
    },
    `Selected Q, true design` = chain(function(run) mean(sub(run)[[paste0("Q_", k)]]), 1),
    `ARI, true design` = chain(function(run) mean(sub(run)[[paste0("ARI_", k)]])),
    check.names = FALSE
  )
}
tab_design <- do.call(rbind, lapply(unlist(families), design_row))

## selected Q of every candidate design, last run only
q_row <- function(t) {
  sub <- runs[[length(runs)]]$design
  sub <- sub[sub$truth == t & is.na(sub$error), ]
  candidates <- families[[if (t %in% families$node) "node" else "dyad"]]
  data.frame(`Simulated with` = t,
             `Fitted with` = paste(candidates, collapse = " / "),
             `Selected Q` = paste(fmt(colMeans(sub[, paste0("Q_", 1:3)]), 1), collapse = " / "),
             `ARI` = paste(fmt(colMeans(sub[, paste0("ARI_", 1:3)])), collapse = " / "),
             check.names = FALSE)
}
tab_q <- do.call(rbind, lapply(unlist(families), q_row))

## ---------------------------------------------------------------------------------------
errors <- sapply(runs, function(run) sum(!is.na(run$grid$error)) + sum(!is.na(run$design$error)))
md <- c(
  sprintf("Each cell: %s, mean over %d replicates.", paste(labels, collapse = " → "), max(runs[[1]]$grid$rep)),
  sprintf("Runs: %s.", paste(sprintf("%s at %s, ICL %s (%.0f min)", labels, sapply(runs, `[[`, "commit"),
                                     sapply(runs, function(r) if (is.null(r$ICL)) "default" else r$ICL),
                                     sapply(runs, `[[`, "elapsed")), collapse = "; ")),
  if (any(errors > 0)) sprintf("Fits that failed: %s.", paste(labels, errors, sep = " ", collapse = ", ")),
  "",
  "### Experiment 1: block-node and block-dyad, Q = 1 to 6 explored, true Q = 3",
  "",
  knitr::kable(tab_grid, format = "pipe", align = "lllrrrrrrrr", row.names = FALSE),
  "",
  "### Experiment 2: choice of the design by the ICL (assortative network, medium rate)",
  "",
  knitr::kable(tab_design, format = "pipe", align = "lrrlrr", row.names = FALSE),
  "",
  sprintf("Selected Q and ARI of each candidate design, %s:", labels[length(labels)]),
  "",
  knitr::kable(tab_q, format = "pipe", align = "llrr", row.names = FALSE)
)
writeLines(md, file.path(dir, "summary.md"))
cat(md, sep = "\n")
