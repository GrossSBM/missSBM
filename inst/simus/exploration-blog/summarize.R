## Tables of the benchmark, written to inst/simus/exploration-blog/results/summary.md
##
## Usage (from the root of the package):
##   Rscript inst/simus/exploration-blog/summarize.R [labels of the runs, in order]

args   <- commandArgs(TRUE)
labels <- if (length(args)) args else c("master", "PR48-completed", "PR48-observed")
dir    <- file.path("inst", "simus", "exploration-blog", "results")
I      <- readRDS(file.path(dir, "instances.rds"))
fmt    <- function(x, d = 1) formatC(x, format = "f", digits = d)

n_local_min <- function(x) sum(sapply(seq_along(x), function(i) all(x[i] < x[intersect(c(i - 1, i + 1), seq_along(x))])))
R <- do.call(rbind, lapply(labels, function(l) {
  b <- readRDS(file.path(dir, paste0(l, ".rds")))
  df  <- b$jobs
  ICL <- t(sapply(b$res, `[[`, "ICL"))
  ## best known ICL for each number of blocks: the best over the settings, for the same version,
  ## instance and design (in practice, the most explored setting)
  key  <- paste(df$instance, df$design)
  best <- do.call(rbind, lapply(key, function(k) apply(ICL[key == k, , drop = FALSE], 2, min)))
  gap  <- ICL - best
  df$version    <- l
  df$time       <- sapply(b$res, `[[`, "time")
  df$gap_mean   <- rowMeans(gap)
  df$gap_max    <- apply(gap, 1, max)
  df$gap_sel    <- apply(ICL, 1, min) - sapply(key, function(k) min(ICL[key == k, ]))
  df$local_min  <- apply(ICL, 1, n_local_min)
  df$Q_sel      <- apply(ICL, 1, which.min)
  df$degenerate <- sapply(b$res, function(r) sum(r$degenerate))
  df$ARI_ref    <- sapply(b$res, function(r) aricode::ARI(r$best, I$ref_cl))
  df
}))
R$setting <- factor(R$setting, c("none", "polish", "default", "iter2", "iter3"))

table_for <- function(sub) {
  sub <- sub[order(match(sub$version, labels), sub$setting), ]
  do.call(rbind, lapply(split(sub, list(sub$setting, factor(sub$version, labels)), drop = TRUE), function(d) data.frame(
    Version = d$version[1], Setting = as.character(d$setting[1]),
    `Time (s)` = fmt(mean(d$time)),
    `ICL gap, mean over Q` = fmt(mean(d$gap_mean)),
    `ICL gap, max` = fmt(max(d$gap_max)),
    `ICL gap, selected model` = fmt(mean(d$gap_sel)),
    `Local minima in Q` = fmt(mean(d$local_min)),
    `Collapsed models` = sum(d$degenerate),
    `Selected Q` = paste(d$Q_sel, collapse = ", "),
    `ARI with the reference` = fmt(mean(d$ARI_ref), 2),
    check.names = FALSE)))
}
groups <- list(
  "Fully observed network, `node`" = R$instance == "full",
  "5 block-node samplings, fitted with `node`" = R$instance != "full" & R$design == "node",
  "5 block-node samplings, fitted with `block-node`" = R$instance != "full" & R$design == "block-node")
md <- c(
  sprintf("Runs: %s. Means over the instances of each group; the ICL gap is measured, for each number of blocks, against the best ICL obtained by the same version over all the settings.", paste(labels, collapse = ", ")),
  "The times are those of runs with 18 jobs in parallel on 20 cores, except `iter2`, run afterwards with 11 jobs: compare them within a run.",
  unlist(lapply(names(groups), function(g) c("", paste("###", g), "",
    knitr::kable(table_for(R[groups[[g]], ]), format = "pipe", row.names = FALSE, align = "llrrrrrrlr")))))
writeLines(md, file.path(dir, "summary.md"))
cat(md, sep = "\n")
