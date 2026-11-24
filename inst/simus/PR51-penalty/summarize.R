## Tables of inst/simus/PR51-penalty/README.md, written to results/summary.md
res <- readRDS("inst/simus/PR51-penalty/results/results.rds")
out <- "inst/simus/PR51-penalty/results/summary.md"
agg <- aggregate(cbind(observed_dyads, right_all = Q_all == Qtrue, right_observed = Q_observed == Qtrue,
                       Qhat_all = Q_all, Qhat_observed = Q_observed) ~ net + Qtrue + design + rate + n,
                 res, mean)
agg <- agg[order(agg$net, agg$design, agg$rate, agg$n), ]
fmt <- data.frame(
  Network = agg$net, `True Q` = agg$Qtrue, Design = agg$design, Rate = agg$rate, n = agg$n,
  `Observed dyads` = sprintf("%.2f", agg$observed_dyads),
  `True Q selected` = sprintf("%.2f → %.2f", agg$right_all, agg$right_observed),
  `Mean selected Q` = sprintf("%.1f → %.1f", agg$Qhat_all, agg$Qhat_observed),
  check.names = FALSE)
diff <- res[res$Q_all != res$Q_observed, c("net", "n", "design", "rate", "Qtrue", "Q_all", "Q_observed")]
sink(out)
cat("Each cell: all the dyads → the observed dyads in the penalty, mean over 20 replicates.\n")
cat("Commit:", attr(res, "commit"), "\n\n")
print(knitr::kable(fmt, align = "r"))
cat("\n\nNetworks whose selected Q differs between the two penalties:", nrow(diff), "out of", nrow(res), "\n\n")
print(knitr::kable(diff, row.names = FALSE))
sink()
