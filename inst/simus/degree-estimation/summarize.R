## Tables of inst/simus/degree-estimation/README.md, written to results/summary.md
out <- "inst/simus/degree-estimation/results/summary.md"
files <- sort(Sys.glob("inst/simus/degree-estimation/results/results-n*.rds"))
sink(out)
cat("Mean over 20 replicates per setting. A = stages 1 and 2; B = coordinate ascent from the usual",
    "initialization; C = A, then the coordinate ascent; package = missSBM (implements C). J = the bound J_zeta.",
    "'degree chosen' = the ICL of degree is below those of node and block-node (at the true Q).\n\n")
for (f in files) {
  r <- readRDS(f)
  r$a_err <- r$a_hat - r$a; r$b_rel_err <- (r$b_hat - r$b) / ifelse(r$b == 0, NA, r$b)
  r$degree_chosen <- r$ICL < pmin(r$ICL_node, r$ICL_block_node)
  r$effect <- factor(r$effect, c("none", "weak", "medium", "strong"))
  r$method <- factor(r$method, c("A", "B", "C", "package"))
  agg <- aggregate(cbind(observed_nodes, ARI, ARI_vs_A, rmse_pi, a_err, J, degree_chosen, iterations, time) ~
                   method + effect, r, mean, na.action = na.pass)
  b_rel <- aggregate(b_rel_err ~ method + effect, r, mean, na.action = na.pass)
  agg <- merge(agg, b_rel, all.x = TRUE)
  agg <- agg[order(agg$effect, agg$method), ]
  ## the gain of the bound with respect to A, network by network
  w <- reshape(r[, c("effect", "rep", "method", "J", "ARI", "ICL")], idvar = c("effect", "rep"),
               timevar = "method", direction = "wide")
  fmt <- data.frame(
    Effect = agg$effect, Method = agg$method,
    `Observed nodes` = sprintf("%.2f", agg$observed_nodes),
    ARI = sprintf("%.3f", agg$ARI), `ARI vs A` = sprintf("%.3f", agg$ARI_vs_A),
    `RMSE pi` = sprintf("%.4f", agg$rmse_pi),
    `a error` = sprintf("%+.2f", agg$a_err),
    `b rel. error` = ifelse(is.na(agg$b_rel_err), "", sprintf("%+.1f%%", 100 * agg$b_rel_err)),
    `J - J(A)` = sprintf("%+.2f", agg$J - ave(ifelse(agg$method == "A", agg$J, NA), agg$effect,
                                             FUN = function(x) x[!is.na(x)])),
    `Degree chosen` = sprintf("%.2f", agg$degree_chosen),
    Iterations = sprintf("%.1f", agg$iterations), `Time (s)` = sprintf("%.2f", agg$time),
    check.names = FALSE)
  cat(sprintf("### n = %d (%s)\n\n", r$n[1], attr(r, "version")))
  print(knitr::kable(fmt, align = "r"))
  cat(sprintf("\n\nJ(C) - J(A): min %.3f, mean %.2f; J(C) - J(B): mean %.2f, J(C) < J(B) in %d networks out of %d.",
              min(w$J.C - w$J.A), mean(w$J.C - w$J.A), mean(w$J.C - w$J.B), sum(w$J.C < w$J.B - 1e-6), nrow(w)))
  cat(sprintf(" ARI of C below that of A by more than 0.01 in %d networks, above in %d.\n\n",
              sum(w$ARI.C < w$ARI.A - .01), sum(w$ARI.C > w$ARI.A + .01)))
  cat(sprintf("Package vs C: ICL lower by more than 1 in %d networks, higher in %d (out of %d); ARI mean difference %+.3f.\n\n",
              sum(w$ICL.package < w$ICL.C - 1), sum(w$ICL.package > w$ICL.C + 1), nrow(w), mean(w$ARI.package - w$ARI.C)))
}
sink()
