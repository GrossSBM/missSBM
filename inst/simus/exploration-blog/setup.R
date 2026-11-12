## Instances of the benchmark: the blog network (frenchblog2007), fully observed, and 5
## samplings with block-node as in the JSS analysis (inst/JSS-analyses/main_analysis.R): the
## nodes of the blocks with less than 10% of the blogs observed with probability 0.2, the
## others with 0.8, the blocks being those of a reference fit on the full network.
##
## Usage (from the root of the package):
##   Rscript inst/simus/exploration-blog/setup.R <path to missSBM source>
## writes inst/simus/exploration-blog/results/instances.rds

pkg <- commandArgs(TRUE)[1]
dir <- file.path("inst", "simus", "exploration-blog", "results")
dir.create(dir, showWarnings = FALSE, recursive = TRUE)
suppressMessages({devtools::load_all(pkg, quiet = TRUE); library(igraph)})
future::plan("sequential")

fb   <- delete_vertices(frenchblog2007, which(degree(frenchblog2007) == 0))
blog <- as_adjacency_matrix(fb)
set.seed(1)
ref <- estimateMissSBM(blog, 1:14, "node", control = missSBM_param(trace = FALSE, iterates = 3))$bestModel
cl  <- ref$fittedSBM$memberships
rho <- ifelse(ref$fittedSBM$blockProp < 0.1, 0.2, 0.8)

inst <- list(full = list(A = blog, designs = "node"))
for (s in c(3052008, 1:4)) { # 3052008: the seed of the JSS analysis
  set.seed(s)
  inst[[paste0("block-node_", s)]] <- list(
    A = observeNetwork(adjacencyMatrix = blog, sampling = "block-node", parameters = rho, clusters = cl),
    designs = c("node", "block-node"))
}
saveRDS(list(inst = inst, ref_cl = cl, rho = rho, party = as.numeric(as.factor(vertex.attributes(fb)$party))),
        file.path(dir, "instances.rds"))
cat("reference fit:", ref$fittedSBM$nbBlocks, "blocks; observed dyads:",
    round(sapply(inst, function(i) mean(!is.na(i$A[upper.tri(i$A)]))), 2), "\n")
