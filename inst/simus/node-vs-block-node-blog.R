## block-node vs node on the blog network: 12 block-node samplings (rates 0.2 for the blocks with
## less than 10% of the blogs, 0.8 for the others, blocks of a fit on the full network), each fitted
## with node and block-node (Q = 1 to 12, 2 passes): ICL, ARI with the full network, AUC of the
## imputations. Usage (root of the package):
##   SIMU_LIB=<lib-v1.2.0> OPENBLAS_NUM_THREADS=1 Rscript inst/simus/node-vs-block-node-blog.R
source(file.path("inst", "simus", "load-missSBM.R"))
library(igraph)
fb <- delete_vertices(frenchblog2007, which(degree(frenchblog2007) == 0))
blog <- as_adjacency_matrix(fb, sparse = FALSE)
set.seed(3052008)
full <- estimateMissSBM(blog, 1:12, "node", control = missSBM_param(trace = FALSE))
cl <- full$bestModel$fittedSBM$memberships
rate <- ifelse(full$bestModel$fittedSBM$blockProp < 0.1, 0.2, 0.8)
auc <- function(truth, prediction) { r <- rank(prediction); n1 <- sum(truth == 1); n0 <- sum(truth == 0); (sum(r[truth == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
res <- do.call(rbind, parallel::mclapply(1:12, function(s) {
  set.seed(s)
  A <- observeNetwork(blog, "block-node", rate, clusters = cl)
  ctl <- missSBM_param(trace = FALSE, iterates = 2)
  n <- estimateMissSBM(A, 1:12, "node", control = ctl); b <- estimateMissSBM(A, 1:12, "block-node", control = ctl)
  miss <- which(is.na(A) & upper.tri(A))
  data.frame(seed = s, dICL = min(n$ICL) - min(b$ICL),
    Q_node = n$bestModel$fittedSBM$nbBlocks, Q_block = b$bestModel$fittedSBM$nbBlocks,
    ARI_node = aricode::ARI(cl, n$bestModel$fittedSBM$memberships), ARI_block = aricode::ARI(cl, b$bestModel$fittedSBM$memberships),
    AUC_node = auc(blog[miss], predict(n$bestModel)[miss]), AUC_block = auc(blog[miss], predict(b$bestModel)[miss]))
}, mc.cores = 12))
print(round(res, 3)); print(round(colMeans(res[,-1]), 3))
