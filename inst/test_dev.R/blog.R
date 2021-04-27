library(missSBM)
library(igraph)
library(ggplot2)

data("frenchblog2007", package = "missSBM")
class(frenchblog2007)
adjacencyMatrix <- frenchblog2007 %>% as_adj(sparse = FALSE)
party <- vertex.attributes(frenchblog2007)$party

vBlocks <- 1:14
control <- list(trace = 1, iterates = 3, cores = 10)

sbm_full <- estimateMissSBM(adjacencyMatrix, vBlocks, "node",  control = control)

samplingParameters <- base::sample(
  x       = c(0.2, 0.8),
  size    = sbm_full$bestModel$fittedSBM$nbBlocks,
  replace = TRUE)

sampledNet <-
  observeNetwork(
    adjacencyMatrix = adjacencyMatrix,
    sampling        = "block-node",
    parameters      = samplingParameters,
    clusters        = sbm_full$bestModel$fittedSBM$memberships
  )

sbm_node <- estimateMissSBM(sampledNet, vBlocks, "node", control = control)

sbm_block <- estimateMissSBM(sampledNet, vBlocks, "block-node", control = control)

ICLs <- rbind.data.frame(
  data.frame(Q = vBlocks, ICL = sbm_node$ICL , sampling = "node"),
  data.frame(Q = vBlocks, ICL = sbm_block$ICL , sampling = "block-node"),
  data.frame(Q = vBlocks, ICL = sbm_full$ICL , sampling = "fully observed")
)

p <- ggplot(ICLs, aes(x = Q, y = ICL, color = sampling)) + theme_bw(base_size = 20) +
  theme(axis.title = element_blank()) + geom_point() + geom_line()

print(p)

aricode::ARI(sbm_block$bestModel$fittedSBM$memberships,
             sbm_full$bestModel$fittedSBM$memberships)
aricode::ARI(sbm_node$bestModel$fittedSBM$memberships ,
             sbm_full$bestModel$fittedSBM$memberships)

