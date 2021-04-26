library(missSBM)
library(igraph)
library(ggplot2)
library(future.apply)

data("frenchblog2007", package = "missSBM")
class(frenchblog2007)
adjacencyMatrix <- frenchblog2007 %>% as_adj(sparse = FALSE)
party <- vertex.attributes(frenchblog2007)$party


vBlocks <- 1:14
control <- list(trace = 0, iterates = 2, cores = 1)
smoothing_type <- "both"

sbm_full <- estimateMissSBM(adjacencyMatrix, vBlocks, "node",  control = control)
smooth(sbm_full, smoothing_type, control)


psi <- 0.8
sampledNet <-
  observeNetwork(
    adjacencyMatrix = adjacencyMatrix,
    sampling        = "block-node",
    parameters      = psi,
  )



sbm_node <- estimateMissSBM(sampledNet, vBlocks, "node", control = control)
smooth(sbm_node, smoothing_type, control)

ICLs <- rbind.data.frame(
  data.frame(Q = vBlocks, ICL = sbm_node$ICL , sampling = "node"),
  data.frame(Q = vBlocks, ICL = sbm_node$ICL , sampling = "block-node"),
  data.frame(Q = vBlocks, ICL = sbm_full$ICL , sampling = "fully observed")
)

p <- ggplot(ICLs, aes(x = Q, y = ICL, color = sampling)) + theme_bw(base_size = 20) +
  theme(axis.title = element_blank()) + geom_point() + geom_line()

print(p)
