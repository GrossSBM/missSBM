library(missSBM)
library(igraph)
library(ggplot2)
library(future.apply)

set.seed(222)
#future::plan("sequential")
future::plan("multisession", workers = 10)

data("frenchblog2007", package = "missSBM")
class(frenchblog2007)
adjacencyMatrix <- frenchblog2007 %>% as_adj(sparse = FALSE)
party <- vertex.attributes(frenchblog2007)$party
vBlocks <- 1:13

sbm_full  <- estimateMissSBM(adjacencyMatrix, vBlocks, "node")

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


sbm_node  <- estimateMissSBM(sampledNet, vBlocks, "node", control = list(trace = 2, iterates = 1))

sbm_block <- estimateMissSBM(sampledNet, vBlocks, "block-node", control = list(trace = 2, iterates = 3))

optimStatus <- rbind(
  data.frame(sbm_full$optimizationStatus, sampling = "full"),
  data.frame(sbm_node$optimizationStatus, sampling = "node"),
  data.frame(sbm_block$optimizationStatus, sampling = "block-node")
)

optimStatus %>% dplyr::group_by(sampling) %>%
ggplot(aes(x = cumsum(iteration), y = elbo, color = sampling)) + theme_bw(base_size = 20) +
  theme(axis.title = element_blank()) + geom_point() + geom_line()

ICLs <- rbind.data.frame(
  data.frame(Q = vBlocks, ICL = sbm_node$ICL , sampling = "node"),
  data.frame(Q = vBlocks, ICL = sbm_block$ICL , sampling = "block-node"),
  data.frame(Q = vBlocks, ICL = sbm_full$ICL , sampling = "fully observed")
)

ggplot(ICLs, aes(x = Q, y = ICL, color = sampling)) + theme_bw(base_size = 20) +
  theme(axis.title = element_blank()) + geom_point() + geom_line()

