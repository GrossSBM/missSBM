library(missSBM)
library(aricode)
library(ggplot2)

data("frenchblog2007", package = "missSBM")
class(frenchblog2007)
adjacencyMatrix <- frenchblog2007 %>% igraph::as_adj(sparse = FALSE)
party <- igraph::vertex.attributes(frenchblog2007)$party
vBlocks <- 1:12

sbm_full  <- estimateMissSBM(adjacencyMatrix, vBlocks, "node", control = list(core = 10))

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

sbm_node  <- estimateMissSBM(sampledNet, vBlocks, "node", control = list(trace = 2, core = 4))

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
  ,
  data.frame(Q = vBlocks, ICL = sbm_block$ICL , sampling = "block-node"),
  data.frame(Q = vBlocks, ICL = sbm_full$ICL , sampling = "fully observed")
)

ggplot(ICLs, aes(x = Q, y = ICL, color = sampling)) + theme_bw(base_size = 20) +
  theme(axis.title = element_blank()) + geom_point() + geom_line()

ARI(sbm_full$bestModel$fittedSBM$memberships, sbm_node$bestModel$fittedSBM$memberships)

ARI(sbm_full$bestModel$fittedSBM$memberships, sbm_block$bestModel$fittedSBM$memberships)
