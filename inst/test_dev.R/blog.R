library(missSBM)
library(igraph)
library(ggplot2)

data("frenchblog2007", package = "missSBM")
adjacencyMatrix <- frenchblog2007 %>%  as_adj(sparse = FALSE)
party <- vertex.attributes(frenchblog2007)$party

vBlocks <- 1:14
control <- list(trace = 1, iterates = 3, prop_swap = c(0, 0.25, 0.5), cores = 10)

sbm_full <- estimateMissSBM(adjacencyMatrix, vBlocks, "node", control = control)


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







blog_subgraph <-
  frenchblog2007 %>%
  igraph::induced_subgraph(V(frenchblog2007)$party %in% c( "right", "left"))
blog_subgraph <-
  delete_vertices(blog_subgraph, which(degree(blog_subgraph) == 0))

dummy_party <- dummy_party <- (V(blog_subgraph)$party == "left") * 1

control <- list(trace = 1, iterates = 2, prop_swap = 0, cores = 10)
sbm_sub_full_cov <- estimateMissSBM(blog_subgraph %>% as_adj(sparse = FALSE), 1:8, "node", covariates = list(dummy_party), control = control)

