library(missSBM)
library(igraph)
library(ggplot2)

set.seed(222)

data("frenchblog2007", package = "missSBM")
adjacencyMatrix <- frenchblog2007 %>%  as_adj(sparse = FALSE)
party <- vertex.attributes(frenchblog2007)$party

vBlocks <- 1:12
control <- list(trace = 1, cores = 10)
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
  data.frame(Q = vBlocks, ICL = sbm_block$ICL, sampling = "block-node"),
  data.frame(Q = vBlocks, ICL = sbm_full$ICL , sampling = "fully observed")
)

p <- ggplot(ICLs, aes(x = Q, y = ICL, color = sampling)) + theme_bw(base_size = 20) +
  theme(axis.title = element_blank()) + geom_point() + geom_line()

print(p)

aricode::NMI(sbm_block$bestModel$fittedSBM$memberships,
             sbm_full$bestModel$fittedSBM$memberships)
aricode::NMI(sbm_node$bestModel$fittedSBM$memberships ,
             sbm_full$bestModel$fittedSBM$memberships)

blog_subgraph <-
  frenchblog2007 %>%
  igraph::induced_subgraph(V(frenchblog2007)$party %in% c( "right", "left"))
blog_subgraph <-
  delete_vertices(blog_subgraph, which(degree(blog_subgraph) == 0))

dummy_party <- dummy_party <- (V(blog_subgraph)$party == "left") * 1
vBlocks <- 1:8

sbm_sub_full_cov <- estimateMissSBM(blog_subgraph %>% as_adj(sparse = FALSE), vBlocks, "node", covariates = list(dummy_party), control = control)

blog_subgraph_obs <-
  blog_subgraph %>% as_adj(sparse = FALSE) %>%
  missSBM::observeNetwork(sampling="covar-node", parameters = 3,
                          covariates = list(dummy_party), intercept =-.5)

control_cov <-  control
control_cov$useCov <- FALSE
sbm_covar2 <- blog_subgraph_obs %>%
  estimateMissSBM(vBlocks, "covar-node", covariates = list(dummy_party), control = control_cov)
sbm_covar2$bestModel$fittedSampling$parameters # sampling parameters

