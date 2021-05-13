library(missSBM)
library(greed)
library(aricode)
library(ggplot2)

data("frenchblog2007", package = "missSBM")
adjacencyMatrix <- igraph::as_adj(frenchblog2007)
party <- igraph::vertex.attributes(frenchblog2007)$party
vBlocks <- 1:20

greed_full <- greed(adjacencyMatrix, model = new("sbm", type = "undirected"))

sbm_full  <- estimateMissSBM(adjacencyMatrix, vBlocks <- 1:20, "node", control = list(core = 10, iterates = 1))
plot(sbm_full, "monitoring")

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

sbm_mar   <- estimateMissSBM(sampledNet, 1:10, "node", control = list(trace = 2, core = 4, iterates = 10))
greed_mar <- greed(sampledNet, model = new("misssbm", type = 'undirected', sampling = 'dyad'))

ARI(sbm_mar$models[[5]]$fittedSBM$memberships, party)
ARI(sbm_mar$models[[5]]$fittedSBM$memberships, sbm_full$bestModel$fittedSBM$memberships)
ARI(greed_mar@cl, party)
ARI(greed_mar@cl, greed_full@cl)
ARI(greed_mar@cl, sbm_full$bestModel$fittedSBM$memberships)
ARI(sbm_full$bestModel$fittedSBM$memberships, party)
