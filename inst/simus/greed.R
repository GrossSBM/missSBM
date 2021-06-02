library(sbm)
library(missSBM)
library(greed)
library(aricode)
library(ggplot2)

data("frenchblog2007", package = "missSBM")
adjacencyMatrix <- igraph::as_adj(frenchblog2007)
party <- igraph::vertex.attributes(frenchblog2007)$party
vBlocks <- 1:20

options(future.fork.enable = TRUE)
future::plan("multicore", workers = 10)

greed_full <- greed(adjacencyMatrix, model = new("sbm", type = "undirected"))

sbm_full  <- estimateMissSBM(adjacencyMatrix, vBlocks <- 1:20, "node")

plot(sbm_full, "monitoring")

sampledNet <-
  observeNetwork(
    adjacencyMatrix = adjacencyMatrix,
    sampling        = "block-dyad",
    parameters      = sbm_full$bestModel$fittedSBM$connectParam$mean,
    clusters        = sbm_full$bestModel$fittedSBM$memberships
  )

R <- as.matrix((!is.na(sampledNet)) * 1)
diag(R) <- NA

apply(cutree(hclust(dist(R, "manhattan")), 1:20), 2, ARI, sbm_full$bestModel$fittedSBM$memberships)

ARI(sbm_mar$models[[5]]$fittedSBM$memberships, party)

sbm_R   <- estimateMissSBM(R, 1:20, "node")

sbm_mar   <- estimateMissSBM(sampledNet, 1:10, "node", control = list(trace = 2, core = 4, iterates = 10))
greed_mar <- greed(sampledNet, model = new("misssbm", type = 'undirected', sampling = 'dyad'))

ARI(sbm_mar$models[[5]]$fittedSBM$memberships, party)
ARI(sbm_mar$models[[5]]$fittedSBM$memberships, sbm_full$bestModel$fittedSBM$memberships)
ARI(greed_mar@cl, party)
ARI(greed_mar@cl, greed_full@cl)
ARI(greed_mar@cl, sbm_full$bestModel$fittedSBM$memberships)
ARI(sbm_full$bestModel$fittedSBM$memberships, party)
