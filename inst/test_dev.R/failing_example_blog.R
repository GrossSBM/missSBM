library(missSBM)
library(aricode)
library(ggplot2)
library(igraph)
set.seed(222)

data("frenchblog2007", package = "missSBM")

### with covariates
blog_subgraph <-
  frenchblog2007 %>%
  igraph::induced_subgraph(V(frenchblog2007)$party %in% c( "right", "left"))
blog_subgraph <-
  delete_vertices(blog_subgraph, which(degree(blog_subgraph) == 0))
dummy_party <- (V(blog_subgraph)$party == "left") * 1

blog_subgraph_obs <-
  blog_subgraph %>% as_adj(sparse = FALSE) %>%
  missSBM::observeNetwork(sampling="covar-node", parameters = 3,
                          covariates = list(dummy_party))
vBlocks <- 1:8

## reference
sbm_covar_full <- blog_subgraph %>% as_adj(sparse = FALSE) %>%
  estimateMissSBM(vBlocks, "node", covariates =  list(dummy_party), control = list(cores = 10))
plot(sbm_covar_full, "icl")
plot(sbm_covar_full, "monitoring")

## not working

net <- missSBM:::partlyObservedNetwork$new(blog_subgraph_obs, covariates = list(dummy_party), similarity = missSBM:::l1_similarity)
cl0 <- net$clustering(vBlocks)
purrr::map_dbl(cl0, ARI, sbm_covar_full$bestModel$fittedSBM$memberships)


ctrl <- list(
  threshold = 1e-2, trace = TRUE, cores = 1, imputation = "median", similarity = missSBM:::l1_similarity, useCov = TRUE,
  maxIter = 50, fixPointIter = 3, iterates = 0, prop_swap = 0, smoothing = "both", clusterInit = NULL
)

my_collection <- missSBM_collection$new(net, sampling = "node", cl0, control = ctrl)
## my_collection$estimate(ctrl)

out <- my_collection$models[[5]]$doVEM(ctrl)

sbm_covar1 <-
  estimateMissSBM(blog_subgraph_obs, vBlocks, "covar-node",
                  covariates =  list(dummy_party),
                  control = list(useCov = TRUE, cores = 10))
plot(sbm_covar1, "icl")
plot(sbm_covar1, "monitoring")

sbm_covar3 <-
  estimateMissSBM(blog_subgraph_obs, vBlocks, "node",
                  covariates =  list(dummy_party),
                  control = list(useCov = TRUE, cores = 10))
plot(sbm_covar3,  "icl")
plot(sbm_covar3, "monitoring")

## Working

sbm_covar2 <-
  estimateMissSBM(blog_subgraph_obs, vBlocks, "covar-node",
                  covariates = list(dummy_party),
                  control = list(useCov = FALSE, cores = 10))
plot(sbm_covar2, "icl")
plot(sbm_covar2, "monitoring")


sbm_covar4 <-
  estimateMissSBM(blog_subgraph_obs, vBlocks, "node",
                  control = list(cores = 10))
plot(sbm_covar4,  "icl")
plot(sbm_covar4, "monitoring")

Blocks  <- sbm_covar_full$bestModel$fittedSBM$memberships
Blocks1 <- sbm_covar1$bestModel$fittedSBM$memberships
Blocks2 <- sbm_covar2$bestModel$fittedSBM$memberships
Blocks3 <- sbm_covar3$bestModel$fittedSBM$memberships
Blocks4 <- sbm_covar4$bestModel$fittedSBM$memberships
ARIs <- matrix(c(ARI(Blocks,Blocks1),ARI(Blocks,Blocks2),
  ARI(Blocks,Blocks3),ARI(Blocks,Blocks4),
  NA,ARI(Blocks1,Blocks2),ARI(Blocks1,Blocks3),
  ARI(Blocks1,Blocks4),
  NA,NA,ARI(Blocks2,Blocks3),ARI(Blocks2,Blocks4),
  NA,NA,NA,ARI(Blocks3,Blocks4)),nrow=4,ncol=4,byrow=TRUE)
ARIs <- as.data.frame(ARIs)
rownames(ARIs) <- c("Full", "i", "ii", "iii")
colnames(ARIs) <- c("i", "ii", "iii", "iv")

save.image(file = "last_run.RData")

