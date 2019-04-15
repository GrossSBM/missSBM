rm(list=ls())
library(missSBM)
library(igraph)

data(war)
beligerent_adjacency <- as_adj(war$beligerent, sparse = FALSE)
nWar = nrow(beligerent_adjacency)
matsumpower = matrix(war$beligerent$power,nrow = nWar,ncol = nWar) + matrix(war$beligerent$power,nrow = nWar,ncol = nWar,byrow=T)
diag(matsumpower) = 0 # diagonal has to be set to 0
sampleNet_cov <- prepare_data(beligerent_adjacency, covariates = list(matsumpower))

vBlocks <- 1:5
collection_sbm_cov_full <- estimate(sampleNet_cov, vBlocks = vBlocks, sampling = "dyad")
res_unsmoothed <- data.frame(
  ICL     = collection_sbm_cov_full$ICL,
  nBlocks = vBlocks,
  type    = "raw"
)
smooth(collection_sbm_cov_full, "both", control = list(iterates = 2))
res_smoothed <- data.frame(
  ICL     = collection_sbm_cov_full$ICL,
  nBlocks = vBlocks,
  type    = "smoothed"
)
rbind(res_unsmoothed, res_smoothed) %>%
  ggplot(aes(x = nBlocks, y = ICL, group = type, color = type)) +
  geom_line() + theme_bw()

# this is working well


# because the covariates are a matrix, a node sampling is not accepted
collection_sbm_cov_full <- estimate(sampleNet_cov, vBlocks = vBlocks, sampling = "node")


# other solution

sampleNet_cov2 <- prepare_data(beligerent_adjacency, covariates = list(war$beligerent$power),similarity = function(v1,v2) {
  v1+v2})
sampleNet_cov2$covarArray
collection_sbm_cov_full2 <- estimate(sampleNet_cov2, vBlocks = vBlocks, sampling = "node")
smooth(collection_sbm_cov_full2)
collection_sbm_cov_full2$bestModel$fittedSBM$covarParam

