rm(list=ls())
library(missSBM)
library(igraph)

data(war)
beligerent_adjacency <- as_adj(war$beligerent, sparse = FALSE)


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





