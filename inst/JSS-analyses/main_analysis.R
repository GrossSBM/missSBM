## required libraries
library("missSBM")
library("aricode")
library("tidyverse"); theme_set(theme_bw())
library("igraph")
library("pROC")
library("parallel")
library("future")

options(future.fork.enable = TRUE) # so it can be used from Rstudio
future::plan("multicore", workers = 10)

## import blog data
data("frenchblog2007", package = "missSBM")
frenchblog2007 <- delete_vertices(frenchblog2007, which(degree(frenchblog2007) == 0))
blog  <- as_adjacency_matrix(frenchblog2007)
party <- vertex.attributes(frenchblog2007)$party
party_num <- as.numeric(as.factor(party))

## SBM with no missing edge or node
blocks <- 1:18
sbm_full <- estimateMissSBM(blog, blocks, "node")
plot(sbm_full)
plot(sbm_full, "monitoring")
plot(sbm_full$bestModel)
coef(sbm_full$bestModel, type = "mixture")
plot(sbm_full$bestModel, dimLabels = list(row = "blogs", col = "blogs"))

set.seed(03052008)

## sampling in the original network to get partly a observed blog network
samplingParameters <- ifelse(sbm_full$bestModel$fittedSBM$blockProp <  0.1, 0.2, 0.8)
blog_obs <- observeNetwork(adjacencyMatrix = blog, sampling = "block-node",
  parameters = samplingParameters,
  clusters = sbm_full$bestModel$fittedSBM$memberships)

## try MAR
sbm_node <- estimateMissSBM(blog_obs, blocks, "node", control = list(iterates = 3))
plot(sbm_node)
plot(sbm_node, type = "monitoring")

## try MNAR with block-dyad sampling
sbm_block <- estimateMissSBM(blog_obs, blocks, "block-node", control = list(iterates = 5))
plot(sbm_block)
plot(sbm_block, type = "monitoring")

## compare the models via ICL
rbind(tibble(Q = blocks, ICL = sbm_node$ICL, sampling  = "node"),
  tibble(Q = blocks, ICL = sbm_block$ICL, sampling = "block-node"),
  tibble(Q = blocks, ICL = sbm_full$ICL, sampling = "fully observed")
)%>% ggplot(aes(x = Q, y = ICL, color = sampling)) +
  geom_line() + geom_point() + ggtitle("Model Selection") +
  labs(x = "# blocks", y = "Integrated Classification Likelihood")

## compare the models via ARI
ARI(sbm_full$bestModel$fittedSBM$memberships, sbm_node$bestModel$fittedSBM$memberships)
ARI(sbm_full$bestModel$fittedSBM$memberships, sbm_block$bestModel$fittedSBM$memberships)

myModel <- sbm_block$bestModel
ARI(party_num, myModel$fittedSBM$memberships)
ARI(party_num, sbm_full$bestModel$fittedSBM$memberships)
plot(myModel, type = "expected", dimLabels = list(row = "blogs", col = "blogs"))

## Alluvial plot
library("alluvial")
cl_party <- as.factor(party)
levels(cl_party)[3] <- "center-right"
new_order <- match(c("far-right", "right", "liberal", "center-right", "analyst",
                     "center-left", "left", "green", "far-left"), levels(cl_party))
cl_party <- factor(cl_party, levels = levels(cl_party)[new_order], ordered = TRUE)
cl_block <- as.factor(myModel$fittedSBM$memberships)
cl_block <- factor(cl_block, levels = levels(cl_block)[c(1, 6, 3, 4, 10, 5, 8, 2, 7, 9)],
                   ordered = TRUE)
d <- data.frame(party = cl_party, block = cl_block, Freq = 1)
d %>% group_by(party, block) %>% summarise(n = sum(Freq)) -> alluv
palettes <- rainbow(14, alpha=.35)
cols <- palettes[rev(as.numeric(alluv$party))]
alluvial(alluv[, c(1,2)], freq = alluv$n, col = cols, border = cols, blocks = FALSE, alpha = 0.35)

## Imputation and AUC
cl0 <- sbm_full$bestModel$fittedSBM$memberships
nBlocks <- sbm_full$bestModel$fittedSBM$nbBlocks

future::plan("sequential")

res_auc <- mclapply(1:500, function(i) {
  subGraph <- observeNetwork(blog, "block-node", runif(nBlocks), cl0)
  missing <- which(as.matrix(is.na(subGraph)))
  true_dyads <- blog[missing]
  sbm_block <- estimateMissSBM(subGraph, nBlocks, "block-node", control = list(cores = 1, trace = 0))
  imputed_dyads <- sbm_block$bestModel$imputedNetwork[missing]
  c(rate = 1 - length(missing)/length(blog), auc  = auc(true_dyads, imputed_dyads, quiet = TRUE))
}, mc.cores = 10)
# purrr::reduce(res_auc, rbind) %>% as.data.frame() %>%
#   ggplot() + aes(x = rate, y = auc) + geom_point(size = 0.25) +
#   geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
#   labs(x = "sampling rate", y = "Area under the ROC curve")

## consider a subgraph and a nodal covariate
future::plan("multicore", workers = 10)
blog_subgraph <- frenchblog2007 %>%
  igraph::induced_subgraph(V(frenchblog2007)$party %in%
  c( "right", "left"))
blog_subgraph <- delete_vertices(blog_subgraph,
  which(degree(blog_subgraph) ==  0))
dummy_party <- (V(blog_subgraph)$party == "left") * 1

## observe this graph according to this covariate
blog_subgraph_obs <- blog_subgraph %>% as_adjacency_matrix() %>%
  missSBM::observeNetwork(sampling="covar-node", parameters = 3,
    covariates = list(dummy_party))

## fit SBM on the fully observed, taking the covariate into account in the SBM
blocks <- 2:8
sbm_covar_full <- blog_subgraph %>% as_adjacency_matrix() %>%
  estimateMissSBM(blocks, "node", covariates =  list(dummy_party),
    control = list(useCov = TRUE, iterates = 2))
# plot(sbm_covar_full, "icl")
# plot(sbm_covar_full, "monitoring")

## fit SBM on the partly observed, taking the covariate into account in the SBM and the sampling
sbm_covar1 <- estimateMissSBM(blog_subgraph_obs, blocks, "covar-node",
  covariates =  list(dummy_party), control = list(useCov = TRUE, iterates = 2))
# plot(sbm_covar1, "icl")
# plot(sbm_covar1, "monitoring")

## fit SBM on the partly observed, taking the covariate into account only in the sampling
sbm_covar2 <- estimateMissSBM(blog_subgraph_obs, blocks, "covar-node",
  covariates = list(dummy_party), control = list(useCov = FALSE, iterates = 2))
# plot(sbm_covar2, "icl")
# plot(sbm_covar2, "monitoring")

## fit SBM on the partly observed, taking the covariate into account only in the SBM
sbm_covar3 <- estimateMissSBM(blog_subgraph_obs, blocks, "node",
  covariates =  list(dummy_party), control = list(useCov = TRUE, iterates = 2))
# plot(sbm_covar3,  "icl")
# plot(sbm_covar3, "monitoring")

## fit SBM on the partly observed, not taking the covariate into account
sbm_covar4 <- estimateMissSBM(blog_subgraph_obs, blocks, "node",
  control = list(useCov = FALSE, iterates = 2))
# plot(sbm_covar4, "icl")
# plot(sbm_covar4, "monitoring")

sbm_covar_full <- as_adjacency_matrix(blog_subgraph) %>%
   estimateMissSBM(blocks, "node", covariates =  list(dummy_party))


### compare all the SBM models in the covariate framework
# rbind(tibble(Q = blocks, ICL = sbm_covar1$ICL, sampling = "covar-node", useCov = "true"),
#       tibble(Q = blocks, ICL = sbm_covar2$ICL, sampling = "covar-node", useCov = "false"),
#       tibble(Q = blocks, ICL = sbm_covar3$ICL, sampling = "node", useCov = "true"),
#       tibble(Q = blocks, ICL = sbm_covar4$ICL, sampling = "node", useCov = "false")
# ) %>% ggplot(aes(x = Q, y = ICL, color = sampling, shape = useCov))  +
#   geom_line() + geom_point() + labs(x = "#blocks", y = "ICLs")


# Blocks <- sbm_covar_full$bestModel$fittedSBM$memberships
# Blocks1 <- sbm_covar1$bestModel$fittedSBM$memberships
# Blocks2 <- sbm_covar2$bestModel$fittedSBM$memberships
# Blocks3 <- sbm_covar3$bestModel$fittedSBM$memberships
# Blocks4 <- sbm_covar4$bestModel$fittedSBM$memberships
# ARIs <- matrix(c(aricode::ARI(Blocks, Blocks1), aricode::ARI(Blocks, Blocks2),
#   aricode::ARI(Blocks, Blocks3), aricode::ARI(Blocks, Blocks4),
#   NA, aricode::ARI(Blocks1, Blocks2), aricode::ARI(Blocks1, Blocks3), aricode::ARI(Blocks1, Blocks4),
#   NA, NA, aricode::ARI(Blocks2, Blocks3), aricode::ARI(Blocks2, Blocks4),
#   NA, NA, NA, aricode::ARI(Blocks3, Blocks4)), nrow = 4, ncol = 4, byrow = TRUE)
# ARIs <- as.data.frame(ARIs)
# rownames(ARIs) <- c("Full", "i", "ii", "iii")
# colnames(ARIs) <- c("i", "ii", "iii", "iv")
#
# ARIs

future::plan("sequential")

## Save output for reproducibility in the Rnw
save.image(file = "main_analysis.RData")
