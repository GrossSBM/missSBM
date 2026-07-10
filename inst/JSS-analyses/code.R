## ----load-libraries, results  = "hide", message = FALSE, warning = FALSE, echo = TRUE----
library("missSBM")
library("aricode")
library("tidyverse"); theme_set(theme_bw())
library("igraph")
library("future")
library("ggplot2")
library("forcats")

set.seed(03052008)

## ----timings, echo = TRUE, fig.cap = "\\label{fig:timings}Timings for adjusting binary SBM with \\textbf{missSBM} on the French political blogosphere with a single core for a varying number of blocks (50 replicated runs per \\# blocks).", fig.pos = "htbp!", fig.height = 3.5, fig.width = 8, out.width = ".7\\textwidth", fig.align = "center", warning = FALSE, message = FALSE----
# source("timings.R") # replicate interim results saved in timings.RData
load("timings.RData")

ggplot(res_timings) + 
  geom_jitter(aes(x = `#group`, y = `time (sec.)`, group = factor(`#group`)), alpha = .75, size = 1, shape = 16) +
  geom_smooth(method = "lm", aes(x = as.numeric(`#group`), y = `time (sec.)`), se = FALSE) +
  facet_grid(. ~ fct_relevel(missingness, c("none", " MAR", "NMAR"))) +
  labs(x = "# blocks") + ylim(c(0, 2)) + theme_bw()

## ----estimate-function, eval = FALSE------------------------------------------
# estimateMissSBM(adjacencyMatrix, vBlocks, sampling = "dyad",
#   covariates = list(), control = list())


## ----available-samplings------------------------------------------------------
missSBM:::available_samplings


## ----sample-function, eval = FALSE--------------------------------------------
## observeNetwork(adjacencyMatrix, sampling, parameters, clusters = NULL, 
##   covariates = list(), similarity = l1_similarity, intercept = 0)


## ----set-future-beginning, eval = FALSE---------------------------------------
## future::plan("multicore", workers = 10)


## ----load results, echo = TRUE, results  = "hide", message = FALSE, warning = FALSE----
# source("main_analysis.R") # replicate interim results saved in main_analysis.RData
load("main_analysis.RData")


## ----load-data, message = FALSE-----------------------------------------------
data("frenchblog2007", package = "missSBM")
frenchblog2007 <- delete_vertices(frenchblog2007, 
  which(degree(frenchblog2007) ==  0))
blog  <- as_adj(frenchblog2007)
party <- vertex.attributes(frenchblog2007)$party


## ----set-blocks1--------------------------------------------------------------
blocks <- 1:18


## ----inference-full-graph, message = FALSE, eval = FALSE----------------------
## sbm_full <- estimateMissSBM(blog, blocks, "node")


## ----plot-monitoring-full, message = FALSE, warning = FALSE, fig.align = "center", fig.height = 5, fig.width = 6, out.width = ".6\\textwidth", fig.pos = "htbp!", fig.cap = "\\label{fig:monitoring}Evolution of the ELBO during the optimization for the successive numbers of blocks"----
plot(sbm_full, type = "monitoring")


## ----print-icl-full-----------------------------------------------------------
which.min(sbm_full$ICL)


## ----print-inference-full-----------------------------------------------------
sbm_full$bestModel


## ----print-coef-full----------------------------------------------------------
coef(sbm_full$bestModel, type = "mixture")


## ----plot-monitoring-best, message = FALSE, warning = FALSE, fig.pos = "htbp!", fig.align = "center", fig.height = 3, fig.width = 4, out.width = ".6\\textwidth", fig.cap = "\\label{fig:imputed}Network data reorganized by the estimated block memberships"----
plot(sbm_full$bestModel, dimLabels = list(row = "blogs", col = "blogs"))


## ----sampling-data, warning = FALSE, eval = FALSE-----------------------------
## samplingParameters <- ifelse(sbm_full$bestModel$fittedSBM$blockProp <  0.1, 0.2, 0.8)
## blog_obs <- observeNetwork(adjacencyMatrix = blog, sampling = "block-node", 
##   parameters = samplingParameters, 
##   clusters = sbm_full$bestModel$fittedSBM$memberships)


## ----inference-block-and-node, message = FALSE, eval = FALSE------------------
## sbm_block <- estimateMissSBM(blog_obs, blocks, "block-node", 
##   control = list(iterates = 5))
## sbm_node <- estimateMissSBM(blog_obs, blocks, "node", 
##   control = list(iterates = 5))


## ----plot-ICLs-blog, echo = TRUE, results  = "hide", warning = FALSE, fig.cap = "\\label{fig:ICL_samplings}ICL for models with fully observed data, block-node sampling and node sampling", fig.align = "center", fig.width = 6, fig.height = 4, out.width = ".6\\textwidth"----
rbind(tibble(Q = blocks, ICL = sbm_node$ICL, sampling  = "node"), 
  tibble(Q = blocks, ICL = sbm_block$ICL, sampling = "block-node"), 
  tibble(Q = blocks, ICL = sbm_full$ICL, sampling = "fully observed")
  )%>% ggplot(aes(x = Q, y = ICL, color = sampling)) + 
  geom_line() + geom_point() + ggtitle("Model Selection") +
  labs(x = "# blocks", y = "Integrated Classification Likelihood")


## ----ARI-party, comment = NA, warning = FALSE-----------------------------------
ARI(sbm_block$bestModel$fittedSBM$memberships, 
    sbm_full$bestModel$fittedSBM$memberships)
ARI(sbm_node$bestModel$fittedSBM$memberships, 
    sbm_full$bestModel$fittedSBM$memberships)


## ----output-missSBMfit, comment = NA------------------------------------------
myModel <- sbm_block$bestModel


## ----output-SBMfit, comment = NA----------------------------------------------
myModel$fittedSBM


## ----output-samplingfit, comment = NA-----------------------------------------
myModel$fittedSampling


## ----connectivity-parameters, comment = NA, fig.align = "center", fig.width = 4, fig.height = 4, out.width = ".6\\textwidth", fig.cap = "\\label{fig:connectivity}Probabilities of connection predicted by the SBM with block-node sampling", fig.pos = "htbp!", warning = FALSE----
plot(myModel, type = "expected", dimLabels = list(row = "blogs", col = "blogs"))


## ----comparison-clustering, comment = NA, warning = FALSE---------------------
ARI(party, myModel$fittedSBM$memberships)
ARI(party, sbm_full$bestModel$fittedSBM$memberships)


## ----Alluvial-plot, echo = TRUE, message = FALSE, results  = "hide", warning = FALSE, fig.cap = "\\label{fig:alluvialPlot}Alluvial plot between block-node sampling clustering and political parties.", fig.width = 9, fig.height = 7, out.width = ".8\\textwidth", fig.pos = "htbp!", fig.align = "center"----
library("alluvial")
cl_party <- as.factor(party)
levels(cl_party)[3] <- "center-right"
new_order <- match(c("far-right", "right", "liberal", "center-right", "analyst", "center-left", "left", "green", "far-left"), levels(cl_party))
cl_party <- factor(cl_party, levels = levels(cl_party)[new_order], ordered = TRUE)

cl_block <- as.factor(myModel$fittedSBM$memberships)
cl_block <- factor(cl_block, levels = levels(cl_block)[c(1, 6, 3, 4, 10, 5, 8, 2, 7, 9)], ordered = TRUE)

d <- data.frame(party = cl_party, 
                block = cl_block, Freq = 1)
d %>% group_by(party, block) %>% summarise(n = sum(Freq)) -> alluv

palettes <- rainbow(14, alpha = .35)
cols <- palettes[rev(as.numeric(alluv$party))]
alluvial(alluv[, c(1, 2)], freq = alluv$n, col = cols, border = cols, blocks = FALSE, alpha = .35)


## ----AUC-sampling, eval = FALSE-------------------------------------------------
## library("pROC")
## library("parallel")
## cl0     <- sbm_full$bestModel$fittedSBM$memberships
## nBlocks <- sbm_full$bestModel$fittedSBM$nbBlocks
## future::plan("sequential")
## res_auc <- mclapply(1:500, function(i) {
##   subGraph   <- observeNetwork(blog, "block-node", runif(nBlocks), cl0)
##   missing    <- which(as.matrix(is.na(subGraph)))
##   true_dyads <- blog[missing]
##   sbm_block  <- estimateMissSBM(subGraph, nBlocks, "block-node", 
##                     control = list(cores = 1, trace = 0))
##   imputed_dyads <- sbm_block$bestModel$imputedNetwork[missing]
##   c(rate = 1 - length(missing)/length(blog), 
##     auc  = auc(true_dyads, imputed_dyads, quiet = TRUE))
## }, mc.cores = 10)


## ----AUC-sampling-show, fig.cap = "\\label{fig:AUC}Area Under the Curve (AUC) of the imputation as a function of the sampling rate.", fig.pos = "htbp!", fig.height = 4, fig.width = 5, out.width = ".6\\textwidth", fig.align = "center", message = FALSE, echo = TRUE----
purrr::reduce(res_auc, rbind) %>% as.data.frame() %>% 
  ggplot() + aes(x = rate, y = auc) + geom_point(size = 0.25) + 
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  labs(x = "sampling rate", y = "Area under the ROC curve")


## ----extract-subgraph, eval = FALSE-------------------------------------------
## blog_subgraph <- frenchblog2007 %>%
##   igraph::induced_subgraph(V(frenchblog2007)$party %in%
##   c( "right", "left"))
## blog_subgraph <- delete_vertices(blog_subgraph, 
##   which(degree(blog_subgraph) ==  0))


## ----plot-subgraph show, echo = TRUE, fig.keep = "none"-------------------------
plot(blog_subgraph, vertex.shape = "none", vertex.label = V(blog_subgraph)$party, 
  vertex.label.color = "steel blue", vertex.label.font = 1.5, 
  vertex.label.cex = .6, edge.color = "gray70", edge.width = 1)


## ----plot-subgraph, echo = TRUE, fig.fig.width = 5, fig.height = 5, out.width = "65%", fig.align = "center", fig.cap = "\\label{fig:subgraph}Subnetwork extracted from the French political blogosphere (only blogs with attribute party in $\\{\\mathrm{left, right}\\}$ were kept)", fig.pos = "htbp!", message = FALSE----
par_old <- par(mar = c(0, 0, 0, 0)+.1)
plot(blog_subgraph, vertex.shape = "none", vertex.label = V(blog_subgraph)$party, 
  vertex.label.color = "steel blue", vertex.label.font = 1.5, 
  vertex.label.cex = .6, edge.color = "gray70", edge.width = 1)
par(par_old)


## ----dummy-party, eval = FALSE------------------------------------------------
## dummy_party <- (V(blog_subgraph)$party = =  "left") * 1


## ----subgraph-sampling, eval = TRUE, message = FALSE, warning = FALSE-------------
blog_subgraph_obs <- blog_subgraph %>% as_adj() %>% 
  observeNetwork(sampling = "covar-node", parameters = 10, 
  covariates = list(dummy_party))
blocks <- 2:8
future::plan("multicore", workers = 10)


## ----covar-party-1, eval = FALSE----------------------------------------------
## sbm_covar1 <- estimateMissSBM(blog_subgraph_obs, blocks, 
##   "covar-node", covariates = list(dummy_party), 
##   control = list(useCov = TRUE, iterates = 2))


## ----covar-party-1-b show, echo = TRUE---------------------------------------
sbm_covar1$bestModel$fittedSampling$parameters
sbm_covar1$bestModel$fittedSBM$covarParam     


## ----covar-party-2, eval = FALSE----------------------------------------------
## sbm_covar2 <- estimateMissSBM(blog_subgraph_obs, blocks, 
##   "covar-node", covariates = list(dummy_party), 
##   control = list(useCov = FALSE, iterates = 2))


## ----covar-party-2-show, echo = TRUE-----------------------------------------
sbm_covar2$bestModel$fittedSampling$parameters 


## ----covar-party-3, eval = FALSE----------------------------------------------
## sbm_covar3 <- estimateMissSBM(blog_subgraph_obs, blocks, 
##   "node", covariates = list(dummy_party), 
##   control = list(useCov = TRUE, iterates = 2))


## ----covar-party-3-show, echo = TRUE-----------------------------------------
sbm_covar3$bestModel$fittedSampling$parameters 
sbm_covar3$bestModel$fittedSBM$covarParam      


## ----covar-party-4, eval = FALSE----------------------------------------------
## sbm_covar4 <- estimateMissSBM(blog_subgraph_obs, blocks, 
##   "node", control = list(useCov = FALSE, iterates = 2))


## ----covar-party-4-show, echo = TRUE-----------------------------------------
sbm_covar4$bestModel$fittedSampling$parameters


## ----plot-ICLs-covar, echo = TRUE, warning = FALSE, fig.cap = "\\label{fig:ICLcovar}ICL criterion for different numbers of blocks under the four models which make different use of the covariate political party.", fig.pos = "htbp!", fig.height = 4, fig.width = 6, out.width = ".6\\textwidth", fig.align = "center"----
rbind(tibble(Q = blocks, ICL = sbm_covar1$ICL, sampling = "covar-node", useCov = "true"), 
  tibble(Q = blocks, ICL = sbm_covar2$ICL, sampling = "covar-node", useCov = "false"), 
  tibble(Q = blocks, ICL = sbm_covar3$ICL, sampling = "node", useCov = "true"), 
  tibble(Q = blocks, ICL = sbm_covar4$ICL, sampling = "node", useCov = "false")
  ) %>% ggplot(aes(x = Q, y = ICL, color = sampling, shape = useCov))  +
  geom_line() + geom_point() + labs(x = "#blocks", y = "ICLs")


## ----covar-full, eval = FALSE-------------------------------------------------
## sbm_covar_full <- as_adj(blog_subgraph) %>%
##    estimateMissSBM(blocks, "node", covariates =  list(dummy_party))


## ----ARI-covar, echo = TRUE--------------------------------------------------
Blocks <- sbm_covar_full$bestModel$fittedSBM$memberships
Blocks1 <- sbm_covar1$bestModel$fittedSBM$memberships
Blocks2 <- sbm_covar2$bestModel$fittedSBM$memberships
Blocks3 <- sbm_covar3$bestModel$fittedSBM$memberships
Blocks4 <- sbm_covar4$bestModel$fittedSBM$memberships
ARIs <- matrix(c(aricode::ARI(Blocks, Blocks1), aricode::ARI(Blocks, Blocks2), 
  aricode::ARI(Blocks, Blocks3), aricode::ARI(Blocks, Blocks4), 
  NA, aricode::ARI(Blocks1, Blocks2), aricode::ARI(Blocks1, Blocks3), 
  aricode::ARI(Blocks1, Blocks4), 
  NA, NA, aricode::ARI(Blocks2, Blocks3), aricode::ARI(Blocks2, Blocks4), 
  NA, NA, NA, aricode::ARI(Blocks3, Blocks4)), nrow = 4, ncol = 4, byrow = TRUE)
ARIs <- as.data.frame(ARIs)
rownames(ARIs) <- c("Full", "i", "ii", "iii")
colnames(ARIs) <- c("i", "ii", "iii", "iv")

round(ARIs, 2)

# kableExtra::kbl(ARIs, digits = 2, booktabs = T,
#  caption = "\\label{tab:covarARIs}Clustering comparison with adjusted Rand indices (ARIs) between models taking the covariate into account (models $i)$, $ii)$, $iii)$, $iv)$ and model with fully observed network).") %>%
#  kableExtra::kable_styling(latex_options = c("hold_position"))

