library(missSBM)
library(greed)
library(aricode)
library(future.apply)
library(tidyverse)
library(viridis)

## network structure: community
## network parameters
n <- 300
pi <-  c(1/3,1/3,1/3)
lambda  <- matrix(.2,3,3) + diag(3) * (.3)
theta   <- list(mean = lambda)

## block-dyad sampling

## dyad (MAR) sampling
## Simulation
nbrSimu <- 20

default_ctrl <-
  list(threshold = 0.01, trace = FALSE, imputation = "median",
      similarity = missSBM:::l1_similarity, useCov = TRUE, maxIter = 50,
      fixPointIter = 3, iterates = 1, exploration = "both",
      clusterInit = NULL)

options(future.fork.enable = TRUE)
future::plan("multicore", workers = 10)

list_psi_mar <- seq(0.3, 1, by = 0.1)
RES <- lapply(1:nbrSimu, function(i) {
  res <- future.apply::future_lapply(list_psi_mar, function(psi) {

  A     <- sbm::sampleSimpleSBM(nbNodes = n, blockProp = pi, connectParam = theta, model = "bernoulli" )
  A_obs <- missSBM::observeNetwork(A$networkData, sampling="dyad", parameters = psi)

  sbm_mar   <- estimateMissSBM(A_obs, 1:6, "node", control = list(iterates = 0, trace = 0))
  missSBM_i0 <- sbm_mar$bestModel
  sbm_mar$explore(default_ctrl)
  missSBM_i1 <- sbm_mar$bestModel
  sbm_mar$explore(default_ctrl)
  missSBM_i2 <- sbm_mar$bestModel

  greed_mar <- greed(A_obs, model = new("misssbm", type = 'undirected', sampling = 'dyad'))

  res <- data.frame(
    ARI = c(ARI(missSBM_i0$fittedSBM$memberships, A$memberships),
            ARI(missSBM_i1$fittedSBM$memberships, A$memberships),
            ARI(missSBM_i2$fittedSBM$memberships, A$memberships),
            ARI(greed_mar@cl, A$memberships)),
    choice = c(missSBM_i0$fittedSBM$nbBlocks,
               missSBM_i1$fittedSBM$nbBlocks,
               missSBM_i2$fittedSBM$nbBlocks,
               greed_mar@K),
    variant = c("ICL_0", "ICL_1", "ICL_2", "greed")
  )
  res$simu <- i
  res$psi  <- psi
  res
  }) %>% purrr::reduce(rbind)
  res
}) %>% reduce(rbind)

p <- RES %>% pivot_longer(c(ARI, choice), names_to = "measure") %>%
  group_by(variant) %>%
  ggplot() + theme_bw() + viridis::scale_colour_viridis(discrete = TRUE) + viridis::scale_fill_viridis(discrete = TRUE) +
  aes(x = factor(psi), y = value, fill = variant) +
  geom_boxplot() + geom_point() + facet_wrap(~ measure, scales = "free")
