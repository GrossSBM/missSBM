library(missSBM)
library(sbm)
library(aricode)
library(viridis)
library(tidyverse)
library(aricode)
library(profvis)

## SBM parameters
N <- 2000 # number of nodes
Q <- 5   # number of clusters
pi <- rep(1,Q)/Q     # block proportion
theta <- list(mean = diag(.25,Q) + .1 ) # connectivity matrix

## generate a undirected binary SBM with no covariate
sbm <- sbm::sampleSimpleSBM(N, pi, theta)

control <- list(threshold = 1e-4, maxIter = 50, fixPointIter = 5, trace = 0)

psi <- seq(from = 1, to = 0.2, by = -0.1)
timings  <- vector("numeric", length(psi))
accuracy <- vector("numeric", length(psi))
iterates <- vector("numeric", length(psi))

res <- lapply(1:20, function(k) {
  for (i in seq_along(psi) ) {
    A_obs <-
      missSBM::observeNetwork(sbm$networkData, sampling="dyad", parameters = psi[i]) %>%
      missSBM:::partlyObservedNetwork$new()

    timings[i] <- system.time({
      cl_init <- A_obs$clustering(5)[[1]]
      MAR    <- missSBM_fit$new(A_obs, "dyad", cl_init)
      out    <- MAR$doVEM(control)
      cl_MAR <- MAR$fittedSBM$memberships
    })[[3]]
    accuracy[i] <- aricode::ARI(cl_MAR, sbm$memberships)
    iterates[i] <- length(out$iteration)
  }
  data.frame(
          psi = psi,
          accuracy = accuracy,
          timings = timings,
          iterates = iterates, simu = k)
}) %>% reduce(rbind)

ggplot(res) + aes(x = iterates, y = accuracy, group = simu, color = factor(psi)) +
  geom_line() + geom_point() + scale_color_viridis(discrete = TRUE) + theme_bw()

