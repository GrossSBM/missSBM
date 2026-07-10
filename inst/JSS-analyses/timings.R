library("missSBM")
library("parallel")
library("tidyverse")

data("frenchblog2007", package = "missSBM")
adjacencyMatrix <- igraph::as_adjacency_matrix(frenchblog2007, sparse = FALSE)
sbm_full <- estimateMissSBM(adjacencyMatrix, 10, "node", control = list(core = 10,
  iterates = 0))
samplingParameters <- base::sample(x = c(0.2, 0.8), size = sbm_full$bestModel$fittedSBM$nbBlocks,
  replace = TRUE)

vBlocks <- 1:14
nSim <- 50

res_timings_noNA <- mclapply(1:nSim, function(i)
  map_dbl(vBlocks, ~system.time(
    estimateMissSBM(adjacencyMatrix, ., "node", control = list(trace = 0, smoothing = "none")))["elapsed"]),
  mc.cores = 10) %>%
  setNames(paste0("sim-", 1:nSim)) %>%
  as_tibble() %>%
  add_column(`#group` = factor(vBlocks)) %>%
  pivot_longer(-`#group`, values_to = "time (sec.)", names_to = "#sim")

res_timings_MAR <- mclapply(1:nSim, function(i) {
  sampledNet <- observeNetwork(adjacencyMatrix = adjacencyMatrix, sampling = "node",
    parameters = 0.75)
  map_dbl(vBlocks, ~system.time(
    estimateMissSBM(sampledNet, ., "node", control = list(trace = 0, smoothing = "none")))["elapsed"])
}, mc.cores = 10) %>%
  setNames(paste0("sim-", 1:nSim)) %>%
  as_tibble() %>%
  add_column(`#group` = factor(vBlocks)) %>%
  pivot_longer(-`#group`, values_to = "time (sec.)", names_to = "#sim")

res_timings_MNAR <- mclapply(1:nSim, function(i) {
  sampledNet <- observeNetwork(adjacencyMatrix = adjacencyMatrix, sampling = "block-node",
    parameters = samplingParameters, clusters = sbm_full$bestModel$fittedSBM$memberships)
  map_dbl(vBlocks, ~system.time(
    estimateMissSBM(sampledNet, ., "node", control = list(trace = 0, smoothing = "none")))["elapsed"])
}, mc.cores = 10) %>%
  setNames(paste0("sim-", 1:nSim)) %>%
  as_tibble() %>%
  add_column(`#group` = factor(vBlocks)) %>%
  pivot_longer(-`#group`, values_to = "time (sec.)", names_to = "#sim")

res_timings_noNA$missingness <- "none"
res_timings_MAR$missingness <- "MAR"
res_timings_MNAR$missingness <- "MNAR"

res_timings <- rbind(res_timings_noNA, res_timings_MAR, res_timings_MNAR)

save(res_timings, file = "timings.RData")
