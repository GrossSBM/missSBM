library(sbm)
library(missSBM)
library(future.apply)
library(tidyverse)
library(viridis)

## network parameters
n <- 200
pi <-  c(1/3,1/3,1/3)
lambda  <- matrix(.25,3,3) + diag(3) * (.15)
theta <- list(mean = lambda)

## sampling parameters
rho_min <- 0.1
rho_max <- 0.9
delta <- seq(0, rho_max-rho_min, 0.1)
list_psi <- map(delta, ~matrix(rho_min + .x, 3, 3) + diag(3) * (rho_max-rho_min - 2 * .x))

control <- list(threshold = 1e-2, maxIter = 50, fixPointIter = 3, trace = 0)

## Simulation
nbrSimu <- 100

future::plan("multisession", workers = 10)

RES <- lapply(1:nbrSimu, function(i) {
  res <- future.apply::future_lapply(list_psi, function(psi) {

  A <- sbm::sampleSimpleSBM(nbNodes = n, blockProp = pi, connectParam = theta, model = "bernoulli" )

  A_obs <-
    missSBM::observeNetwork(A$networkData, sampling="block-dyad", clusters = A$memberships, parameters = psi) %>%
    missSBM:::partlyObservedNetwork$new()

  cl_init <- A_obs$clustering(1:3)[[3]]

  MAR <- missSBM_fit$new(A_obs, "dyad", cl_init)
  out <- MAR$doVEM(control)

  MNAR <- missSBM_fit$new(A_obs, "block-dyad", cl_init)
  out <- MNAR$doVEM(control)

  res <- data.frame(
    ARI = c(aricode::ARI(cl_init,A$memberships), aricode::ARI(MNAR$fittedSBM$memberships,A$memberships), aricode::ARI(MAR$fittedSBM$memberships,A$memberships)),
    MSE = c(NA, sqrt(sum((MNAR$fittedSBM$connectParam$mean-theta$mean)^2)), sqrt(sum((MAR$fittedSBM$connectParam$mean-theta$mean)^2))),
    Psi = c(NA, sqrt(sum((MNAR$fittedSampling$parameters-psi)^2)), NA),
    variant = c("init", "MNAR", "MAR")
  )
  res$simu  <- i
  res$delta <- psi[2,1] - psi[1,1]
  res
  }) %>% purrr::reduce(rbind)
  res
}) %>% reduce(rbind)


p <- RES %>% dplyr::select(-Psi) %>% pivot_longer(c(ARI, MSE), names_to = "measure") %>%
  group_by(variant) %>%
  ggplot() + theme_bw() + viridis::scale_colour_viridis(discrete = TRUE) + viridis::scale_fill_viridis(discrete = TRUE) +
  aes(x = factor(delta), y = value, fill = variant) +
  geom_boxplot() + geom_point() + facet_wrap(~ measure, scales = "free")

print(p)
