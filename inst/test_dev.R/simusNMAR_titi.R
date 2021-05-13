
library(sbm)
library(missSBM)
library(future.apply)
library(dplyr)
library(tidyr)
library(purrr)
library(viridis)
library(greed)

## network parameters
n <- 200
pi <-  c(1/3,1/3,1/3)
epsilon <- .25
delta   <- .4
lambda  <- matrix(.25,3,3) + diag(3) * (.15)
theta <- list(mean = lambda)

## sampling parameters
rho_min <- 0.1
rho_max <- 0.9
delta <- seq(0, rho_max-rho_min, 0.1)
list_psi <- map(delta, ~matrix(rho_min + .x, 3, 3) + diag(3) * (rho_max-rho_min - 2 * .x))

psi=list_psi[[1]]


A <- sbm::sampleSimpleSBM(nbNodes = n, blockProp = pi, connectParam = theta, model = "bernoulli" )

A_obs <-
  missSBM::observeNetwork(A$networkData, sampling="block-dyad", clusters = A$memberships, parameters = psi) 


sol=greed(A_obs)


sol=greed(A_obs,model=new("missbm",type="undirected"),K=25)