library(missSBM)
library(greed)
library(aricode)
library(future.apply)
library(tidyverse)
library(viridis)

## network structure: community
## network parameters
n <- 200
pi <-  c(1/3,1/3,1/3)
lambda  <- matrix(.25,3,3) + diag(3) * (.15)
theta <- list(mean = lambda)

## block-dyad sampling

## dyad (MAR) sampling

## double-standard sampling

