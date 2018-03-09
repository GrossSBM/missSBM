rm(list=ls())
library(missSBM)
library(pbmcapply)

n <- 200
Q <- 4
pi <- diag(0.45,Q) + .05
alpha <- rep(1,Q)/Q
family <- "Bernoulli"
directed <- FALSE
mySBM <- simulateSBM(n, alpha, pi, directed)
adjacencyMatrix <- mySBM$adjacencyMatrix
samplingParameters <- .3
sampling <- "dyad"
sampSBM <- samplingSBM(adjacencyMatrix, sampling, samplingParameters)
sampledAdjMatrix <- sampSBM$adjacencyMatrix
vBlocks <- 1:8
sbm <- inferSBM(sampledAdjMatrix, vBlocks, "dyad", nIterSmoothingICL = 1)

# smoothed_fwrd_half <- smoothingForward_half(sbm$models, vBlocks, sampledAdjMatrix, sampling)
# smoothed_fwrd_SpCl <- smoothingForward_SpCl(sbm$models, vBlocks, sampledAdjMatrix, sampling)
# smoothed_back      <- smoothingBackward(sbm$models, vBlocks, sampledAdjMatrix, sampling)
# smoothed_fb_half   <- smoothingForBackWard_half(sbm$models, vBlocks, sampledAdjMatrix, sampling)
# smoothed_fb_SpCl   <- smoothingForBackWard_SpCl(sbm$models, vBlocks, sampledAdjMatrix,sampling)

vICL           <- sapply(sbm$models, function(model) model$vICL); plot(vICL, type = "l")
# vICL_back      <- sapply(smoothed_back, function(model) model$vICL); lines(vICL_back, col="blue")
# vICL_fwrd_half <- sapply(smoothed_fwrd_half, function(model) model$vICL); lines(vICL_fwrd_half, col="red")
# vICL_fwrd_SpCl <- sapply(smoothed_fwrd_SpCl, function(model) model$vICL); lines(vICL_fwrd_SpCl, col="yellow")
# vICL_fb_half   <- sapply(smoothed_fb_half, function(model) model$vICL); lines(vICL_fb_half  , col="orange")
# vICL_fb_SpCl   <- sapply(smoothed_fb_SpCl, function(model) model$vICL); lines(vICL_fb_SpCl  , col="green")
