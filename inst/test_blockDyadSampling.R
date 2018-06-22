rm(list = ls())
library(missSBM)
library(mclust)

n <- 200
Q <- 3
alpha <- rep(1,Q)/Q
pi <- diag(.45,Q) + .05
directed <- FALSE

sbm <- simulateSBM(n, alpha, pi, directed)
adjacencyMatrix <- sbm$adjacencyMatrix
samplingParameters <- diag(.3,Q) + .3
sampling <- "block_dyad"
sampledAdjMatrix <- samplingSBM(adjacencyMatrix, sampling, samplingParameters, clusters = sbm$memberships)

# Sampling rate
(sum(sampledAdjMatrix$samplingMatrix)-n)/(n*(n-1))

# Inference
vBlocks <- 1:5
infer <- inferSBM(adjacencyMatrix, vBlocks, sampling)

# ICL
ICL <- sapply(infer$models, function(x) x$vICL)
plot(ICL)

# Error
sum((infer$models[[3]]$fittedSampling$parameters - samplingParameters)^2)/sum((infer$models[[3]]$fittedSampling$parameters)^2)
adjustedRandIndex(sbm$memberships, infer$models[[3]]$fittedSBM$memberships)

