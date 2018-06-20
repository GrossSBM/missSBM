library(missSBM)

n <- 100
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
sbm <- inferSBM(sampledAdjMatrix$adjacencyMatrix, vBlocks, sampling)

# ICL
ICL <- sapply(sbm$models, function(x) x$vICL)
plot(ICL)

# Error
frobenius(sbm$models[[3]]$fittedSampling$parameters - samplingParameters)/frobenius(sbm$models[[3]]$fittedSampling$parameters)
