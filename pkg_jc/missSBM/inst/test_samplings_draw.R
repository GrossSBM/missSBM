library(missSBM)

### A SBM model : ###
n <- 300
Q <- 3
alpha <- rep(1,Q)/Q                                                                # mixture parameter
pi <- diag(.45,Q) + .05                                                            # connectivity matrix
family <- "Bernoulli"                                                              # the emmission law
directed <- FALSE                                                                  # if the network is directed or not

### Draw a SBM model
mySBM <- simulateSBM(n, alpha, pi, family, directed)                               # simulation of ad Bernoulli non-directed SBM

adjacencyMatrix <- mySBM$adjacencyMatrix                                           # the adjacency matrix

## Draw
dyad <- samplingSBM(adjacencyMatrix, "dyad", 0.1)
par(mfrow = c(1,2))
image_NA(dyad$samplingMatrix , main = "sampling matrix")
image_NA(dyad$adjacencyMatrix, main = "adjacency matrix")
title(main = paste("\ndyad sampling, sampling rate:", signif(dyad$samplingRate,3)), outer = TRUE)

node <- samplingSBM(adjacencyMatrix, "node", 0.1)
par(mfrow = c(1,2))
image_NA(node$samplingMatrix , main = "sampling matrix")
image_NA(node$adjacencyMatrix, main = "adjacency matrix")
title(main = paste("\nnode sampling, sampling rate:", signif(node$samplingRate,3)), outer = TRUE)

block <- samplingSBM(adjacencyMatrix,"block", c(.1, .2, .7), mySBM$clusters)
par(mfrow = c(1,2))
image_NA(block$samplingMatrix[order(mySBM$clusters),order(mySBM$clusters)] , main = "sampling matrix")
image_NA(block$adjacencyMatrix[order(mySBM$clusters),order(mySBM$clusters)], main = "adjacency matrix")
title(main = paste("\nblock sampling, sampling rate:", signif(block$samplingRate,3)), outer = TRUE)

double_standard <- samplingSBM(adjacencyMatrix,"double_standard", c(0.1, 0.5))
par(mfrow = c(1,2))
image_NA(double_standard$samplingMatrix , main = "sampling matrix")
image_NA(double_standard$adjacencyMatrix, main = "adjacency matrix")
title(main = paste("\ndouble_standard sampling, sampling rate:", signif(double_standard$samplingRate,3)), outer = TRUE)

degree <- samplingSBM(adjacencyMatrix,"degree", c(0.01,0.01))
par(mfrow = c(1,2))
image_NA(degree$samplingMatrix , main = "sampling matrix")
image_NA(degree$adjacencyMatrix, main = "adjacency matrix")
title(main = paste("\ndegree sampling, sampling rate:", signif(degree$samplingRate,3)), outer = TRUE)

snowball <- samplingSBM(adjacencyMatrix,"snowball", .3)
par(mfrow = c(1,2))
image_NA(snowball$samplingMatrix , main = "sampling matrix")
image_NA(snowball$adjacencyMatrix, main = "adjacency matrix")
title(main = paste("\nsnowball sampling, sampling rate:", signif(snowball$samplingRate, 3)), outer = TRUE)
