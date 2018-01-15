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
my_sampling_dyad <- sampling_model$new("dyad", 0.1)
A_dyad <- mySampling_Dyad$rSampling(mySBM$adjacencyMatrix)

my_sampling_node <- sampling_model$new("node", .3)
A_node <- my_sampling_node$rSampling(mySBM$adjacencyMatrix)

my_sampling_block <- sampling_model$new("block", .3)
A_block <- my_sampling_block$rSampling(mySBM$adjacencyMatrix, mySBM$clusters)

my_sampling_double_Standard <- sampling_model$new("double_standard", c(0.1, 0.5))
A_double_standard <- my_sampling_double_Standard$rSampling(mySBM$adjacencyMatrix)

my_sampling_degree <- sampling_model$new("degree", c(-.1,.2))
A_degree <- my_sampling_degree$rSampling(mySBM$adjacencyMatrix)

my_sampling_snowball <- sampling_model$new("snowball", .3)
A_snowball <- my_sampling_snowball$rSampling(mySBM$adjacencyMatrix)

samplingSBM(adjacencyMatrix, "dyad", 0.1)
## Sampling of the data : ##
# samplingParameters <- .5                                                           # the sampling rate
# sampling <- "MAREdge"                                                              # the sampling design
# sampledAdjMatrix <- samplingSBM(adjacencyMatrix, sampling, samplingParameters)     # the sampled adjacency matrix


# ## Inference :
# vBlocks <- 1:5                                                                     # number of classes
# sbm <- inferSBM(sampledAdjMatrix, vBlocks, sampling)             # the inference
