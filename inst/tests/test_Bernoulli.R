library(missSBM)

#### DEFINE A SBM MODEL
alpha <- c(1/3, 1/3, 1/3)
pi <- diag(.45,3) + .05
n <- 100
mySBM <- SBM_BernoulliUndirected$new(n, alpha, pi)

## SAMPLE SOME NETWORK DATA
SBMdata       <- mySBM$rSBM() # full graph
sampling_rate <- 0.5
mySampled  <- sampling_randomPairMAR$new(n, sampling_rate, FALSE)
sample     <- mySampled$rSampling(SBMdata$adjacencyMatrix)

sbm <- SBM_collection$new(sample$adjacencyMatrix, 1:8, "MAREdge", "Bernoulli", FALSE)

plot(sbm$vICLs, type = "l")

