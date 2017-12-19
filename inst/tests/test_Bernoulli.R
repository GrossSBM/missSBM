library(missSBM)

#### DEFINE A SBM MODEL
# alpha <- c(1/3, 1/3, 1/3)
# pi <- diag(.22,3) + .05
n <- 300

## Another SBM model :
Q <- 3
alpha <- rep(1,3)/3
pir <- .05
pia <- 5*pir
pi <- matrix(c(pia, pir ,pir, pia, pir, pir, pir, pia, pir),3,3)

mySBM <- SBM_BernoulliDirected$new(n, alpha, pi)

## SAMPLE SOME NETWORK DATA
SBMdata       <- mySBM$rSBM() # full graph
sampling_rate <- 1
mySampled  <- sampling_randomNodesMAR$new(n, sampling_rate, TRUE)
sample     <- mySampled$rSampling(SBMdata$adjacencyMatrix)

# adj_test <- (sample$adjacencyMatrix | t(sample$adjacencyMatrix))*1

sbm <- SBM_collection$new(sample$adjacencyMatrix, 3, "MARNode", "Bernoulli", TRUE)

t(alpha) %*% pi %*% alpha
sbm$models[[1]]$sampledNetwork$samplingRate
adjustedRandIndex(apply(sbm$models[[1]]$blockVarParam, 1, which.max), SBMdata$blocks %*% 1:3)

