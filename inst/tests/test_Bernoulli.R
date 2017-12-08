library(missSBM)

#### DEFINE A SBM MODEL
alpha <- c(1/3, 1/3, 1/3)
pi <- diag(.05,3) + .03
n <- 300
mySBM <- SBM_BernoulliDirected$new(n, alpha, pi)

## SAMPLE SOME NETWORK DATA
SBMdata       <- mySBM$rSBM() # full graph
sampling_rate <- .01
mySampled  <- sampling_snowball$new(n, sampling_rate, TRUE)
sample     <- mySampled$rSampling(SBMdata$adjacencyMatrix)

sbm <- SBM_collection$new(sample$adjacencyMatrix, 3, "snowball", "Bernoulli", TRUE)

# plot(sbm$vICLs, type = "l")

# Xmis <- sample$adjacencyMatrix; Xmis[is.na(Xmis)] <- 0
# a <- func.ICL2_SBM_Poisson(Xmis,sample$samplingMatrix,qmin=1,qmax = 10)
# plot(a$ICL, type = "l")
