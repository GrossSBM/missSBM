library(missSBM)

#### DEFINE A SBM MODEL
alpha <- c(1/3, 1/3, 1/3)
pi <- diag(.45,3) + .05
n <- 100
mySBM <- SBM_BernoulliUndirected$new(n, alpha, pi)

## SAMPLE SOME NETWORK DATA
SBMdata       <- mySBM$rSBM() # full graph
sampling_rate <- 0.5
mySampled  <- sampling_snowball$new(n, rep(sampling_rate, n), FALSE)
sample     <- mySampled$rSampling(SBMdata$adjacencyMatrix)

sbm <- SBM_collection$new(sample$adjacencyMatrix, 1:10, "MAREdge", "Poisson", FALSE)

plot(sbm$vICLs, type = "l")

# Xmis <- sample$adjacencyMatrix; Xmis[is.na(Xmis)] <- 0
# a <- func.ICL2_SBM_Poisson(Xmis,sample$samplingMatrix,qmin=1,qmax = 10)
# plot(a$ICL, type = "l")
