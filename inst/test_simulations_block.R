rm(list=ls())
set.seed(12345)

library(missSBM)
library(aricode)

old_rep <- getwd()
setwd("../svn_oldies/Code/code_these/functions/")
source("func_missSBM.class.R")
source("func_missSBM.R")
setwd(old_rep)

## ----------------------------------------------
## SIMULATION

# SBM parameters
N  <- 300
pr <- .2
pi <- matrix(pr,2,2)
pi[1,1] <- 3*pr
diri <- FALSE

# sampling parameters
a    <- 1/2
rho  <- .7
samp <- "block"
parameters <- c(rho*a,rho)

alpha1 <- .15
alpha  <- c(alpha1,alpha1*(3*a*rho-2*a^2*rho^2-rho+a*rho^2)/(-a*rho+a*rho^2+rho-rho^2))

control <- list(threshold = 1e-3, maxIter = 200, fixPointIter = 10, trace = TRUE)

## simuler un nouveau réseau
Net            <- simulateSBM(N,alpha,pi,diri)
cl_star        <- as.vector(Net$blocks%*%1:2)
Netsamp        <- samplingSBM(Net$adjacencyMatrix,samp,parameters,clusters = cl_star)
sampledNet     <- Netsamp$adjacencyMatrix
samplingVector <- rep(0,N); samplingVector[!is.na(rowSums(sampledNet))] <- 1

# ancien code -------------------------------------------------------------
block_old <- func_missSBM.class(sampledNet, seq.Q = 2, cl.init = "CAH")@models[[1]]

# package -----------------------------------------------------------------
mar <-
  inferSBM(
    adjacencyMatrix = sampledNet,
    vBlocks = 2,
    sampling = "node",
    clusterInit = "hierarchical")$models[[1]]

new <-
    inferSBM(
      adjacencyMatrix = sampledNet,
      vBlocks = 2,
      sampling = "block",
      clusterInit = "hierarchical", control_VEM = control)

block_new <- new$models[[1]]

res <- data.frame(version  = c("new","new", "old"),
          sampling  = c("MAR", "block", "block"),
          ARI       = c(ARI(cl_star, mar$fittedSBM$memberships),
                        ARI(cl_star, block_new$fittedSBM$memberships),
                        ARI(cl_star, getClusters(block_old))))
print(res)

plot(-new$monitor$objective, type = "l", log = "y")
plot(-block_old@crit, type = "l", log = "y")

psi_new <- block_new$fittedSampling$parameters
psi_old <- block_old@psi

pi_new <- block_new$fittedSBM$connectParam
pi_old <- block_old@theta[[1]]$pi

