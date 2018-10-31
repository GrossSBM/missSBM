rm(list=ls())
setwd("~/work/students/timothee/svn_oldies/Code/code_these/functions/")

set.seed(12345)

library(devtools)
library(missSBM)
#ancien code
source("func_missSBM.class.R")
source("func_missSBM.R")

library(igraph)
library(missSBM)
library(Matrix)
library(ape)

## ----------------------------------------------
## SIMULATION

# SBM parameters
N = 300
pr = .2
pi = matrix(pr,2,2)
pi[1,1]=3*pr
diri = FALSE
#pi[2,2]=.04

# sampling parameters
a=1/2
rho = .8
samp = "block"
parameters = c(rho*a,rho)

alpha1=.15
alpha = c(alpha1,alpha1*(3*a*rho-2*a^2*rho^2-rho+a*rho^2)/(-a*rho+a*rho^2+rho-rho^2))
alpha*N

## simuler un nouveau reseau
Net = simulateSBM(N,alpha,pi,diri)
Netsamp <- samplingSBM(Net$adjacencyMatrix,samp,parameters,clusters = Net$blocks%*%1:2)
sampledNet <- Netsamp$adjacencyMatrix
samplingVector <- rep(0,N); samplingVector[!is.na(rowSums(sampledNet))] <- 1
cl_star <- as.vector(Net$blocks%*%1:2) ## true classifaction

## ----------------------------------------------
## CAH
CAH.init <- graphCAH(sampledNet,2)
table(CAH.init,samplingVector)
table(CAH.init,as.vector(Net$blocks%*%1:2))

# ancien code -------------------------------------------------------------
outPut_class_old <- func_missSBM.class(sampledNet, seq.Q = 2, cl.init = "CAH")
getPsi(getModel(outPut_class_old, 2))
old_VEM <- getModel(outPut_class_old, 2)

# package -----------------------------------------------------------------
outPut_mar  <-
  inferSBM(
    adjacencyMatrix = sampledNet,
    vBlocks = 2,
    sampling = "node",
    clusterInit = "hierarchical")

outPut_block  <-
  inferSBM(
    adjacencyMatrix = sampledNet,
    vBlocks = 2,
    sampling = "block",
    clusterInit = "hierarchical")

mar_missSBM    <- outPut_mar$models[[1]]
block_missSBM  <- outPut_block$models[[1]]

## performance
aricode::ARI(cl_star, CAH.init)
aricode::ARI(cl_star, mar_missSBM$fittedSBM$memberships)
aricode::ARI(cl_star, getClusters(old_VEM)  )
aricode::ARI(cl_star, block_missSBM$fittedSBM$memberships)

