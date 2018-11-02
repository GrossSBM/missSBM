rm(list=ls())
setwd("~/GitLab/svn_oldies/Code/code_these/functions/")

library(parallel)
library(pbmcapply)
library(devtools)
library(missSBM)
library(igraph)
library(Matrix)
library(ape)

source("func_missSBM.class.R")
source("func_missSBM.R")
source("~/Git/missSBM/inst/initializations.R")

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

# fonctions à reproduire
simulations <- function(i){

  ## simuler un nouveau réseau
  Net            <- simulateSBM(N,alpha,pi,diri)
  cl_star        <- as.vector(Net$blocks%*%1:2)
  Netsamp        <- samplingSBM(Net$adjacencyMatrix,samp,parameters,clusters = cl_star)
  sampledNet     <- Netsamp$adjacencyMatrix
  samplingVector <- rep(0,N); samplingVector[!is.na(rowSums(sampledNet))] <- 1

  # Initialisation spectral clustering pour réseau sparse -------------------
  Init_SpClSparse <- SparseSpectralClustering_NAisMean(sampledNet,2)

  # ancien code -------------------------------------------------------------
  outPut_class_old <- func_missSBM.class(sampledNet, seq.Q = 2, cl.init = "CAH")
  outPut_class_old_SpClSparse <- func_missSBM.class(sampledNet, seq.Q = 2, CL.init = Init_SpClSparse)

  # package -----------------------------------------------------------------
  outPut_mar <-
    inferSBM(
      adjacencyMatrix = sampledNet,
      vBlocks = 2,
      sampling = "node",
      clusterInit = "hierarchical",
      smoothing = "none")

  outPut_mar_smoothed <-
    inferSBM(
      adjacencyMatrix = sampledNet,
      vBlocks = 2,
      sampling = "node",
      clusterInit = Init_SpClSparse,
      smoothing = "none",
      iter_both = 2)

  outPut_block <-
    inferSBM(
      adjacencyMatrix = sampledNet,
      vBlocks = 2,
      sampling = "block",
      clusterInit = "hierarchical",
      smoothing = "none")

  outPut_block_smoothed <-
    inferSBM(
      adjacencyMatrix = sampledNet,
      vBlocks = 2,
      sampling = "block",
      clusterInit = Init_SpClSparse,
      smoothing = "none",
      iter_both = 2)

  ICL_mar                     <- sapply(outPut_mar$models  , function(model) model$vICL)
  ICL_mar_smoothed            <- sapply(outPut_mar_smoothed$models  , function(model) model$vICL)
  ICL_block_missSBM           <- sapply(outPut_block$models, function(model) model$vICL)
  ICL_block_missSBM_smoothed  <- sapply(outPut_block_smoothed$models, function(model) model$vICL)

  mar_missSBM             <- outPut_mar$models[[which.min(ICL_mar)]]
  mar_missSBM_smoothed    <- outPut_mar_smoothed$models[[which.min(ICL_mar_smoothed)]]
  block_missSBM           <- outPut_block$models[[which.min(ICL_block_missSBM)]]
  block_missSBM_smoothed  <- outPut_block_smoothed$models[[which.min(ICL_block_missSBM_smoothed)]]
  block_old               <- getBestModel(outPut_class_old)
  block_old_SpClSparse    <- getBestModel(outPut_class_old_SpClSparse)

  return(data.frame(simu  = i,
                    algo  = c("MAR",
                             "MAR_SpClSparse",
                             "block_new",
                             "block_new_SpClSparse",
                             "block_old",
                             "block_old_SpClSparse"),
                    ICL   = c(ICL_mar[which.min(ICL_mar)],
                             ICL_mar_smoothed[which.min(ICL_mar_smoothed)],
                             ICL_block_missSBM[which.min(ICL_block_missSBM)],
                             ICL_block_missSBM_smoothed[which.min(ICL_block_missSBM_smoothed)],
                             getBestModel(outPut_class_old)@ICL,
                             getBestModel(outPut_class_old_SpClSparse)@ICL),
                    Q_est = c(which.min(ICL_mar), which.min(ICL_mar_smoothed), which.min(ICL_block_missSBM), which.min(ICL_block_missSBM_smoothed), length(getAlpha(getBestModel(outPut_class_old))),length(getAlpha(getBestModel(outPut_class_old_SpClSparse)))),
                    ARI   = c(aricode::ARI(cl_star, mar_missSBM$fittedSBM$memberships),
                             aricode::ARI(cl_star, mar_missSBM_smoothed$fittedSBM$memberships),
                             aricode::ARI(cl_star, block_missSBM$fittedSBM$memberships),
                             aricode::ARI(cl_star, block_missSBM_smoothed$fittedSBM$memberships),
                             aricode::ARI(cl_star, getClusters(block_old)),
                             aricode::ARI(cl_star, getClusters(block_old_SpClSparse))),
                    sampRate = Netsamp$samplingRate))

}

results <- do.call(rbind, pbmclapply(1:50, simulations, mc.cores = 1))
results$interval <- cut(as.numeric(results$sampRate), seq(min(as.numeric(results$sampRate))-0.001,max(as.numeric(results$sampRate)), len=2))

library(ggplot2)
ggplot(results, aes(x=interval, y=ARI, color=algo)) + geom_boxplot()
save(results, file = "~/Git/missSBM/inst/Rdata_simu3.RData")
