rm(list=ls())

library(parallel)
library(pbmcapply)
library(devtools)
library(missSBM)
library(igraph)
library(Matrix)
library(ape)
source("inst/initializations.R")

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

control <- list(threshold = 1e-3, maxIter = 100, fixPointIter = 10, trace = FALSE)

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
  outPut_class_old_hierarchical <- func_missSBM.class(sampledNet, seq.Q = 2, cl.init = "CAH")
  outPut_class_old_spectral     <- func_missSBM.class(sampledNet, seq.Q = 2, CL.init = Init_SpClSparse)

  # package -----------------------------------------------------------------
  outPut_mar_hierachical <-
    inferSBM(
      adjacencyMatrix = sampledNet,
      vBlocks = 2,
      sampling = "node",
      clusterInit = "hierarchical",
      smoothing = "none")

  outPut_mar_spectral <-
    inferSBM(
      adjacencyMatrix = sampledNet,
      vBlocks = 2,
      sampling = "node",
      clusterInit = Init_SpClSparse,
      smoothing = "none")

  outPut_block_hierarchical <-
    inferSBM(
      adjacencyMatrix = sampledNet,
      vBlocks = 2,
      sampling = "block",
      clusterInit = "hierarchical",
      smoothing = "none", control_VEM = control)

  outPut_block_spectral <-
    inferSBM(
      adjacencyMatrix = sampledNet,
      vBlocks = 2,
      sampling = "block",
      clusterInit = Init_SpClSparse,
      smoothing = "none", control_VEM = control)

  mar_missSBM_hierarchical   <- outPut_mar_hierachical$models[[1]]
  mar_missSBM_spectral       <- outPut_mar_spectral$models[[1]]
  block_missSBM_hierarchical <- outPut_block_hierarchical$models[[1]]
  block_missSBM_spectral     <- outPut_block_spectral$models[[1]]

  ICL_mar_hierarchical   <- outPut_mar_hierachical$models[[1]]$vICL
  ICL_mar_spectral       <- outPut_mar_spectral$models[[1]]$vICL
  ICL_block_hierarchical <- outPut_block_hierarchical$models[[1]]$vICL
  ICL_block_spectral     <- outPut_block_spectral$models[[1]]$vICL

  block_old_hierarchical <- getBestModel(outPut_class_old_hierarchical)
  block_old_spectral     <- getBestModel(outPut_class_old_spectral)

  return(data.frame(simu  = i,
                    algo  = c("MAR_hierarchical",
                             "MAR_spectral",
                             "block_new_hierarchical",
                             "block_new_spectral",
                             "block_old_hierarchical",
                             "block_old_spectral"),
                    version  = c("new",
                             "new",
                             "new",
                             "new",
                             "old",
                             "old"),
                    sampling  = c("MAR",
                             "MAR",
                             "block",
                             "block",
                             "block",
                             "block"),
                    init  = c("hierarchical",
                             "spectral",
                             "hierarchical",
                             "spectral",
                             "hierarchical",
                             "spectral"),
                    ICL   = c(ICL_mar_hierarchical,
                             ICL_mar_spectral,
                             ICL_block_hierarchical,
                             ICL_block_spectral,
                             getBestModel(outPut_class_old_hierarchical)@ICL,
                             getBestModel(outPut_class_old_spectral)@ICL),
                    ARI   = c(aricode::ARI(cl_star, mar_missSBM_hierarchical$fittedSBM$memberships),
                             aricode::ARI(cl_star, mar_missSBM_spectral$fittedSBM$memberships),
                             aricode::ARI(cl_star, block_missSBM_hierarchical$fittedSBM$memberships),
                             aricode::ARI(cl_star, block_missSBM_spectral$fittedSBM$memberships),
                             aricode::ARI(cl_star, getClusters(block_old_hierarchical)),
                             aricode::ARI(cl_star, getClusters(block_old_spectral))),
                    sampRate = Netsamp$samplingRate))

}

n_sim <- 50
results <- do.call(rbind, pbmclapply(1:n_sim, simulations, mc.cores = 10))

library(ggplot2)
p <- ggplot(results, aes(x = sampling, y = ARI, color = init, fill = version)) + geom_boxplot()
print(p)
