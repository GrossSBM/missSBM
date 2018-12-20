### ====================== ###
### Tests anciens codes :  ###
### ====================== ###
rm(list=ls())
if(Sys.info()["user"] == "tabouy") {
  setwd("~/GitLab/svn_oldies/Code/code_these/functions/")
}
source("func_missSBM.R")
source("func_missSBM.twoStd.R")
source("func_missSBM.class.R")
source("func_missSBM.degree.R")
source("func_sampling.R")

## A SBM model : ##
N <- 100
Q <- 3
alpha <- c(0.25,0.5,0.25) # rep(1,Q)/Q     # mixture parameter
pi <- diag(.45,Q) + .05 # connectivity matrix
directed <- FALSE
vBlocks <- 1:5 # number of classes

## Simulation of an Bernoulli non-directed SBM
mySBM <- simulateSBM(N, alpha, pi, directed)


### ==================
### MAR edge :
### ==================

## Sampling of the data : ##
samplingParameters <- .5
sampling <- "dyad"
sampledNetwork <-
  samplingSBM(
    mySBM$adjacencyMatrix,
    sampling,
    samplingParameters
  )

## Inference : ##
sbm_dyad <- missSBM(sampledNetwork$adjacencyMatrix,  seq.Q = Q, missingness = "mar.edge", cl.init = "spectral")

## Estimated parameters : ##
alpha_dyad <- getAlpha(getBestModel(sbm_dyad))
pi_dyad    <- getPi(getBestModel(sbm_dyad))
cl_dyad    <- getClusters(getBestModel(sbm_dyad))

aricode::ARI(getClusters(getBestModel(sbm_dyad)), mySBM$memberships)

### ==================
### MAR node :
### ==================

## Sampling of the data : ##
samplingParameters <- .5
sampling <- "node"
sampledNetwork <-
  samplingSBM(
    mySBM$adjacencyMatrix,
    sampling,
    samplingParameters
  )

## Inference : ##
sbm_node <- missSBM(sampledNetwork$adjacencyMatrix,  seq.Q = Q, missingness = "mar.node", cl.init = "spectral")

## Estimated parameters : ##
alpha_node <- getAlpha(getBestModel(sbm_node))
pi_node    <- getPi(getBestModel(sbm_node))
cl_node    <- getClusters(getBestModel(sbm_node))

aricode::ARI(getClusters(getBestModel(sbm_node)), mySBM$memberships)


### ===========================
### NMAR double standard :
### ===========================

## Sampling of the data : ##
samplingParameters <- c(.3,.7)
sampling <- "double_standard"
sampledNetwork <-
  samplingSBM(
    mySBM$adjacencyMatrix,
    sampling,
    samplingParameters
  )

## Inference : ##
sbm_twoStd <- func_missSBM.twoStd(sampledNetwork$adjacencyMatrix, Q)

## Estimated parameters : ##
alpha_twoStd <- getAlpha(getBestModel(sbm_twoStd))
pi_twoStd    <- getPi(getBestModel(sbm_twoStd))
cl_twoStd    <- getClusters(getBestModel(sbm_twoStd))

aricode::ARI(getClusters(getBestModel(sbm_twoStd)), mySBM$memberships)

### ==================
### NMAR block :
### ==================

## Sampling of the data : ##
samplingParameters <- c(0.75,0.5,0.1)
sampling <- "block"
sampledNetwork <-
  samplingSBM(
    mySBM$adjacencyMatrix,
    sampling,
    samplingParameters,
    mySBM$memberships
  )
## Inference : ##
sbm_block <- func_missSBM.class(sampledNetwork$adjacencyMatrix, Q, cl.init = "spectral")

## Taux d'échantillonnage : ##
cat("\n proportion of observed edges: ",sum(!is.na(sampledNetwork$adjacencyMatrix))/length(sampledNetwork$adjacencyMatrix))

## Estimated parameters : ##
alpha_block <- getAlpha(getBestModel(sbm_block))
pi_block    <- getPi(getBestModel(sbm_block))
cl_block    <- getClusters(getBestModel(sbm_block))

### ==================
### NMAR degree :
### ==================

## Sampling of the data : ##
samplingParameters <- c(-3.6,0.1)
sampling <- "degree"
sampledNetwork <-
  samplingSBM(
    mySBM$adjacencyMatrix,
    sampling,
    samplingParameters
  )

## Taux d'échantillonnage : ##
cat("\n proportion of observed edges: ",sum(!is.na(sampledNetwork$adjacencyMatrix))/length(sampledNetwork$adjacencyMatrix))


## Inference : ##
sbm_degree <- func_missSBM.degree(sampledNetwork$adjacencyMatrix, Q)

## Estimated parameters : ##
alpha_degree <- getAlpha(getBestModel(sbm_degree))
pi_degree    <- getPi(getBestModel(sbm_degree))
cl_degree    <- getClusters(getBestModel(sbm_degree))

aricode::ARI(getClusters(getBestModel(sbm_degree)), mySBM$memberships)

### ======================
### MCAR|X covariables :
### ======================

mySBM <- simulateSBM(N, alpha, pi, directed, covariates, covarParam)

## Sampling of the data : ##
samplingParameters <- c(0.1,-5)
sampling <- "covariates"
sampledNetwork <-
  samplingSBM(
    mySBM$adjacencyMatrix,
    sampling,
    samplingParameters
  )

## Inference : ##
sbm_covariable <- func_missSBM.degree(sampledNetwork$adjacencyMatrix, Q)

# # =====
#
# sbm_node <- missSBM(sampledNetwork$adjacencyMatrix,  seq.Q = Q, missingness = "mar.node", cl.init = "spectral")
#
# ## Estimated parameters : ##
# alpha_node <- getAlpha(getBestModel(sbm_node))
# pi_node    <- getPi(getBestModel(sbm_node))
# cl_node    <- getClusters(getBestModel(sbm_node))
#
# aricode::ARI(getClusters(getBestModel(sbm_node)), mySBM$memberships)
# distPI(pi, getPi(getModel(sbm_node, Q)))
# aricode::ARI(getClusters(getBestModel(sbm_block)), mySBM$memberships)
# distPI(pi, getPi(getModel(sbm_block, Q)))
