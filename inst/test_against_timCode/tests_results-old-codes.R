### ====================== ###
### Tests anciens codes :  ###
### ====================== ###
rm(list=ls())
library(missSBM)

if(Sys.info()["user"] == "tabouy") {
  setwd("~/GitLab/svn_oldies/Code/code_these/functions/")
} else if (Sys.info()["user"] == "tim") {
  setwd("~/Documents/svn_oldies/Code/code_these/functions/")
}
source("func_missSBM.R")
source("func_missSBM.twoStd.R")
source("func_missSBM.class.R")
source("func_missSBM.degree.R")
source("func_sampling.R")

source("~/GitLab/covariables-SBM/code/drawCovSBM.R")
source("~/GitLab/covariables-SBM/code/VEMCovSBM_model2.R")

## A SBM model : ##
N <- 100
Q <- 3
alpha <- c(0.25,0.5,0.25) # rep(1,Q)/Q     # mixture parameter
pi <- diag(.45,Q) + .05 # connectivity matrix
directed <- FALSE
vBlocks <- 1:5 # number of classes

## Simulation of an Bernoulli non-directed SBM
mySBM <- simulateSBM(N, alpha, pi, directed)

l.sbm <- list(Y=mySBM$adjacencyMatrix,n=N,Q=Q,alpha=alpha,pi=pi,directed=directed,cl=mySBM$memberships)

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

# aricode::ARI(getClusters(getBestModel(sbm_dyad)), mySBM$memberships)

l.dyad <- list(Y=sampledNetwork$adjacencyMatrix,alpha=alpha_dyad,pi=pi_dyad,cl=cl_dyad)


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

# aricode::ARI(getClusters(getBestModel(sbm_node)), mySBM$memberships)

l.node <- list(Y=sampledNetwork$adjacencyMatrix,alpha=alpha_node,pi=pi_node,cl=cl_node)

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
psi_twoStd   <- getPsi(getBestModel(sbm_twoStd))

# aricode::ARI(getClusters(getBestModel(sbm_twoStd)), mySBM$memberships)

l.twoStd <- list(Y=sampledNetwork$adjacencyMatrix,alpha=alpha_twoStd,pi=pi_twoStd,cl=cl_twoStd,psi_vrai=samplingParameters,psi_chap=psi_twoStd)

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
    clusters = mySBM$memberships
  )
## Inference : ##
sbm_block <- func_missSBM.class(sampledNetwork$adjacencyMatrix, Q, cl.init = "spectral")

## Taux d'échantillonnage : ##
cat("\n proportion of observed edges: ",sum(!is.na(sampledNetwork$adjacencyMatrix))/length(sampledNetwork$adjacencyMatrix))

## Estimated parameters : ##
alpha_block <- getAlpha(getBestModel(sbm_block))
pi_block    <- getPi(getBestModel(sbm_block))
cl_block    <- getClusters(getBestModel(sbm_block))
psi_block   <- getPsi(getBestModel(sbm_block))

l.block <- list(Y=sampledNetwork$adjacencyMatrix,alpha=alpha_block,pi=pi_block,cl=cl_block,psi_vrai=samplingParameters,psi_chap=psi_block)

### ==================
### NMAR degree :
### ==================

## Sampling of the data : ##
samplingParameters <- c(-3.6,0.12)
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
alpha_degree <- getAlpha(getModel(sbm_degree,Q))
pi_degree    <- getPi(getModel(sbm_degree,Q))
cl_degree    <- getClusters(getModel(sbm_degree,Q))
psi_degree   <- getPsi(getModel(sbm_degree,Q))

# aricode::ARI(getClusters(getBestModel(sbm_degree)), mySBM$memberships)

l.degree <- list(Y=sampledNetwork$adjacencyMatrix,alpha=alpha_degree,pi=pi_degree,cl=cl_degree,psi_vrai=samplingParameters,psi_chap=psi_degree)

### ==========================================
### Modèle SBM avec Covariables :
### ==========================================

## Paramètres d'échantillonnage et du SBM : ##
n <- 100
N <- 1
beta  <- 1

epsilon <- .25
pi      <- diag(epsilon,Q) + .05
gamma   <- log(pi)-log(1-pi)
directed <- FALSE
## SBM : ##
cov  <- drawCov(N,n)
sbm  <- drawCovSBM_mod2(N,n,Q,alpha,gamma,X=phi(cov),beta,directed)

l.sbmCov <- list(Y=sbm$Y,n=n,N=N,Q=Q,alpha=alpha, beta=beta,gamma=gamma,directed=directed, cl=sbm$cl, cov=phi(cov))


### ==========================================
### MCAR|X covariables  sur les dyades :
### ==========================================

## Sampling of the data : ##
sampMat <- sampleCovSBM(N,sbm$Y,phi(cov),beta,-.2, Node=FALSE)
## Matrice d'adjacence avec des NA : ##
Y <- sbm$Y ; Y[sampMat == 0] <- NA

## Taux d'échantillonnage : ##
print(sum(sampMat)/n^2)

## Inférence MCAR covariable
infer <- func_missSBM.CovMod2(Y, Q, phi(cov))

## Estimated parameters : ##
alpha_covDyad <- getAlpha(getBestModel(infer))
pi_covDyad    <- getPi(getBestModel(infer))
cl_covDyad    <- getClusters(getBestModel(infer))
gamma_covDyad    <- getGamma(getBestModel(infer))


l.covDyad <- list(Y=sbm$Y,alpha=alpha_covDyad,pi=pi_covDyad,gamma_covDyad=gamma_covDyad,cl=cl_covDyad)


### ==========================================
### MCAR|X covariables  sur les noeuds :
### ==========================================

## Sampling of the data : ##
sampMat <- sampleCovSBM(N,sbm$Y,cov,5*beta,-3.5, Node = TRUE)
## Matrice d'adjacence avec des NA : ##
Y <- sbm$Y ; Y[sampMat == 0] <- NA

## Taux d'échantillonnage : ##
print(sum(sampMat)/n^2)

## Inférence MCAR covariable
infer <- func_missSBM.CovMod2(Y, Q, phi(cov))

## Estimated parameters : ##
alpha_covNode <- getAlpha(getBestModel(infer))
pi_covNode    <- getPi(getBestModel(infer))
cl_covNode    <- getClusters(getBestModel(infer))
gamma_covNode    <- getGamma(getBestModel(infer))


l.covNode <- list(Y=sbm$Y,alpha=alpha_covNode,pi=pi_covNode,gamma_covDyad=gamma_covNode,cl=cl_covNode)


### ==========================
### Sauvegarde des données :
### ==========================


save(l.sbm, l.dyad, l.node, l.twoStd, l.block, l.degree, l.sbmCov, l.covDyad, l.covNode, file = "~/Git/missSBM/inst/test_against_timCode/codesTim.RData")
