rm(list = ls())
library(mclust)
source("~/Git/missSBM/covariables/drawCovSBM.R")
source("~/Git/missSBM/covariables/VEMCovSBM_typeII.R")
source("~/Git/missSBM/covariables/VEMCovSBM_typeI.R")

N <- 3
n <- 200
Q <- 3
pi <- diag(.45, Q) +.05

########################################################################################
### Test VEM with covariates type II :
########################################################################################


### Simulation d'un graphe :

X     <- drawCov(N,n)
X <- phi(X)
beta  <- drawTheta(N)
# alpha <- drawAlpha(n,Q,beta,X)
alpha <- rep(1,3)/3

sbm <- drawCovSBM_typeI(N,n,Q,alpha,pi,X,beta)
sampMat <- sampleCovSBMII(sbm$Y,X)

Y <- sbm$Y; Y[sampMat == 0] <- NA

### Inférence :

infer <- func_missSBM.CovI(Y, 3, X)
# getBeta(infer@models[[1]])
plot(infer@ICLs)
adjustedRandIndex(getClusters(infer@models[[3]]), sbm$cl)
sum(sampMat)/n^2


########################################################################################
### Test VEM with covariates type II :
########################################################################################


