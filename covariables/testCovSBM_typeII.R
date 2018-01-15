rm(list = ls())
library(mclust)
source("~/Git/missSBM/covariables/drawCovSBM_typeII.R")
source("~/Git/missSBM/covariables/VEM_covSBM_typeII.R")

N <- 3
n <- 200
Q <- 3
pi <- diag(.45, Q) +.05

########################################################################################
### Test VEM with covariates type II :
########################################################################################


### Simulation d'un graphe :

X     <- drawCov(N,n)
beta  <- drawBeta(N,Q)
alpha <- drawAlpha(n,Q,beta,X)

sbm <- drawCovSBM_typeII(N,n,Q,pi,X,beta)
sampMat <- sampleCovSBM(sbm$Y,X)

Y <- sbm$Y; Y[sampMat == 1] <- NA

### Inférence :

infer <- func_missSBM.CovII(Y, 1:4, X)
# getBeta(infer@models[[1]])
plot(infer@ICLs)
adjustedRandIndex(getClusters(infer@models[[3]]), sbm$cl)
sum(sampMat)/n^2


########################################################################################
### Test VEM with covariates type II :
########################################################################################


