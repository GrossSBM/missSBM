
source("~/Git/missSBM/covariables/drawCovSBM_typeII.R")
source("~/Git/missSBM/covariables/VEM_covSBM_typeII.R")

### Simulation d'un graphe :
N <- 3
n <- 100
Q <- 2
pi <- diag(.45, Q) +.05

X     <- drawCov(N,n)
beta  <- drawBeta(N,Q)
alpha <- drawAlpha(n,Q,beta,X)

sbm <- drawCovSBM(N,n,Q,pi,X,beta)

### Inférence :

infer <- func_missSBM.CovII(sbm$Y, 2, X)
adjustedRandIndex(getClusters(infer@models[[1]]), sbm$cl)