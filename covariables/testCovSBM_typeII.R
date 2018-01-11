
library(mclust)
source("~/Git/missSBM/covariables/drawCovSBM_typeII.R")
source("~/Git/missSBM/covariables/VEM_covSBM_typeII.R")

### Simulation d'un graphe :
N <- 3
n <- 100
Q <- 2
pi <- diag(.45, Q) +.05

X     <- rbind(1,drawCov(N,n))
beta  <- drawBeta(N+1,Q)
alpha <- drawAlpha(n,Q,beta,X)

sbm <- drawCovSBM(N+1,n,Q,pi,X,beta,TRUE)
sampMat <- sampleCovSBM_Node(sbm$Y,X, directed=T)

Y <- sbm$Y; Y[sampMat == 1] <- NA

### Inférence :

infer <- func_missSBM.CovII(Y, 1:4, X, directed=T)
# getBeta(infer@models[[1]])
plot(infer@ICLs)
adjustedRandIndex(getClusters(infer@models[[2]]), sbm$cl)
sum(sampMat)/10000
