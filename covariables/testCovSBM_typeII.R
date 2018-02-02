rm(list = ls())
library(mclust)
source("~/Git/missSBM/covariables/drawCovSBM.R")
source("~/Git/missSBM/covariables/VEMCovSBM_typeII.R")
source("~/Git/missSBM/covariables/VEMCovSBM_typeI.R")

N <- 3
n <- 100
Q <- 3
pi <- diag(.45, Q) +.05

########################################################################################
### Test VEM with covariates type II :
########################################################################################


# Simulation d'un graphe
X     <- drawCov(N,n)
cov   <- phi(X)

# Simulation d'un SBM et d'une matrice d'échantillonnage
sbm <- drawCovSBM_typeII(N,n,Q,pi,X)
sampMat <- sampleCovSBM(sbm$Y,cov)

# Données
Y <- sbm$Y; Y[sampMat == 0] <- NA

# Taux d'échantillonnage
sum(sampMat)/n^2

# Inférence
infer <- func_missSBM.CovII(Y, 1:4, X)

# ICL
plot(infer@ICLs)

# Erreur d'estimation
cat(paste("\n ARI : ", adjustedRandIndex(getClusters(infer@models[[1]]), sbm$cl), "\n", "Error Pi : ", sum((getPi(infer@models[[1]])-pi)^2)/sum((getPi(infer@models[[1]]))^2)))


########################################################################################
### Test VEM with covariates type I :
########################################################################################

# Simulation d'un graphe
X     <- drawCov(N,n)
cov   <- phi(X)
alpha <- rep(1,Q)/Q
gamma <- 1/(1+exp(-pi))

# Simulation d'un SBM et d'une matrice d'échantillonnage
sbm <- drawCovSBM_typeI(N,n,Q,alpha,gamma,cov)
sampMat <- sampleCovSBM(sbm$Y,cov)

# Données
Y <- sbm$Y; Y[sampMat == 0] <- NA

# Taux d'échantillonnage
sum(sampMat)/n^2

# Inférence
infer <- func_missSBM.CovI(Y, 3, cov)

# ICL
plot(infer@ICLs)

# Erreur d'estimation
cat(paste("\n ARI : ", adjustedRandIndex(getClusters(infer@models[[1]]), sbm$cl), "\n", "Error Pi : ", sum((getPi(infer@models[[1]])-sbm$Z%*%pi%*%t(sbm$Z))^2)/sum((getPi(infer@models[[1]]))^2)))

