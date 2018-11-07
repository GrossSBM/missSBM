rm(list=ls())
set.seed(12345)

library(missSBM)
library(aricode)

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
pi <- matrix(pr,2,2); pi[1,1] <- 3*pr
diri <- FALSE

# sampling parameters
a    <- 1/2; rho  <- .7
parameters <- c(rho*a,rho)

alpha1 <- .15
alpha  <- c(alpha1,alpha1*(3*a*rho-2*a^2*rho^2-rho+a*rho^2)/(-a*rho+a*rho^2+rho-rho^2))


## SBML parameters
Net        <- simulateSBM(N, alpha , pi, diri)
cl_star    <- as.vector(Net$blocks %*% 1:2)
Netsamp    <- samplingSBM(Net$adjacencyMatrix, "block", parameters, clusters = cl_star)
sampledNet <- Netsamp$adjacencyMatrix

# ancien code -------------------------------------------------------------
old <- func_missSBM.class(sampledNet, seq.Q = 2,  cl.init = "CAH")@models[[1]]

# package -----------------------------------------------------------------
control <- list(threshold = 1e-3, maxIter = 200, fixPointIter = 5, trace = TRUE)
new <- missingSBM_fit$new(Netsamp, 2, "block", clusterInit = "hierarchical")
optim_new <- new$doVEM(control)


plot(-optim_new$objective, type = "l", log = "y")
plot(-old@crit, type = "l", log = "y")

psi_new <- new$fittedSampling$parameters
psi_old <- old@psi[[1]]

pi_new <- new$fittedSBM$connectParam
pi_old <- old@theta[[1]]$pi

res <- data.frame(
  version  = c("new", "old"), sampling  = c("block", "block"),
  ARI = c(ARI(cl_star, new$fittedSBM$memberships), ARI(cl_star, getClusters(old))),
  err_pi  = c(frobenius(pi_new - matrix(rev(pi), 2, 2)), frobenius(pi_old - matrix(rev(pi), 2, 2))),
  err_psi = c(frobenius(rev(psi_new) - parameters), frobenius(rev(psi_old) - parameters))
)
print(res)
