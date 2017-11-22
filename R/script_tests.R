
##########################################################################
############################# Script de tests ############################
##########################################################################
rm(list=ls())
library(R6)

# Ancien :
setwd("~/SVN/Code/code_these/functions/")
source("func_missSBM.R")
source("func_missSBM.twoStd.R")
source("func_missSBM.degree.R")
source("func_missSBM.class.R")


# Nouveau :
setwd("~/Git/missSBM/R")
source("SBM-R6Class.R")
source("sampling-R6Class.R")
source("sampledNetwork-R6Class.R")
source("SBM_VEMfit-R6Class.R")
source("SBM_collection-R6Class.R")


####################
### Test général ###
####################


# double standard : check
# MAR edge        : check
# MAR node        : check
# star degree     : check
# class           : check

#### SBM : ####
mySBM <- SBM_BernoulliUndirected.fit$new(400, rep(1, 5)/5, diag(.45,5)+.05)

#### Sampled SBM : ####
mySampledSBM   <- sampling_starDegree$new(400, c(-2.5, 0.05), FALSE)
SBMdata        <- mySBM$rSBM()
Y              <- SBMdata$adjacencyMatrix

#### Pour le class sampling : ####
# Z              <- SBMdata$blocks
# Znum           <- Z %*% c(1:2)
# sample         <- mySampledSBM$rSampling(Y, Z)

#### Pour les autres : ####
sample         <- mySampledSBM$rSampling(Y)
sample
### Sampled SBM (MAR) :
# mySampledSBM  <- sampling_randomNodesMAR$new(100, 0.5, FALSE)
# Y             <- mySBM$rSBM()$adjacencyMatrix
# sample        <- mySampledSBM$rSampling(Y)

#### VEM  nouveau : ####
sbm <- SBM_collection$new(sample$adjacencyMatrix, 1:7, "starDegree", "Bernoulli", FALSE)

plot(sbm$vICLs)
sbm$getBestModel()

# ## Smoothing : ##
sbm$smoothingBackward()

#### VEM  ancien : ####
# res.twoStd   <- func_missSBM.twoStd(sample$adjacencyMatrix, 1:10)
# res.class    <- func_missSBM.class(sample$adjacencyMatrix, 1:5)
res.degree    <- func_missSBM.degree(sample$adjacencyMatrix, 1:10)
# res.mar      <- missSBM(sample$adjacencyMatrix, 2, missingness  = "class")
plot(res.degree@ICLs)

