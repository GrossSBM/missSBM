library(missSBM)


### Reprise code timothee 150219
## Paramètres SBM
# Q <- 3
n <- 100
epsilon <- .1
# alpha   <- rep(1,Q)/Q
# pi      <- diag(epsilon,Q) + .05
directed <- FALSE

## Paramètres Covariabes
N <- 1
beta  <- 1

## Paramètres de simulation
nbrSimu <- 1
res1    <- data.frame()
seq.p <- seq(.4,1,length=5)


type_graph= c("affiliation")

alpha <- switch(type_graph,
                "affiliation" = c(1/3,1/3,1/3),
                "bipartite"   = c(1/4,1/4,1/4,1/4),
                "star"        = c(.15,.35,.15,.35))
Q <- length(alpha)


epsilon = .1
delta = .5
  pi <- switch(type_graph,
               "affiliation" = matrix(c(delta,epsilon,epsilon,epsilon,delta,epsilon,epsilon,epsilon,delta),3,3),
               "bipartite"   = matrix(c(epsilon,1-epsilon,epsilon,epsilon,1-epsilon,epsilon,epsilon,epsilon,epsilon,epsilon, epsilon,1-epsilon,epsilon,epsilon,1-epsilon,epsilon),4,4),
               "star"        = matrix(c(1-epsilon,1-epsilon,epsilon,0,1-epsilon,epsilon,0,0,epsilon,0,1-epsilon,1-epsilon,0,0,1-epsilon,epsilon),4,4))

  gamma <- log(pi)-log(1-pi)


library(sbm)
# A$memberships
# A$networkData

#sampPar = matrix(.2,3,3) + rbind(c(.4,.4,.4),c(.4,.2,.2),c(.4,.2,0))
sampPar = matrix(rbeta(9,1,2),3,3)
curve(dbeta(x,1,2))

sampPar = matrix(.5,3,3)+diag(3)*.3

library(parallel)
RES = mclapply(1:50,function(i){


  A=sbm::sampleSimpleSBM(nbNodes = n,blockProp = alpha,connectParam = list(mean=pi),model = "bernoulli",directed = F)

Aobs = missSBM::observeNetwork(A$networkData,sampling="block-dyad",clusters = A$memberships,parameters = sampPar)
# Aobs
# sum(is.na(Aobs))/prod(dim(Aobs))
# a comparer a
#3*.33^2 *.9 + 6*.33^2 *.1^2

#class(A$networkData)

ResNMAR = missSBM::estimateMissSBM(Aobs,1:10,sampling="block-dyad")
missSBM::smooth(ResNMAR,"both")
#plot(ResNMAR$ICL)
# ResNMAR$bestModel$fittedSBM$blockProp
# ResNMAR$bestModel$fittedSBM$connectParam
# ResNMAR$bestModel$fittedSampling$parameters


a=aricode::ARI(ResNMAR$models[[3]]$fittedSBM$memberships,A$memberships)
b=sqrt(sum((ResNMAR$models[[3]]$fittedSBM$connectParam$mean-pi)^2))


ResMAR = missSBM::estimateMissSBM(Aobs,1:10,sampling="dyad")
missSBM::smooth(ResMAR,"both")
#plot(ResMAR$ICL)
# ResMAR$bestModel$fittedSBM$blockProp
# ResMAR$bestModel$fittedSBM$connectParam
# ResMAR$bestModel$fittedSampling$parameters
c=aricode::ARI(ResMAR$models[[3]]$fittedSBM$memberships,A$memberships)
d=sqrt(sum((ResMAR$models[[3]]$fittedSBM$connectParam$mean-pi)^2))
return(c(a,c,b,d))
},mc.cores = 10)
RES2 = do.call("rbind",RES)
RES2 = as.data.frame(RES2)
names(RES2) = rep(c("block-dyad","dyad"),2)
boxplot(RES2[,1:2])
title("ARI, un ami qui vous veut du bien")
boxplot(RES2[,1]-RES2[,2])
abline(h=0)
boxplot(RES2[,3:4])
boxplot(RES2[,3]-RES2[,4],ylim=c(-.1,.1)/8)
abline(h=0)
#save.image(file="~/Documents/Recherche/Timothee/simublockdyad110321.Rdata")

