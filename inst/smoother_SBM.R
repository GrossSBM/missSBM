rm(list=ls())
library(missSBM)

## Model :
n <- 200
Q <- 4
pi <- diag(0.45,Q)+.05
alpha <- rep(1,Q)/Q
family <- "Bernoulli"
directed <- FALSE
mySBM <- simulateSBM(n, alpha, pi, directed)
adjacencyMatrix <- mySBM$adjacencyMatrix
samplingParameters <- 1
sampling <- "dyad"
sampSBM <- samplingSBM(adjacencyMatrix, sampling, samplingParameters)
sampledAdjMatrix <- sampSBM$adjacencyMatrix
vBlocks <- 1:6
sbm <- inferSBM(sampledAdjMatrix, vBlocks, "double_standard")


smoothingBackward <- function(models,vBlocks) {
  for(i in rev(vBlocks[-1])){
    comb <- combn(i, 2, simplify = FALSE)
    for(j in 1:length(comb)){
      cl_fusion <- factor(models[[i]]$fittedSBM$blocks)
      if(length(levels(cl_fusion)) == i){
        levels(cl_fusion)[which(levels(cl_fusion) == paste(comb[[j]][1]))] <- paste(comb[[j]][2])
        levels(cl_fusion) <- as.character(1:(i-1))
        clone             <- missingSBM_fit$new(sampledNet, i-1, sampling, cl_fusion)
        clone$doVEM()
        if(clone$vICL < models[[i-1]]$vICL){
          models[[i-1]] <- clone
          cat('+')
        }
      }
    }
  }
}

smoothingForward <- function(models) {
  for(i in vBlocks[-length(vBlocks)]){
    cl_split <- factor(models[[i]]$fittedSBM$blocks)
    tab <- ceiling(tabulate(cl_split)/2)
    levels(cl_split) <- c(levels(cl_split),as.character(i+1))
    for (j in 1:i) {
      if(length(levels(cl_split))-1 == i){
        cl    <- cl_split; cl[which(cl==j)][1:tab[j]] <- i+1
        clone <- missingSBM_fit$new(sampledNet, i+1, sampling, cl)
        clone$doVEM()
        if(clone$vICL < models[[i+1]]$vICL){
          models[[i+1]] <- clone
          cat('+')
        }
      }
    }
  }
}
