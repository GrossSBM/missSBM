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
samplingParameters <- .5
sampling <- "dyad"
sampSBM <- samplingSBM(adjacencyMatrix, sampling, samplingParameters)
sampledAdjMatrix <- sampSBM$adjacencyMatrix
vBlocks <- 1:10
sbm <- inferSBM(sampledAdjMatrix, vBlocks, "dyad")

smoothingBackward <- function(models,vBlocks,sampledAdjMatrix,sampling) {
  for(i in rev(vBlocks[-1])){
    comb <- combn(i, 2, simplify = FALSE)
    for(j in 1:length(comb)){
      cl_fusion <- factor(apply(models[[i]]$fittedSBM$blocks,1, which.max))
      if(length(levels(cl_fusion)) == i){
        levels(cl_fusion)[which(levels(cl_fusion) == paste(comb[[j]][1]))] <- paste(comb[[j]][2])
        levels(cl_fusion) <- as.character(1:(i-1))
        clone <- inferSBM(sampledAdjMatrix, i-1, sampling, cl_fusion)
        clone$models[[1]]$doVEM()
        if(clone$models[[1]]$vICL < models[[i-1]]$vICL){
          models[[i-1]] <- clone$models[[1]]
          cat('+')
        }
      }
    }
  }
  models
}

smoothingForward <- function(models, vBlocks, sampledAdjMatrix,sampling) {
  for(i in vBlocks[-length(vBlocks)]){
    cl_split <- factor(apply(models[[i]]$fittedSBM$blocks,1,which.max))
    tab <- ceiling(tabulate(cl_split)/2)
    levels(cl_split) <- c(levels(cl_split),as.character(i+1))
    for (j in 1:i) {
      if(length(levels(cl_split))-1 == i){
        cl <- as.numeric(cl_split); cl[which(cl==j)][1:tab[j]] <- i+1
        clone <- inferSBM(sampledAdjMatrix, i+1, sampling, cl)
        clone$models[[1]]$doVEM()
        if(clone$models[[1]]$vICL < models[[i+1]]$vICL){
          models[[i+1]] <- clone$models[[1]]
          cat('+')
        }
      }
    }
  }
  models
}

smoothingForward_2 <- function(models, vBlocks, sampledAdjMatrix,sampling) {
  for(i in vBlocks[-length(vBlocks)]){
    cl_split <- factor(apply(models[[i]]$fittedSBM$blocks,1,which.max))
    levels(cl_split) <- c(levels(cl_split),as.character(i+1))
    for (j in 1:i) {
      if(length(levels(cl_split))-1 == i){
        # browser()
        cl <- as.numeric(cl_split); indices <- which(cl==j)
        cut <- as.numeric(init_spectral(sampledAdjMatrix[indices, indices],2))
        cl[which(cl==j)][which(cut==1)] <- j; cl[which(cl==j)][which(cut==2)] <- i+1
        clone <- inferSBM(sampledAdjMatrix, i+1, sampling, cl)
        clone$models[[1]]$doVEM()
        if(clone$models[[1]]$vICL < models[[i+1]]$vICL){
          models[[i+1]] <- clone$models[[1]]
          cat('+')
        }
      }
    }
  }
  models
}

smoothingForBackWard <- function(models, vBlocks, sampledAdjMatrix,sampling, nIter=2){
  out <- models
  for (i in 1:nIter) {
    out <- smoothingBackward(out, vBlocks, sampledAdjMatrix,sampling)
    out <- smoothingForward(out, vBlocks, sampledAdjMatrix,sampling)
  }
  out
}

out <- smoothingBackward(sbm$models, vBlocks, sampledAdjMatrix,sampling)
out2 <- smoothingForward(sbm$models, vBlocks, sampledAdjMatrix,sampling)
out3 <- smoothingForBackWard(sbm$models, vBlocks, sampledAdjMatrix,sampling)
out4 <- smoothingForward_2(sbm$models, vBlocks, sampledAdjMatrix,sampling)

# par(mfrow=c(2,2))
vICL <- sapply(sbm$models, function(model) model$vICL); plot(vICL)
vICL_smoothed <- sapply(out, function(model) model$vICL); plot(vICL_smoothed)
vICL_smoothed2 <- sapply(out2, function(model) model$vICL); plot(vICL_smoothed2)
vICL_smoothed3 <- sapply(out3, function(model) model$vICL); plot(vICL_smoothed3)



