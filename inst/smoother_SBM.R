rm(list=ls())
library(missSBM)
library(pbmcapply)

## Model :
n <- 200
Q <- 4
pi <- diag(0.45,Q) + .05
alpha <- rep(1,Q)/Q
family <- "Bernoulli"
directed <- FALSE
mySBM <- simulateSBM(n, alpha, pi, directed)
adjacencyMatrix <- mySBM$adjacencyMatrix
samplingParameters <- .3
sampling <- "dyad"
sampSBM <- samplingSBM(adjacencyMatrix, sampling, samplingParameters)
sampledAdjMatrix <- sampSBM$adjacencyMatrix
vBlocks <- 1:10
sbm <- inferSBM(sampledAdjMatrix, vBlocks, "dyad")

smoothingBackward <- function(models, vBlocks, sampledAdjMatrix, sampling, mc.cores = 4) {
  sampledNet <- sampledNetwork$new(sampledAdjMatrix)
  for (i in rev(vBlocks[-1])) {
    cl0 <- factor(models[[i]]$fittedSBM$memberships)
    if (nlevels(cl0) == i) {
      candidates <- pbmclapply(combn(i, 2, simplify = FALSE), function(couple) {
        cl_fusion <- cl0
        levels(cl_fusion)[which(levels(cl_fusion) == paste(couple[1]))] <- paste(couple[2])
        levels(cl_fusion) <- as.character(1:(i - 1))
        model <- missingSBM_fit$new(sampledNet, i - 1, sampling, cl_fusion)
        model$doVEM(trace = FALSE)
        model
      }, mc.cores = mc.cores)
      vICLs <- sapply(candidates, function(candidate) candidate$vICL)
      best_one <- candidates[[which.min(vICLs)]]
      if (best_one$vICL < models[[i - 1]]$vICL)
        models[[i - 1]] <- best_one
    }
  }
  models
}

smoothingForward <- function(models, vBlocks, sampledAdjMatrix, sampling, mc.cores = 4) {
  sampledNet <- sampledNetwork$new(sampledAdjMatrix)
  for(i in vBlocks[-length(vBlocks)]){
    cl_split <- factor(models[[i]]$fittedSBM$memberships)
    tab <- ceiling(tabulate(cl_split)/2)
    levels(cl_split) <- c(levels(cl_split), as.character(i+1))
    if (nlevels(cl_split) - 1 == i) {
      candidates <- pbmclapply(1:i, function(j) {
        cl <- as.numeric(cl_split); cl[which(cl==j)][1:tab[j]] <- i+1
        model <- missingSBM_fit$new(sampledNet, i + 1, sampling, cl)
        model$doVEM(trace = FALSE)
        model
      }, mc.cores = mc.cores)
      vICLs <- sapply(candidates, function(candidate) candidate$vICL)
      best_one <- candidates[[which.min(vICLs)]]
      if (best_one$vICL < models[[i + 1]]$vICL)
        models[[i + 1]] <- best_one
    }
  }
  models
}

smoothingForward_2 <- function(models, vBlocks, sampledAdjMatrix,sampling) {
  sampledNet <- sampledNetwork$new(sampledAdjMatrix)
  for(i in vBlocks[-length(vBlocks)]){
    cl_split <- factor(apply(models[[i]]$fittedSBM$blocks,1,which.max))
    levels(cl_split) <- c(levels(cl_split),as.character(i+1))
    for (j in 1:i) {
      if(length(levels(cl_split))-1 == i){
        # browser()
        cl <- as.numeric(cl_split); indices <- which(cl==j)
        cut <- as.numeric(init_spectral(sampledAdjMatrix[indices, indices],2))
        cl[which(cl==j)][which(cut==1)] <- j; cl[which(cl==j)][which(cut==2)] <- i+1
        clone <- missingSBM_fit$new(sampledNet, i + 1, sampling, cl)
        clone$doVEM(trace = FALSE)
        if(clone$vICL < models[[i+1]]$vICL){
          models[[i+1]] <- clone$models[[1]]
          cat('+')
        }
      }
    }
  }
  models
}

smoothingForBackWard <- function(models, vBlocks, sampledAdjMatrix, sampling, nIter = 2){
  out <- models
  for (i in 1:nIter) {
    out <- smoothingBackward(out, vBlocks, sampledAdjMatrix,sampling)
    out <- smoothingForward(out, vBlocks, sampledAdjMatrix,sampling)
  }
  out
}

smoothed_fwrd <- smoothingForward(sbm$models, vBlocks, sampledAdjMatrix, sampling)
smoothed_back <- smoothingBackward(sbm$models, vBlocks, sampledAdjMatrix, sampling)
smoothed_fb   <- smoothingForBackWard(sbm$models, vBlocks, sampledAdjMatrix, sampling)
## out4 <- smoothingForward_2(sbm$models, vBlocks, sampledAdjMatrix,sampling)

vICL <- sapply(sbm$models, function(model) model$vICL); plot(vICL, type = "l")
vICL_back <- sapply(smoothed_back, function(model) model$vICL); lines(vICL_back, col="blue")
vICL_fwrd <- sapply(smoothed_fwrd, function(model) model$vICL); lines(vICL_fwrd, col="red")
vICL_fb   <- sapply(smoothed_fb  , function(model) model$vICL); lines(vICL_fb  , col="green")



