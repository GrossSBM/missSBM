#' @export
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
      if(is.na(models[[i - 1]]$vICL)){
        models[[i - 1]] <- best_one
      } else if(best_one$vICL < models[[i - 1]]$vICL) {
        models[[i - 1]] <- best_one
      }
    }
  }
  models
}

#' @export
smoothingForward_half <- function(models, vBlocks, sampledAdjMatrix, sampling, mc.cores = 4) {
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

#' @export
smoothingForward_SpCl <- function(models, vBlocks, sampledAdjMatrix, sampling, mc.cores = 4) {
  sampledNet <- sampledNetwork$new(sampledAdjMatrix)
  for(i in vBlocks[-length(vBlocks)]){
    cl_split <- factor(models[[i]]$fittedSBM$memberships)
    levels(cl_split) <- c(levels(cl_split), as.character(i+1))
    if (nlevels(cl_split) - 1 == i) {
      candidates <- pbmclapply(1:i, function(j) {
        # if(i==6) browser()
        cl <- as.numeric(cl_split); indices <- which(cl==j)
        if(length(cl[indices]) > 10){
          cut <- as.numeric(init_spectral(sampledAdjMatrix[indices, indices],2))
          cl[which(cl==j)][which(cut==1)] <- j; cl[which(cl==j)][which(cut==2)] <- i + 1
          model <- missingSBM_fit$new(sampledNet, i + 1, sampling, cl)
          model$doVEM(trace = FALSE)
          model
        } else {
          models[[i + 1]]
        }
      }, mc.cores = mc.cores)
      vICLs <- sapply(candidates, function(candidate) candidate$vICL)
      best_one <- candidates[[which.min(vICLs)]]
      if(is.na(models[[i + 1]]$vICL)){
        models[[i + 1]] <- best_one
      } else if(best_one$vICL < models[[i + 1]]$vICL) {
        models[[i + 1]] <- best_one
      }
    }
  }

  models
}

#' @export
smoothingForBackWard_half <- function(models, vBlocks, sampledAdjMatrix, sampling, nIter = 2){
  out <- models
  for (i in 1:nIter) {
    out <- smoothingBackward(out, vBlocks, sampledAdjMatrix,sampling)
    out <- smoothingForward_half(out, vBlocks, sampledAdjMatrix,sampling)
  }
  out
}

#' @export
smoothingForBackWard_SpCl <- function(models, vBlocks, sampledAdjMatrix, sampling, nIter = 2){
  out <- models
  for (i in 1:nIter) {
    out <- smoothingBackward(out, vBlocks, sampledAdjMatrix,sampling)
    out <- smoothingForward_SpCl(out, vBlocks, sampledAdjMatrix,sampling)
  }
  out
}
