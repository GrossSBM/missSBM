#' @export
smoothingBackward <- function(models, vBlocks, sampledNet, sampling, mc.cores, control) {
  cat("   Going backward ")
  for (i in rev(vBlocks[-1])) {
    cat('+')
    cl0 <- factor(models[[i]]$fittedSBM$memberships)
    if (nlevels(cl0) == i) {
      candidates <- mclapply(combn(i, 2, simplify = FALSE), function(couple) {
        cl_fusion <- cl0
        levels(cl_fusion)[which(levels(cl_fusion) == paste(couple[1]))] <- paste(couple[2])
        levels(cl_fusion) <- as.character(1:(i - 1))
        model <- missingSBM_fit$new(sampledNet, i - 1, sampling, cl_fusion)
        model$doVEM(control)
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
  cat("\r                                                                                                        \r")
  models
}

#' @export
smoothingForward_half <- function(models, vBlocks, sampledNet, sampling, mc.cores, control) {
  cat("   Going forward ")
  for(i in vBlocks[-length(vBlocks)]){
    cat('+')
    cl_split <- factor(models[[i]]$fittedSBM$memberships)
    tab <- ceiling(tabulate(cl_split)/2)
    levels(cl_split) <- c(levels(cl_split), as.character(i+1))
    if (nlevels(cl_split) - 1 == i) {
      candidates <- mclapply(1:i, function(j) {
        cl <- as.numeric(cl_split); cl[which(cl==j)][1:tab[j]] <- i+1
        model <- missingSBM_fit$new(sampledNet, i + 1, sampling, cl)
        model$doVEM(control)
        model
      }, mc.cores = mc.cores)
      vICLs <- sapply(candidates, function(candidate) candidate$vICL)
      best_one <- candidates[[which.min(vICLs)]]
      if (best_one$vICL < models[[i + 1]]$vICL)
        models[[i + 1]] <- best_one
    }
  }
  cat("\r                                                                                                        \r")
  models
}

#' @export
smoothingForward_SpCl <- function(models, vBlocks, sampledNet, sampling, mc.cores, control) {
  cat("   Going forward ")
  for(i in vBlocks[-length(vBlocks)]){
    cat("+")
    cl_split <- factor(models[[i]]$fittedSBM$memberships)
    levels(cl_split) <- c(levels(cl_split), as.character(i+1))
    if (nlevels(cl_split) - 1 == i) {
      candidates <- mclapply(1:i, function(j) {
        cl <- as.numeric(cl_split); indices <- which(cl == j)
        if(length(cl[indices]) > 10){
          cut <- as.numeric(init_spectral(sampledNet$adjacencyMatrix[indices, indices],2))
          cl[which(cl==j)][which(cut==1)] <- j; cl[which(cl==j)][which(cut==2)] <- i + 1
          model <- missingSBM_fit$new(sampledNet, i + 1, sampling, cl)
          model$doVEM(control)
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
  cat("\r                                                                                                        \r")
  models
}

#' @export
smoothingForBackWard_half <- function(models, vBlocks, sampledNet, sampling, mc.cores, control){
  out <- models
  out <- smoothingBackward(out, vBlocks, sampledNet, sampling, mc.cores, control)
  out <- smoothingForward_half(out, vBlocks, sampledNet, sampling, mc.cores, control)
  out
}

#' @export
smoothingForBackWard_SpCl <- function(models, vBlocks, sampledNet, sampling, mc.cores, control){
  out <- models
  out <- smoothingBackward(out, vBlocks, sampledNet, sampling, mc.cores, control)
  out <- smoothingForward_SpCl(out, vBlocks, sampledNet, sampling, mc.cores, control)
  out
}
