#' @importFrom utils combn
#' @importFrom parallel mclapply
smoothingBackward <- function(models, vBlocks, sampledNet, sampling, split_fn, mc.cores, iter_both, control) {
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

smoothingForward <- function(models, vBlocks, sampledNet, sampling, split_fn, mc.cores, iter_both, control) {
  cat("   Going forward ")
  for(i in vBlocks[-length(vBlocks)]){
    cat("+")
    cl_split <- factor(models[[i]]$fittedSBM$memberships)
    levels(cl_split) <- c(levels(cl_split), as.character(i+1))
    if (nlevels(cl_split) - 1 == i) {
      candidates <- mclapply(1:i, function(j) {
        cl <- as.numeric(cl_split); indices <- which(cl == j)
        if(length(cl[indices]) > 10){
          # sampledNetBis <- sampledNet$adjacencyMatrix
          # sampledNetBis[is.na(sampledNetBis)] <- mean(sampledNetBis[!is.na(sampledNetBis)])
          cut <- as.numeric(split_fn(sampledNet$adjacencyMatrix[indices, indices],2))
          # cut <- c(rep(1, length = floor(length(indices)/2)), rep(2, length = ceiling(length(indices)/2)))
          cl[which(cl==j)][which(cut==1)] <- j; cl[which(cl==j)][which(cut==2)] <- i + 1
          model <- missingSBM_fit$new(sampledNet, i + 1, sampling, cl)
          model$doVEM(control)
          model
        } else {
          models[[i + 1]]
        }
      }, mc.cores = 1)
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

smoothingForBackWard <- function(models, vBlocks, sampledNet, sampling, split_fn, mc.cores, iter_both, control){
  out <- models
  for (i in 1: iter_both) {
    out <- smoothingForward(out, vBlocks, sampledNet, sampling, split_fn, mc.cores, iter_both, control)
    out <- smoothingBackward(out, vBlocks, sampledNet, sampling, split_fn, mc.cores, iter_both, control)
  }
  out
}

SpectralClustering_NAisMean <- function (A,Q) {
  if(Q > 1){
    ## basic handling of missing values
    if (anyNA(A)) A[is.na(A)] <- mean(A[!is.na(A)])

    ## handling lonely souls
    cl.final <- rep(NA, ncol(A))
    unconnected <- which(rowSums(A) == 0)
    connected <- setdiff(1:ncol(A), unconnected)

    A <- A[connected,connected]

    ## Normalized Laplacian
    D <- colSums(A)
    L <- diag(rep(1,ncol(A))) -
      diag(D^(-1/2)) %*% A %*% diag(D^(-1/2))

    ## Absolute eigenvalue in order
    E <- order(-abs(eigen(L)$values))

    ## Go into eigenspace
    U <- eigen(L)$vectors[,E]
    U <- U[,c((ncol(U)-Q+1):ncol(U))]
    U <- U / rowSums(U^2)^(1/2)

    ## Applying the K-means in the eigenspace
    cl <- kmeans(U, Q, nstart = 10, iter.max = 30)$cluster

    ## handing lonely souls
    cl.final[connected] <- cl
    cl.final[unconnected] <- which.min(rowsum(D,cl))
    ##   cl.final[unconnected] <- sample(levels(factor(cl)), length(unconnected), rep=TRUE, prob=table(cl)/sum(table(cl)))

    return(as.factor(cl.final))
  } else {
    return(rep(1, ncol(A)))
  }
}

