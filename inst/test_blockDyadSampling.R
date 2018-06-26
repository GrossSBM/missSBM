rm(list = ls())
library(missSBM)
library(mclust)

n <- 200
Q <- 3
alpha <- rep(1,Q)/Q
pi <- diag(.5,Q) + .1
directed <- FALSE

sbm <- simulateSBM(n, alpha, pi, directed)
adjacencyMatrix <- sbm$adjacencyMatrix
samplingParameters <- diag(.35,Q) + .05
samplingParameters[1,] = samplingParameters[,1] <-c(.1,.8,.8)
samplingParameters[2,] = samplingParameters[,2] <-c(.1,.8,.8)

sampling <- "block_dyad"
sampledAdjMatrix <- samplingSBM(adjacencyMatrix, sampling, samplingParameters, clusters = sbm$memberships)
R <- sampledAdjMatrix$samplingMatrix*1


dim(R)
dim(sbm$adjacencyMatrix)
prop.table(table(as.vector(R),as.vector(sbm$adjacencyMatrix)),1)

# Sampling rate
(sum(sampledAdjMatrix$samplingMatrix)-n)/(n*(n-1))

# Inference
vBlocks <- 1:5
infer_nmar <- inferSBM(sampledAdjMatrix$adjacencyMatrix, vBlocks, sampling)
infer_mar  <- inferSBM(sampledAdjMatrix$adjacencyMatrix, vBlocks, "dyad")

# ICL
ICL_nmar <- sapply(infer_nmar$models, function(x) x$vICL)
plot(ICL_nmar)

ICL_mar <- sapply(infer_mar$models, function(x) x$vICL)
plot(ICL_mar)


## Error rho
sum((infer_nmar$models[[3]]$fittedSampling$parameters - samplingParameters)^2)/sum((infer_nmar$models[[3]]$fittedSampling$parameters)^2)

## Error pi
sum((infer_nmar$models[[3]]$fittedSBM$connectParam - pi)^2)/sum((infer_nmar$models[[3]]$fittedSBM$connectParam)^2)
sum((infer_mar$models[[3]]$fittedSBM$connectParam - pi)^2)/sum((infer_mar$models[[3]]$fittedSBM$connectParam)^2)

## ARI
adjustedRandIndex(sbm$memberships, infer_nmar$models[[3]]$fittedSBM$memberships)
adjustedRandIndex(sbm$memberships, infer_mar$models[[3]]$fittedSBM$memberships)



############# Test #############
# logB <- function(x,a){
#   zero <- .Machine$double.eps
#   if(a<zero) a <- zero
#   if(a>1-zero) a <- 1-zero
#   return(x*log(a) + (1-x)*log(1-a))
# }
#
# iclCalculation <- function(X, R, infer, vBlocks){
#   zero <- .Machine$double.eps
#   K <- length(vBlocks)
#   n <- nrow(X)
#   icl <- vector("numeric", length = K)
#   icl2 <- vector("numeric", length = K)
#   icl3 <- vector("numeric", length = K)
#   PEN  <- vector("numeric", length = K)
#   for (p in 1:K) {
#     cat('',p)
#     loglik <- 0
#     pen    <- p*(p+1)*log(n*(n-1)/2) + (p-1)*log(n)
#     Tau    <- infer$models[[p]]$fittedSBM$blocks
#     # for (i in 1:n) {
#     #   for (j in 1:n) {
#     #     if(!is.na(X[i,j])){
#     #       for (q in 1:p) {
#     #         for (l in 1:p) {
#     #           loglik <- loglik + Tau[i,q]*Tau[j,l]*( logB(X[i,j],infer$models[[p]]$fittedSBM$connectParam[q,l]) +
#     #                                                   logB(R[i,j],infer$models[[p]]$fittedSampling$parameters[q,l]) )
#     #         }
#     #       }
#     #     }
#     #   }
#     # }
#     Xo                     <- X; Xo[R==0] <- 0
#     prob                   <- Tau%*% infer$models[[p]]$fittedSBM$connectParam %*%t(Tau)
#     adjMatrix_zeroDiag     <- Xo ; diag(adjMatrix_zeroDiag) <- 0
#     adjMatrix_zeroDiag_bar <- 1 - Xo ; diag(adjMatrix_zeroDiag_bar) <- 0
#     loglikSBM              <- sum(Tau %*% log(infer$models[[p]]$fittedSBM$mixtureParam)) +  .5 * sum( adjMatrix_zeroDiag * log(prob) + adjMatrix_zeroDiag_bar *  log(1 - prob))
#
#     diag(R)        <- 0
#     R_bar          <- 1-R; diag(R_bar) <- 0
#     probSamp       <- Tau%*% infer$models[[p]]$fittedSampling$parameters %*%t(Tau)
#     probSamp[probSamp < zero] <- zero; probSamp[probSamp > 1 - zero] <- 1 - zero
#     loglikSampling <- .5 * sum( R * log(probSamp) + (R_bar) *  log(1 - probSamp))
#
#     icl[p]  <- -2*(loglikSBM + loglikSampling) + pen
#     icl2[p] <-  loglikSampling
#     icl3[p] <-  loglikSBM
#     PEN[p]  <-  pen
#   }
#   return(list(icl=icl, icl2=icl2, icl3=icl3, pen=PEN))
# }
#
# Icl <- iclCalculation(adjacencyMatrix, R, infer, vBlocks)
#
# plot(Icl$icl)
# points(Icl$icl2)
# points(Icl$icl3)


