#' @import R6
SBM_VEMfit <-
  R6Class(classname = "SBM_VEMfit",
    public = list(
      completedNetwork = NULL, # the completed adjacency matrix of the initial network
      init             = NULL, # kind of classification for the beginning
      SBM              = NULL, #
      sampling         = NULL, #
      sampledNetwork   = NULL, #
      blockVarParam    = NULL, # variational parameters for latent blocks, (a.k.a. tau)
      taylorVarParam   = NULL, # variational parameters for Taylor expansion (a.k.a. ksi)
      controlVEM       = NULL, # VEM paramater
      maxIterVEM       = NULL, # VEM paramater
      lowerBound       = NULL, # variational lower bound (a.k.a. J)
      compLogLik       = NULL, # variationnal complete log-likelihood
      vICL             = NULL, # compute the (variational) integrated complete likelihood
      blockInit        = NULL, # Initialization clustering8
      VEstep           = function() {

        if (!(class(self$sampling)[1] %in% c("sampling_randomPairMAR", "sampling_randomNodesMAR", "sampling_snowball"))) {

          self$blockVarParam    <- self$SBM$fixPoint(self$SBM, self$blockVarParam, self$completedNetwork)
          self$completedNetwork <- self$SBM$updateNu(self$SBM, self$sampling, self$sampledNetwork, self$blockVarParam, self$completedNetwork, self$taylorVarParam)
          if(class(self$sampling)[1] == "sampling_starDegree"){
            self$taylorVarParam   <- self$SBM$updateKsi(self$sampling, self$completedNetwork, self$sampledNetwork)
          }
          self$lowerBound       <- c(self$lowerBound, self$SBM$completeLogLik(self$completedNetwork, self$blockVarParam)
                                     - sum(self$blockVarParam*log(self$blockVarParam + 1*(self$blockVarParam==0))))
          self$compLogLik       <-  c(self$compLogLik, self$SBM$completeLogLik(self$completedNetwork, self$blockVarParam))

          self$sampling$missingParam <- self$sampling$updatePsi(self$completedNetwork, self$sampledNetwork,  self$blockVarParam, self$taylorVarParam)
        } else {
          ## MAR SAMPLING DESIGNS
          self$blockVarParam    <- self$SBM$fixPoint_MAR(self$SBM, self$blockVarParam, self$completedNetwork, self$sampledNetwork$samplingMatrix)
          self$lowerBound       <- c(self$lowerBound, self$SBM$completeLogLik_MAR(self$blockVarParam, self$sampledNetwork)
                                     - sum(self$blockVarParam*log(self$blockVarParam + 1*(self$blockVarParam==0))))
          self$compLogLik       <- c(self$compLogLik, self$SBM$completeLogLik_MAR(self$blockVarParam, self$sampledNetwork))

          if(!(class(self$sampling)[1] == "sampling_snowball")){
            self$sampling$missingParam <- self$sampling$updatePsi(self$completedNetwork, self$sampledNetwork,  self$blockVarParam, self$taylorVarParam)
          }
        }
      },
      Mstep            = function() {
        if(!(class(self$sampling)[1] %in% c("sampling_randomPairMAR", "sampling_randomNodesMAR", "sampling_snowball"))){
          self$SBM <- self$SBM$maximization(self$SBM, self$completedNetwork, self$blockVarParam)
        } else {
          self$SBM <- self$SBM$maximization_MAR(self$SBM, self$completedNetwork, self$blockVarParam, self$sampledNetwork$samplingMatrix)
        }
      },
      initialize       = function(SBM, sampledNetwork, sampling, init = "spectralC", blockInit = NULL, controlVEM = 1e-5, maxIterVEM = 1000) {
        self$sampledNetwork   <- sampledNetwork
        self$init             <- init
        self$SBM              <- SBM
        self$sampling         <- sampling
        self$completedNetwork <- sampledNetwork$adjacencyMatrix
        self$controlVEM       <- controlVEM
        self$maxIterVEM       <- maxIterVEM
        self$blockInit        <- blockInit
      }
    )
  )

SBM_VEMfit$set("public", "SpectralClustering",
  function() {
    if (self$SBM$nBlocks > 1) {
      ## basic handling of missing values
      if (anyNA(self$completedNetwork)) self$completedNetwork[is.na(self$completedNetwork)] <- 0

      ## handling lonely souls
      cl.final <- rep(NA, ncol(self$completedNetwork))
      unconnected <- which(rowSums(self$completedNetwork) == 0)
      connected <- setdiff(1:ncol(self$completedNetwork), unconnected)

      self$completedNetwork <- self$completedNetwork[connected,connected]
      if (self$SBM$nBlocks > 1) {

        ## Normalized Laplacian
        D <- colSums(self$completedNetwork)
        L <- diag(rep(1,ncol(self$completedNetwork))) -
          diag(D^(-1/2)) %*% self$completedNetwork %*% diag(D^(-1/2))

        ## Absolute eigenvalue in order
        E <- order(-abs(eigen(L)$values))

        ## Go into eigenspace
        U <- eigen(L)$vectors[,E]
        U <- U[,c((ncol(U)-self$SBM$nBlocks + 1):ncol(U))]
        U <- U / rowSums(U^2)^(1/2)

        ## Applying the K-means in the eigenspace
        cl <- kmeans(U, self$SBM$nBlocks, nstart = 10, iter.max = 30)$cluster
      } else {
        cl <- as.factor(rep(1,ncol(self$completedNetwork)))
      }

      ## handing lonely souls
      cl.final[connected] <- cl
      cl.final[unconnected] <- which.min(rowsum(D,cl))

      return(as.factor(cl.final))
    } else {
      return(rep(1, self$SBM$nNodes))
    }
  }
)

SBM_VEMfit$set("public", "initialization",
  function() {
    if (self$init == "CAH") {
      networkWithZero <- self$completedNetwork; networkWithZero[is.na(networkWithZero)] <- 0
      D  <- as.matrix(dist(networkWithZero, method = "manhattan"))
      D[networkWithZero == 1] <- D[networkWithZero == 1] - 2
      classif <- cutree(hclust(as.dist(D), method = "ward.D"), self$SBM$nBlocks)
      self$blockVarParam <- matrix(0,self$SBM$nNodes,self$SBM$nBlocks)
      self$blockVarParam[cbind(1:self$SBM$nNodes, classif)] <- 1
    } else {
      cl0 <- self$SpectralClustering()
      self$blockVarParam <- matrix(0,self$SBM$nNodes,self$SBM$nBlocks) ; self$blockVarParam[cbind(1:self$SBM$nNodes, cl0)] <- 1
    }

    if (class(self$sampling)[1] == "sampling_starDegree") {
      self$completedNetwork[is.na(self$completedNetwork)] <- mean(self$completedNetwork[which(!is.na(self$completedNetwork))])
      networkWithZeros    <- self$completedNetwork
      networkWithZeros[self$sampledNetwork$missingDyads] <- 0
      Dtilde              <- rowSums(self$completedNetwork)
      Dchap               <- rowSums((self$completedNetwork - networkWithZeros)*(1 - (self$completedNetwork - networkWithZeros))) + Dtilde^2
      self$taylorVarParam <- sqrt(self$sampling$missingParam[1]^2 + (self$sampling$missingParam[2]^2)*Dchap + 2*self$sampling$missingParam[1]*self$sampling$missingParam[2]*Dtilde)
    }
    self$completedNetwork[is.na(self$completedNetwork)] <- 0
    theta <- (t(self$blockVarParam) %*% self$completedNetwork %*% self$blockVarParam) / (t(self$blockVarParam) %*% ((1 - diag(self$SBM$nNodes))) %*% self$blockVarParam)
    if (!(class(self$SBM)[2] %in% c("SBM_PoissonDirected", "SBM_PoissonUndirected")) & !(class(self$sampling)[1] == "sampling_snowball")) {
      self$completedNetwork[self$sampledNetwork$missingDyads] <- ((self$blockVarParam) %*% theta %*% t(self$blockVarParam))[self$sampledNetwork$missingDyads]
    }
  }
)


SBM_VEMfit$set("public", "doVEM",
               function() {

                 conv    <- vector("numeric", self$maxIterVEM) ; conv[1] <- NA
                 theta   <- vector("list", length = self$maxIterVEM)

                 self$initialization()

                 i <- 0; cond <- FALSE
                 while(!cond){
                   i <- i+1

                   self$Mstep()
                   self$VEstep()

                   theta[[i]] <- self$SBM$connectParam

                   if (i > 1) {
                     conv[i] <- sqrt(sum((theta[[i]]-theta[[i-1]])^2)) / sqrt(sum((theta[[i-1]])^2))
                     cond    <- (i > self$maxIterVEM) |  (conv[i] < self$controlVEM)
                   }
                 }
                 if(!(class(self$sampling)[1] == "sampling_snowball")){
                   self$vICL <- -2 * (self$compLogLik[length(self$compLogLik)] + self$sampling$samplingLogLik(self$sampledNetwork, self$completedNetwork, self$blockVarParam)) +
                     self$sampling$penality(self$SBM$nBlocks)
                 } else {
                   self$vICL <- -2 * self$compLogLik[length(self$compLogLik)] + self$sampling$penality(self$SBM$nBlocks)

                 }
               }
)

SBM_VEMfit$set("public", "doVEMPoisson",
  function() {

    conv    <- vector("numeric", self$maxIterVEM) ; conv[1] <- NA
    theta   <- vector("list"   , self$maxIterVEM)

    self$initialization()

    i <- 0; cond <- FALSE
    while (!cond) {
      i <- i + 1

      self$Mstep()
      self$VEstep()

      theta[[i]] <- self$SBM$connectParam

      if (i > 1) {
        conv[i] <- sqrt(sum((theta[[i]]-theta[[i-1]])^2)) / sqrt(sum((theta[[i-1]])^2))
        cond    <- (i > self$maxIterVEM) |  (conv[i] < self$controlVEM)
      }
    }
    self$vICL <- -2 * self$compLogLik[length(self$compLogLik)] + self$sampling$penalityPoisson(self$SBM$nBlocks, self$sampledNetwork$adjacencyMatrix)
  }
)

