rm(list=ls()) 
library(missSBM)
library(parallel)
library(ggplot2)

netSize <- c(100, 200, 400, 600, 800, 1000)
Q <- 3
alpha <- rep(1,Q)/Q                                        
pi <- diag(.45,Q) + .05                                                           
family <- "Bernoulli"                                                              
directed <- FALSE                                                                  

samp <- c("MAREdge", "doubleStandard")
res <- data.frame()

for(n in netSize){
  mySBM <- simulateSBM(n, alpha, pi, family, directed)                           
  adjacencyMatrix <- mySBM$adjacencyMatrix 
  for(sampling in samp){
    if(sampling == "MAREdge"){
      samplingParameters <- .5
    } else {
      samplingParameters <- c(.5, .5)
    }
    res <- rbind(res,do.call(rbind, mclapply(1:40, function(i){
      sampledAdjMatrix <- samplingSBM(adjacencyMatrix, sampling, samplingParameters)
      T1 <- Sys.time() 
      sbm <- inferSBM(sampledAdjMatrix, Q, sampling, family, directed, plot = FALSE)
      T2 <- Sys.time() 
      
      return(data.frame(nodesNumber = n, time = difftime(T2, T1, units = "secs"), sampling = sampling))
    }
    , mc.cores = 8)))
  }
}

save(res, file = "~/Git/missSBM/inst/tempsDeCalcul.RData")
# ggplot(res, aes(x=factor(nodesNumber), y=as.numeric(time))) + geom_boxplot() + facet_grid(.~sampling)



