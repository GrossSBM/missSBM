library(igraph)
library(ggplot2)
library(missSBM)
source("~/Git/missSBM/montpellier2017/sampling_function.R")

#### Simulation des graphes :

graph <- function(pir,dens,top){
  n <- 300
  Q_com <- 3
  Q_hub <- 6
  alpha_com <- rep(1,3)/3
  alpha_hub <- rep(c(.2, .8)/3, 3)
  if(!is.null(dens)){
    pir <- switch(top,
                  "1" = dens,
                  "2" = (9/15)*dens,
                  "3" = (9/15)*dens,
                  "4" = (25/31)*dens,
                  "5" = (25/31)*dens)
  }
  pia <- 3*pir
  pi_com <- matrix(c(pia, pir ,pir, pir, pia, pir, pir, pir, pia),3,3)
  pi_hub <- matrix(c(pia,pia,pir,pir,pir,pir, pia,pir,pir,pir,pir,pir, pir,pir,pia,pia,pir,pir, pir,pir,pia,pir,pir,pir,
                     pir,pir,pir,pir,pia,pia, pir,pir,pir,pir,pia,pir),6,6)


  if(top=="1"){
    Z <- 1
    pi <- pir
  } else if(top=="2"){
    z <- c(rep(1, 100), rep(2, 100), rep(3, 100))
    Z <- matrix(0, n, Q_com); Z[cbind(1:n,z)] <- 1
    pi <- pi_com
  } else if(top=="3"){
    Z <- t(rmultinom(n, size = 1, prob = alpha_com))
    pi <- pi_com
  } else if(top=="4"){
    z <- c(rep(1, 20), rep(2, 80), rep(3, 20), rep(4, 80), rep(5, 20), rep(6, 80))
    Z <- matrix(0, n, Q_hub); Z[cbind(1:n,z)] <- 1
    pi <- pi_hub
  } else if(top == "5"){
    Z <- t(rmultinom(n, size = 1, prob = alpha_hub))
    pi <- pi_hub
  }

  Yvec    <- rbinom(n^2,1,Z %*% pi %*% t(Z))
  X       <- matrix(Yvec,n)
  diag(X) <- 0

  return(list(matAdj = X, Z=Z))
}

A <- graph(dens=.02, top="2")
G=graph_from_adjacency_matrix(A$matAdj,mode="directed")
plot(G)
summary(degree(G))

#### Simulations 1 :

npv=100
nv=3

densities <- seq(.0005, .05, length = 50)
topologies <- as.character(1:5)
res <- data.frame()



for(top in topologies){
  for(dens in densities){
    for (k in 1:1){
    matAdj <- graph(dens=dens,top=top)$matAdj
    n=nrow(matAdj)
    SN1=snowball_village(npv,nv,1,10,matAdj)
    SN2=snowball_village(npv,nv,2,5,matAdj)
    res <- rbind.data.frame(res, data.frame(topology = top, density = dens,empdensity = sum(matAdj)/(n*n-n), samplingRate = c(length(SN1)/300, length(SN2)/300), step = c("One step","Two steps")))
    }
      }
}

#### Représentation graphique :
boxplot(res$density-res$empdensity~res$topology)
ggplot(res, aes(x=density, y=samplingRate, color=topology, linetype = step)) + geom_smooth(se=FALSE) #+ geom_point() #+ facet_grid(. ~ step) + theme_bw(base_size = 20)



#### Simulation 2 :

npv=100
nv=3

samplingRate <- c(.2,.4,.6)
densities <- c(.01) # .005
topologies <- as.character(2:5)
res <- data.frame()

for(dens in densities){
  for(top in topologies){
    Q <- ifelse(top %in% c("2", "3"), Q <- 3, Q <- 6)
    for(sampR in samplingRate){
      res <- rbind(res,do.call(rbind, lapply(1:1, function(i){
        g <- graph(dens=dens,top=top)
        matAdj <- g$matAdj

        ### SN0 ###
        nbreI0 <- sampR*300
        SN0    <- sample(1:300, nbreI0, replace = F)
        matAdj_N0 <- matAdj ; matAdj_N0[SN0,] <- NA ; matAdj_N0[,SN0] <- NA

        ### SN1 ###
        nbreI1 <- snowball_samplingrate(npv,nv,1,nbreI0,dens)
        SN1    <- snowball_village(npv,nv,1,ceiling(nbreI1),matAdj)
        matAdj_N1 <- matAdj ; matAdj_N1[SN1,] <- NA ; matAdj_N1[,SN1] <- NA

        ### SN2 ###
        nbreI2 <- snowball_samplingrate(npv,nv,2,nbreI0,dens)
        SN2    <- snowball_village(npv,nv,2,floor(nbreI2),matAdj)
        matAdj_N2 <- matAdj ; matAdj_N2[SN2,] <- NA ; matAdj_N2[,SN2] <- NA

        ### Inference classification ###
        VEM_SN0 <- SBM_collection$new(sample$adjacencyMatrix, Q, "MARNode", "Bernoulli", TRUE)
        VEM_SN1 <- SBM_collection$new(sample$adjacencyMatrix, Q, "MARNode", "Bernoulli", TRUE)
        VEM_SN2 <- SBM_collection$new(sample$adjacencyMatrix, Q, "MARNode", "Bernoulli", TRUE)

        return(data.frame(density = dens,
                          topology = top,
                          samplingRate = sampR,
                          sampling = c("SN0", "SN1", "SN2"),
                          ARI=c(adjustedRandIndex(apply(VEM_SN0$models[[1]]$blockVarParam, 1, which.max), g$Z %*% (1:Q)),
                                adjustedRandIndex(apply(VEM_SN1$models[[1]]$blockVarParam, 1, which.max), g$Z %*% (1:Q)),
                                adjustedRandIndex(apply(VEM_SN2$models[[1]]$blockVarParam, 1, which.max), g$Z %*% (1:Q)))))
      })))
    }
  }
}

