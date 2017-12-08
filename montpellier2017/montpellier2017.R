library(igraph)
library(ggplot2)
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
                  "2" = (9/15)*dens,     # (36/54)*dens,
                  "3" = (9/15)*dens,   #(36/54)*dens,
                  "4" = (25/31)*dens,  #(9/15)*dens,
                  "5" = (25/31)*dens)  # (9/15)*dens)
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

  return(X)
}

A <- graph(dens=.005, top="5")
G=graph_from_adjacency_matrix(A,mode="directed")
plot(G)
summary(degree(G))

#### Simulations :

npv=100
nv=3

densities <- seq(.005, .1, length = 50)
topologies <- as.character(1:5)
res <- data.frame()



for(top in topologies){
  for(dens in densities){
    for (k in 1:20){
    matAdj <- graph(dens=dens,top=top)
    n=nrow(matAdj)
    SN1=snowball_village(npv,nv,1,10,matAdj)
    SN2=snowball_village(npv,nv,2,1,matAdj)
    res <- rbind.data.frame(res, data.frame(topology = top, density = dens,empdensity = sum(matAdj)/(n*n-n), samplingRate = c(length(SN1)/300, length(SN2)/300), step = c("One step","Two steps")))
    }
      }
}

#### Représentation graphique :
boxplot(res$density-res$empdensity~res$topology)
ggplot(res, aes(x=density, y=samplingRate, color=topology, linetype = step)) + geom_smooth(se=FALSE) #+ geom_point() #+ facet_grid(. ~ step) + theme_bw(base_size = 20)

























