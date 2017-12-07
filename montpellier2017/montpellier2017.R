library(igraph)

# Simulation des graphes :

graph <- function(pir,top){
  n <- 300
  Q_com <- 3
  Q_hub <- 6
  alpha_com <- rep(1,3)/3
  alpha_hub <- rep(c(.2, .8)/3, 3)
  pia <- 3*pir
  pi_com <- matrix(c(pia, pir ,pir, pir, pia, pir, pir, pir, pia),3,3)
  pi_hub <- matrix(c(pia,pia,pir,pir,pir,pir, pia,pir,pir,pir,pir,pir, pir,pir,pia,pia,pir,pir, pir,pir,pia,pir,pir,pir,
                     pir,pir,pir,pir,pia,pia, pir,pir,pir,pir,pia,pir),6,6)
  
  
  if(top==1){
    Z <- 1
    pi <- pir
  } else if(top==2){
    z <- c(rep(1, 100), rep(2, 100), rep(3, 100))
    Z <- matrix(0, n, Q_com); Z[cbind(1:n,z)] <- 1
    pi <- pi_com
  } else if(top==3){
    Z <- t(rmultinom(n, size = 1, prob = alpha_com))
    pi <- pi_com
  } else if(top==4){
    z <- c(rep(1, 20), rep(2, 80), rep(3, 20), rep(4, 80), rep(5, 20), rep(6, 80))
    Z <- matrix(0, n, Q_hub); Z[cbind(1:n,z)] <- 1
    pi <- pi_hub
  } else if(top == 5){
    Z <- t(rmultinom(n, size = 1, prob = alpha_hub))
    pi <- pi_hub
  }
  
  Yvec    <- rbinom(n^2,1,Z %*% pi %*% t(Z))
  X       <- matrix(Yvec,n)
  diag(X) <- 0
  
  return(X)
}

g <- graph(.15, 5)




