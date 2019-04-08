data("war_graphs")
B2 = as.matrix(get.adjacency(war_graphs$beligerent))


# miss SBM full data
vBlocks <- 1:6
collection_sbm_full <-  missSBM::estimate(B2,vBlocks = vBlocks, "node", clusterInit = "hierarchical", cores = 2, trace = TRUE,control_VEM = list(threshold = 1e-1,fixPointIter=1,maxIter=20))
smooth(collection_sbm_full, "both", iterates = 1)
# Note that control_VEM does not seem to change anything...


# Blockmodels
library(blockmodels)
mod2b = BM_bernoulli("SBM_sym",B2)
mod2b$estimate()


mod2b$ICL*(-2) # different value of ICL considered
collection_sbm_full$ICL # but sampling considered...
# if i use this to focus on ICL coming from the SBM
sapply(collection_sbm_full$models,function(mod) mod$fittedSBM$vICL)




# If I choose 3 or 4 clusters i get different results for ICL computed values issued by missSBM
Q=4 # or Q=3
pimiss = collection_sbm_full$models[[Q]]$fittedSBM$connectParam
alphamiss = collection_sbm_full$models[[Q]]$fittedSBM$mixtureParam
Zmiss = collection_sbm_full$models[[Q]]$fittedSBM$blocks

pibloc = mod2b$model_parameters[[Q]]$pi
alphabloc = colMeans(mod2b$memberships[[Q]]$Z)
Zbloc = mod2b$memberships[[Q]]$Z


iclalamano = function(pi,alpha,Z,adj)
{
  adj0 = 1-adj
  diag(adj0) = 0
  vA = 0
  for (i in 2:nrow(adj))
    for (j in 1:(i-1))
      for (q in 1:ncol(Z))
        for (l in 1:ncol(Z))
          vA = vA + Z[i,q]*Z[j,l]*(log(pi[q,l])*adj[i,j]+log(1-pi[q,l])*adj0[i,j])
  return(vA + sum(Z%*%log(alpha)))
}

# ICL from estimated clusters model in missSBM
-2*iclalamano(pimiss,alphamiss,Zmiss,B2) + log(nrow(B2)) * (ncol(Zmiss)-1) + log(nrow(B2)*(nrow(B2)-1)/2)* (ncol(Zmiss)+1)*ncol(Zmiss)/2
# different from
collection_sbm_full$models[[Q]]$fittedSBM$vICL

# ICL from estimated clusters model 
-2*iclalamano(pibloc,alphabloc,Zbloc,B2) + log(nrow(B2)) * (ncol(Zbloc)-1) + log(nrow(B2)*(nrow(B2)-1)/2)* (ncol(Zbloc)+1)*ncol(Zbloc)/2


### This icl alamano will choose Q=3 clusters !!!

