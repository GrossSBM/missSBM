snowball_village=function(npv,nv,S,nfb,A)
#npv : number of nodes per village
#nv : number of village
#S : number of step in snowball
#nfb : number of nodes sampled in the first batch per village
#A : adjacency matrix of the network
{
  #first batch
  FB=sapply(1:nv,
            function(i){
                  return(sample((1+(i-1)*npv):(i*npv),nfb,replace = F))
            }
          )
  Nsampled=as.vector(FB)


  #snowball
  steptodo=S

  while (steptodo>0)
  {
    BN=which(rowSums(A[,Nsampled])>0)
    Nsampled=unique(c(Nsampled,BN))
    steptodo=steptodo-1
  }

  return(Nsampled)
}


# snowball_village(100,3,2,10,A)


#garantir taux d echantillonnage
snowball_samplingrate=function(npv,nv,S,NS,density)
#npv : number of nodes per village
#nv : number of village
#S : number of step in snowball
#NS :desired number of sampling nodes
#density : density of the network
{
  n=nv*npv #total number of nodes

  expect=function(nfirst)
  {
    return(n-(n-nfirst*nv)*(1-density)^(nfirst*nv))
  }

  if (S==1)
  {
   return( optimize(function(x) abs(NS-expect(x)),interval = c(1,npv))$minimum )
  }


  if (S==2)
  {
    expect2=function(nfirst)
    {
      NB1=expect(nfirst)-nfirst*nv
      if (n-nfirst*nv-NB1<0) stop("one wave is enough")
      return((n-nfirst*nv-NB1)*(1-(1-density)^NB1) )
    }
    return( optimize(function(x) abs(NS-expect(x)-expect2(x)),interval = c(1,npv))$minimum )
  }

}

# npv=100
# nv=3
#
# snowball_samplingrate(npv,nv,1,100,sum(A)/((npv*nv)^2-npv*nv))
# SN=snowball_village(npv,nv,1,6,A)
# length(SN)
#
# snowball_samplingrate(npv,nv,2,100,sum(A)/((npv*nv)^2-npv*nv))
# SN=snowball_village(npv,nv,2,1,A)
# length(SN)


