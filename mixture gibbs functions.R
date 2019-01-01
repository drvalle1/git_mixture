#' Samples the cluster assignment for each location
#' 
#' This function samples the cluster assignment for each location
#' 
#' @param dat this matrix has L rows (e.g., locations) and S columns (e.g., species)
#'            and contains the presence-absence data 
#' @param one.minus.dat this matrix has L rows (e.g., locations) and S columns (e.g., species)
#'            and is calculated as 1-dat           
#' @param phi K x S matrix with the probability of observing each species in each group
#' @param theta vector of length K with the proportion of each location group
#' @param ngroup  maximum number of location groups
#' @param nloc number of locations
#' @param nspp number of species
#' @param return this function returns a vector with the cluster assignment of each location
#' @export

update.z=function(dat,one.minus.dat,phi,theta,ngroup,nloc,nspp,z){
  #pre-calculate some useful quantities
  log.theta=log(theta)
  log.phi=log(phi)
  log.one.minus.phi=log(1-phi)
  
  tab=rep(0,ngroup)
  tmp=table(z)
  tab[as.numeric(names(tmp))]=tmp
  
  #calculate the log probability for each possible group
  tmp=matrix(NA,nloc,ngroup)
  for (i in 1:ngroup){
    rasc=dat*matrix(log.phi[i,],nloc,nspp,byrow=T)+
         one.minus.dat*matrix(log.one.minus.phi[i,],nloc,nspp,byrow=T)
    tmp[,i]=rowSums(rasc)+log.theta[i] #sum log of prior probability
  }
  
  for (i in 1:nloc){
    tab[z[i]]=tab[z[i]]-1
    prob=rep(NA,ngroup)
    cond=tab==0
    prob[ cond]=-nspp*log(2)+log.theta[cond]
    prob[!cond]=tmp[i,!cond]

    #get normalized probs
    tmp1=prob-max(prob) #for numerical stability
    tmp2=exp(tmp1) #exponentiate log probability
    prob=tmp2/sum(tmp2) #normalize to sum to 1

    #draw from multinomial distrib
    ind=rmultinom(1,size=1,prob=prob)
    ind1=which(ind==1)
    z[i]=ind1
    tab[ind1]=tab[ind1]+1
  }
  z  
}
#--------------------------------------------

#' Samples theta parameters
#' 
#' This function samples the v parameters, which are then used to calculate the theta parameters
#' 
#' @param z vector with cluster assignment of each location 
#' @param ngroup maximum number of location groups
#' @param gamma1 this is the truncated stick-breaking prior parameter for the 
#'                number of location groups. This value should be between 0 and 1, and
#'                small values enforce more parsimonius results (i.e., fewer groups)
#' @param return this function returns a vector with the theta parameters
#' @export
#' 
update.theta=function(z,ngroup,gamma1,burnin,gibbs.step,theta,phi){
  #re-order thetas. Based on that, re-order z's
  if(gibbs.step<burnin & gibbs.step%%50==0){
    ind=order(theta,decreasing=T)
    theta=theta[ind]
    phi=phi[ind,]
    
    #get z.new
    z.new=z; z.new[]=NA
    for (i in 1:ngroup){
      cond=z==ind[i]
      z.new[cond]=i
    }
    z=z.new
  }

  #calculate the number of locations assigned to each group
  tmp=table(z)
  nk=rep(0,ngroup)
  nk[as.numeric(names(tmp))]=tmp

  #sample v from a beta distribution
  n.greater.k=cumsum(nk[ngroup:1])[ngroup:1]
  v=rbeta(ngroup-1,nk[-ngroup]+1,n.greater.k[-1]+gamma1)
  v1=c(v,1)
  
  #get theta from v1 using the stick-breaking equation 
  theta=rep(NA,ngroup)
  tmp=1
  for (i in 1:ngroup){
    theta[i]=v1[i]*tmp
    tmp=tmp*(1-v1[i])
  }
  
  #to avoid numerical issues
  cond=v>0.99999
  v[cond]=0.99999
  
  list(theta=theta,z=z,v=v,phi=phi)
}

#----------------------------
sample.gamma=function(v,ngroup,gamma.possib){
  ngamma=length(gamma.possib)
  soma=sum(log(1-v[-ngroup]))
  k=(ngroup-1)*(lgamma(1+gamma.possib)-lgamma(gamma.possib))
  res=k+(gamma.possib-1)*soma
  # sum(dbeta(v[-ngroup],1,gamma.possib[5],log=T))
  res=res-max(res)
  res1=exp(res)
  res2=res1/sum(res1)
  tmp=rmultinom(1,size=1,prob=res2)
  ind=which(tmp==1)
  gamma.possib[ind]
}