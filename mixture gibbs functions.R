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

update.z=function(dat,one.minus.dat,phi,theta,ngroup,nloc,nspp){
  #pre-calculate some useful quantities
  log.theta=log(theta)
  log.phi=log(phi)
  log.one.minus.phi=log(1-phi)
  
  #calculate the log probability for each possible group
  tmp=matrix(NA,nloc,ngroup)
  for (i in 1:ngroup){
    rasc=dat*matrix(log.phi[i,],nloc,nspp,byrow=T)+
         one.minus.dat*matrix(log.one.minus.phi[i,],nloc,nspp,byrow=T)
    tmp[,i]=rowSums(rasc)+log.theta[i] #sum log of prior probability
  }
  tmp1=tmp-apply(tmp,1,max) #for numerical stability
  tmp2=exp(tmp1) #exponentiate log probability
  prob=tmp2/rowSums(tmp2) #normalize to sum to 1
  
  #sample cluster assignments based on my customized multinomial distribution function
  rmultinom1(prob=prob,randu=runif(nloc))+1
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
update.theta=function(z,ngroup,gamma1){
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
  theta
}