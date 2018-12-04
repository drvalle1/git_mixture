#' Main function of the Mixture Model
#' 
#' Runs the Gibbs sampler and returns samples from the posterior distribution
#' 
#' @param dat this matrix has L rows (e.g., locations) and S columns (e.g., species)
#'            and contains the presence-absence data 
#' @param ngroup this is the maximum number of groups for locations
#' @param gamma1 this is the truncated stick-breaking prior parameter for the 
#'                number of location groups. This value should be between 0 and 1, and
#'                small values enforce more parsimonius results (i.e., fewer groups)
#' @param ngibbs  number of Gibbs sampler iterations              
#' @param return this function returns a list containing several matrices, all of which have ngibbs rows:
#'               - phi:   matrix with the estimated probability of observing each species in each group
#'               - theta: matrix nwith the estimated proportion of each location group
#'               - llk:   vector with the log-likelihood for each iteration
#'               - z:     matrix with the cluster assignment of each location
#' @export

mixture.gibbs.main.func=function(dat,ngroup,gamma1,ngibbs){
one.minus.dat=1-dat

#useful settings
nloc=nrow(dat)
nspp=ncol(dat)

#initial parameter values
z=sample(1:ngroup,size=nloc,replace=T)
theta=rep(1/ngroup,ngroup)
tmp=runif(ngroup*nspp)
phi=matrix(tmp,ngroup,nspp)

#to store results from gibbs sampler
store.phi=matrix(NA,ngibbs,nspp*ngroup)
store.theta=matrix(NA,ngibbs,ngroup)
store.z=matrix(NA,ngibbs,nloc)
store.logl=rep(NA,ngibbs)

#run gibbs sampler
for (i in 1:ngibbs){
  # print(i)
  z=update.z(dat=dat,one.minus.dat=one.minus.dat,
             phi=phi,theta=theta,
             ngroup=ngroup,nloc=nloc,nspp=nspp)
  
  #summarize data
  tmp=ncs(dat=dat,z=z-1,nspp=nspp,nloc=nloc,ngroup=ngroup)
  
  phi=matrix(rbeta(ngroup*nspp,tmp$ncs1+1,tmp$ncs0+1),ngroup,nspp)
  
  theta=update.theta(z=z,ngroup=ngroup,gamma1=gamma1)
  
  #get loglikelihood
  prob=phi[z,]
  logl=dbinom(dat,size=1,prob=prob,log=T)
  
  #store results
  store.phi[i,]=phi
  store.theta[i,]=theta
  store.logl[i]=sum(logl)
  store.z[i,]=z
}

#output results
list(phi=store.phi,theta=store.theta,logl=store.logl,z=store.z)
}
