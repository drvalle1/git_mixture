# rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(1)

setwd('U:\\GIT_models\\git_mixture')
source('mixture gibbs functions.R')
sourceCpp('aux1.cpp')
dat=data.matrix(read.csv('fake data.csv',as.is=T))
one.minus.dat=1-dat

#useful metrics
nloc=nrow(dat)
nspp=ncol(dat)
ngroup=10
gamma1=0.1

#initial values
z=sample(1:ngroup,size=nloc,replace=T)
theta=rep(1/ngroup,ngroup)
tmp=runif(ngroup*nspp)
phi=matrix(tmp,ngroup,nspp)

#for gibbs sampler
ngibbs=1000
store.phi=matrix(NA,ngibbs,nspp*ngroup)
store.theta=matrix(NA,ngibbs,ngroup)
store.logl=rep(NA,ngibbs)
for (i in 1:ngibbs){
  print(i)
  z=update.z(dat=dat,one.minus.dat=one.minus.dat,
             phi=phi,theta=theta,
             ngroup=ngroup,nloc=nloc,nspp=nspp)
  # z=z.true
  
  tmp=ncs(dat=dat,z=z-1,nspp=nspp,nloc=nloc,ngroup=ngroup)
  phi=matrix(rbeta(ngroup*nspp,tmp$ncs1+1,tmp$ncs0+1),ngroup,nspp)
  
  theta=update.theta(z=z,ngroup=ngroup,gamma1=gamma1)

  #get logl
  prob=phi[z,]
  logl=dbinom(dat,size=1,prob=prob,log=T)
  
  #store results
  store.phi[i,]=phi
  store.theta[i,]=theta
  store.logl[i]=sum(logl)
}

