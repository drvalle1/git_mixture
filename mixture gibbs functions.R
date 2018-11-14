update.z=function(dat,one.minus.dat,phi,theta,ngroup,nloc){
  log.theta=log(theta)
  log.phi=log(phi)
  log.one.minus.phi=log(1-phi)
  
  tmp=matrix(NA,nloc,ngroup)
  for (i in 1:ngroup){
    rasc=dat*matrix(log.phi[i,],nloc,nspp,byrow=T)+
         one.minus.dat*matrix(log.one.minus.phi[i,],nloc,nspp,byrow=T)
    tmp[,i]=rowSums(rasc)+log.theta[i]
  }
  tmp1=tmp-apply(tmp,1,max)
  tmp2=exp(tmp1)
  prob=tmp2/rowSums(tmp2)
  
  #sample cs
  rmultinom1(prob=prob,randu=runif(nloc))+1
}
#--------------------------------------------
update.theta=function(z,ngroup,gamma1){
  tmp=table(z)
  nk=rep(0,ngroup)
  nk[as.numeric(names(tmp))]=tmp
  n.greater.k=cumsum(nk[ngroup:1])[ngroup:1]
  v=rbeta(ngroup-1,nk[-ngroup]+1,n.greater.k[-1]+gamma1)
  v1=c(v,1)
  
  #get theta
  theta=rep(NA,ngroup)
  tmp=1
  for (i in 1:ngroup){
    theta[i]=v1[i]*tmp
    tmp=tmp*(1-v1[i])
  }
  theta
}