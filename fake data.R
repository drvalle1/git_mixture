rm(list=ls(all=TRUE))
set.seed(1)

nloc=1000
nspp=200
ngroup=5

theta=rep(1/ngroup,ngroup)
tmp=rmultinom(nloc,size=1,prob=theta)
z.true=z=apply(tmp==1,2,which)

phi.true=phi=matrix(rbeta(ngroup*nspp,1,1),ngroup,nspp)

y=matrix(NA,nloc,nspp)
for (i in 1:nloc){
  phi1=phi[z[i],]
  y[i,]=rbinom(nspp,size=1,prob=phi1)
}

setwd('U:\\GIT_models\\git_mixture')
write.csv(y,'fake data.csv',row.names=F)

