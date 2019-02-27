rm(list=ls(all=TRUE))
set.seed(1)

#general settings
nloc=100
nspp=20
ngroup=5

#set parameters
theta=rep(1/ngroup,ngroup)
tmp=rmultinom(nloc,size=1,prob=theta)
z.true=z=apply(tmp==1,2,which)
phi.true=phi=matrix(rbeta(ngroup*nspp,1,1),ngroup,nspp)

#generate data
y=matrix(NA,nloc,nspp)
for (i in 1:nloc){
  phi1=phi[z[i],]
  y[i,]=rbinom(nspp,size=1,prob=phi1)
}
colnames(y)=paste0('spp',1:nspp)
rownames(y)=paste0('loc',1:nloc)

#export results
setwd('U:\\GIT_models\\git_mixture')
write.csv(y,'data_mixture.csv')

