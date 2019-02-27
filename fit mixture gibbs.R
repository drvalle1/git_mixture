rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(1)

setwd('U:\\GIT_models\\git_mixture')
source('mixture gibbs functions.R')
sourceCpp('aux1.cpp')
source('mixture gibbs main function.R')

#get data
dat=read.csv('data_mixture.csv',as.is=T)
rownames(dat)=dat$X
dat1=data.matrix(dat[,-1])

ngibbs=1000
res=mixture.gibbs.main.func(dat=dat1,ngroup=50,ngibbs=ngibbs,burnin=ngibbs/2)
    
