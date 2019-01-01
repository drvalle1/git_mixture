rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(1)

setwd('U:\\GIT_models\\git_mixture')
source('mixture gibbs functions.R')
sourceCpp('aux1.cpp')
source('mixture gibbs main function.R')

ngibbs=1000
dat=data.matrix(read.csv('fake data.csv',as.is=T))
    
res=mixture.gibbs.main.func(dat=dat,ngroup=50,ngibbs=ngibbs,burnin=ngibbs/2)
    
