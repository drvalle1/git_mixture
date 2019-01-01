theta=res$theta[nrow(res$theta),]
z=res$z[nrow(res$z),]
phi=matrix(res$phi[nrow(res$phi),],50,nspp)

plot(theta,type='h')
sum(theta>0.01)

table(z.true,z)

ind=c(2,4,9,7,3)
plot(phi.true,phi[ind,])