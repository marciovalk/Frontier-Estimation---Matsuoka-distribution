## Code used to build figure 5

library(expint)
require(lattice)


fnr<-function(x,p,n,r){
  f<-2*r*choose(n,r)*sqrt((-p^3*log(x))/pi)*x^(p-1)*( 2*gammainc(3/2,-p*log(x))/sqrt(pi))^(r-1)*( 
    1- 2*gammainc(3/2,-p*log(x))/sqrt(pi))^(n-r)
  f
}

#  x\in (0,1)

#############################################
# Minimum fixed n
r=1
pp=c(1,2,5,10)
#par(mfrow=c(2,2))
for(p in pp){
x=seq(0.03,0.98,0.02)
#p=seq(0.1,10,length.out =length(x))
n=1:length(x)
f <- outer(x, p=p, r=1,n,fnr)
sum(is.na(f))
op <- par(bg = "white")
#persp(x,n, f, theta = 30, phi = 30, expand = 0.5, col = "gray",main=paste("p=",p))
persp(x,n, f, theta = 30, phi = 30, expand = 0.5, col = "gray")
}

main1<- bquote(paste("Pdf of ", X[(1)], " for n") == .(paste(length(x), " samples")))

title(main=main1,outer=T,cex.main=1.5, line=-1)



#############################################
# Maximum fixed n

pp=c(1,2,5,10)
par(mfrow=c(2,2))
for(p in pp){
  x=seq(0.05,0.95,0.07)
  r=length(x)
  #p=seq(0.1,10,length.out =length(x))
  n=1:length(x)
  f <- outer(x, p=p, r=r,n,fnr)
  sum(is.na(f))
  op <- par(bg = "white")
  persp(x,n, f, theta = 20, phi = 10, expand = 0.5, col = "gray",main=paste("p=",p))
}
mainn<- bquote(paste("Pdf of ", X[(n)], " for n") == .(paste(length(x), " samples")))

title(main=mainn,cex.main=1.5,outer=T, line=-0.8)





#####################################################################
#####################################################################
#####################################################################
#####################################################################

#############################################
# Minimum fixed p
r=1
p=2
nn=c(10,20,40,70)
#par(mfrow=c(2,2))
for(n in nn){
  x=seq(0.05,0.95,length.out =n)
  ns=1:length(x)
  f <- outer(x, p=p, r=1,ns,fnr)
  sum(is.na(f))
  op <- par(bg = "white",mar = c(0.1, 0.1, 0.1, 0.1),cex=3.8)
  persp(x,ns, f, theta = 10, phi = 30, expand = 0.8, col = "gray",ylab="n",zlab="f")
}


main1<- bquote(paste("Pdf of ", X[(1)], " for p") == .(paste(p, ".")))
title(main=main1,outer=T,cex.main=1.5, line=-1)



#############################################
# Maximum fixed p

p=2
nn=c(5,10,20,40)
#par(mfrow=c(2,2))
for(n in nn){
  x=seq(0.05,0.95,length.out =n)
  r=length(x)
  ns=1:length(x)
  f <- outer(x, p=p, r=r,ns,fnr)
  sum(is.na(f))
  op <- par(bg = "white")
  persp(x,ns, f, theta =-10, phi = 20, expand = 0.5, col = "gray",main=paste("n=",n),ylab="n")
}

main1<- bquote(paste("Pdf of ", X[(n)], " for p") == .(paste(p, ".")))

title(main=main1,outer=T,cex.main=1.5, line=-1)

