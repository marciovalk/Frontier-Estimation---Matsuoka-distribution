## Code used to build figure 3


source("functions.R")

# Standard normal example
# https://e-archivo.uc3m.es/bitstream/handle/10016/34168/expectile_JOMA_2021.pdf?sequence=1&isAllowed=y
#page 4
#x=0.43633
#alpha=0.75
#alpha/(1-alpha)
#pnorm2<-function(y)(1-pnorm(y))
#fn1<-function(y) integrate(pnorm,-Inf,y)$value
#fn2<-function(y) integrate(pnorm2,y,Inf)$value
#fn1(x)
#3*fn2(x)

###################################
require(expint)
mcdf1<-function(x,p) (2/sqrt(pi))*gammainc(3/2,-(p*log(x)))
#mcdf1(0.99,2)
mcdf2<-function(x,p) (1-(2/sqrt(pi))*gammainc(3/2,-(p*log(x))))
#mcdf1(0.99,2)
#mcdf2(0.99,2)

#x<-seq(0.001,0.99,length.out=1000)
#y<-mapply(mcdf,x,p=5)
#plot(x,y)

il_mcdf<-function(y,p) integrate(mcdf1,p=p,0,y)$value #integrate left mcdf
ir_mcdf<-function(y,p) integrate(mcdf2,p=p,y,1)$value #integrate right mcdf

#mcdf(0.3,2)
#il_mcdf(0.3,2)
#ir_mcdf(0.3,2)


alpha_expec<-function(alpha,p){
  f<-function(x,alpha,p){return(abs((1-alpha)*il_mcdf(x,p)-alpha*ir_mcdf(x,p)))}
  o = optimize(f,alpha=alpha, p=p, interval=c(0.0015, 0.998))
  return(o$minimum)
}

#
vp=c(0.5,1,2,5,10)
alpha<-seq(0.001,0.99,length.out=100*length(vp))
ma<-matrix(0,ncol=3,nrow=1)
maaux<-matrix(0,ncol=3,nrow=length(alpha))
for(p in vp){
y<-mapply(alpha_expec,alpha,p=p)
maaux[,1]<-alpha
maaux[,2]<-y
maaux[,3]<-rep(p,length(alpha))
ma<-rbind(ma,maaux)
}
ma<-ma[-1,]
colnames(ma)<-c("alpha","expec","p")
madf<-data.frame(ma)

library(ggplot2)

ggplot(madf,aes(y=expec,x=alpha,colour=as.factor(p),group=p,cex=4))+geom_line(size=1.2) + 
  labs(colour='p',x=bquote(~alpha),y=bquote(~epsilon[alpha](X)))+ theme(axis.text = element_text(size = 20),
                          legend.text = element_text(size = 24),
                          legend.title = element_text(size = 24),
                          axis.title.x = element_text( size=26),
                          axis.title.y = element_text( size=26))
