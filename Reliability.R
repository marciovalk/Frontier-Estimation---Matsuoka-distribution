## Code used to build figure 4
source("functions.R")

#######################################
###### Estimating R=P(X>Y)
library(ggplot2)
library(reshape2)
# X ~ M(p)
# Y ~ M(q)
#p=5

qq=c(1,3,4,7,10)
R=c()
Rh=c()
for(i in 1:length(qq)){
  q=qq[i]
  Rhaux=c()
  for(j in 1:5){
    x<-mrdf(2000,p)
    y<-mrdf(2000,q)
    R_hat<-mean(x>y)
    Rhaux[j]=R_hat
  }
  
  s=p/(p+q)
  R[i]<-(2/pi)*((2*s-1)*sqrt(s*(1-s))+asin(sqrt(s)))
  Rh[i]=mean(Rhaux)
}

dataf2<-data.frame(q=qq,R=R,Rhat=Rh)
mdataf2 <- melt(dataf2,id.vars="q")

ggplot(mdataf2, aes(x = q,y=value, colour=variable)) + 
  geom_line(aes(linetype=variable))+geom_point()+
  ggtitle(" P(X<Y), X~M(q) and Y~ M(5) ") 


#####################################
v=c(0.5,1,2,5,10,15,30,50)
R=matrix(0,ncol=3,nrow=length(v)^2)
i=0
for(q in v){
for(p in v){
  i=i+1
  s=p/(p+q)
r<-(2/pi)*((2*s-1)*sqrt(s*(1-s))+asin(sqrt(s)))
R[i,]<-c(p,q,r)

}
}
colnames(R)<-c("p","q","R")
dfr<-data.frame(R)

ggplot(dfr,aes(y=R,x=q,colour=as.factor(p),group=p,cex=4))+geom_line(size=1.2) + 
  labs(colour='p')+ theme(axis.text = element_text(size = 20),
                          legend.text = element_text(size = 24),
                          legend.title = element_text(size = 24),
                          axis.title.x = element_text( size=26),
                          axis.title.y = element_text( size=26)
  )

ggplot(dfr,aes(y=R,x=p,colour=as.factor(q),group=q))+geom_line(size=1.2) + 
  labs(colour='q')+ theme(axis.text = element_text(size = 20),
                          legend.text = element_text(size = 24),
                          legend.title = element_text(size = 24),
                          axis.title.x = element_text( size=26 ),
                          axis.title.y = element_text( size=26)
  )

  
