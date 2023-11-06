## Code used to build figure 2
#skewness

Mskew<-function(p){
a1<-(p/(p+3))^(3/2)-3*p^3/(((p+1)*(p+2))^(3/2))+2*(p/(p+1))^(9/2)
b1<- ((p/(p+2))^(3/2)-(p/(p+1))^3)^(3/2)                           
return(a1/b1)
}
x<-seq(0.05,5,0.1)
df<-data.frame(x)
library(ggplot2)
ggplot(df,aes(x))+
  stat_function(fun=function(x) Mskew(x),size=1.3)+
  labs(x="p",y="")+ theme_bw(base_size = 22)


#######################################




 

Mkurt<-function(p){
  a2<-(p/(p+4))^(3/2)-4*p^3/(((p+1)*(p+3))^(3/2))+6*p^(9/2)/(((p+1)*(p+2))^(3/2))-3*((p/(p+1))^6)
  b2<- ((p/(p+2))^(3/2)-(p/(p+1))^3)^(2)                          
  return(a2/b2)
}

x<-seq(0.05,2,0.1)
df<-data.frame(x)

ggplot(df,aes(x))+
  stat_function(fun=function(x) Mkurt(x),size=1.3)+
  labs(x="p",y="")+ theme_bw(base_size = 22)
