library(zipfR)
library(foreach);library(doParallel)
library(xtable)
library(matrixStats)

###### Functions ############
###########################
F.mv.i=function(q,p){                                         #quantil
  inversa=exp(-1/p*Igamma.inv(a=1.5, y=q*sqrt(pi)/2, lower=FALSE))
  return(inversa)
}                                   #Quantile function Matsuoka

rmv = function(n,p) return(F.mv.i(q = runif(n), p=p))       #Simulator Matsuoka

ll.x=function(y,xd,x0,h,k){
  n=length(y)
  u=xd-x0
  s=apply(cbind(1,u,u^2)*k(u/h),2,sum)/(n*h)
  m_hat=sum((s[3]-s[2]*u)*k(u/h)*y)/(n*h*(s[3]*s[1]-s[2]^2))
  return(m_hat)
}                             #Auxiliary function for local linear estimator

ll=function(y,xd,x0,h,k){
  m=sapply(x0,function(x) ll.x(y=y,xd=xd,x0=x,h=h,k=k))
  return(m)
}                               #Local linear estimator

k.epa = function(x) 0.75*(1-x^2)*(abs(x) <= 1)              #Epanechnikov kernel function

#Loss function to be minimized for selecting the bandwidth
cv.aux = function(hi,x,y,k=k.epa){
  ng = length(x)
  mh.loo = rep(NA,ng)
  for(i in 1:ng){
    mh.loo[i] = ll.x(y=y[-i],xd=x[-i],x0=x[i],h = hi,k=k)
  }
  res = sum((y - mh.loo)^2)
  return(res)
}   

m1 <- function(x) {
  res <- -x^2+4*x
  return(res)
}                                     #f(x)                         

destino="D:R_codes"         #### setting! ####
setwd(dir=destino)

############## Parallelization ###################
noc=5                                                            #n. cores
cl=makeCluster(noc)
registerDoParallel(cl)
clusterEvalQ(cl,c(library(zipfR),library(matrixStats)))          #load packages       
clusterExport(cl,c("k.epa","F.mv.i","rmv","ll",
                   "ll.x","as.band","cv.aux","m1"))              #load custom functions
#############################################################
################### MONTE CARLO #############################
ns=c(100,150,250)                                                #sample sizes
ps=c(1,2,8)                                                      #parameters p
nsim=1000                                                        #number of replications

final=final2=final3=final4=final5=vector("list",length(3))       #initializing lists
for(v in seq_along(ps)){                                         
  resultados2=vector("list",length(ps))
  p0=ps[v]                                                       #selected parameter p
  
  for(q in seq_along(ns)){                                         
    print(q)
    flush.console()
    
    n=ns[q]                                                      #selected sample size n
    
    result=foreach(j=1:nsim,.combine=rbind)%dopar%{              
      xc=runif(n,min=1,max=2)             #generating design points (inputs X)
      w=m1(xc)                            #production frontier f()
      r=rmv(n=n,p=p0)                     #generating efficiency values 
      y=r*w                               #Model: Y=f(X)*R
      ly=-log(y)                          #log-scale: Z=g(X)+eps
      
      h1=optimize(interval = c(0.001,1),f = cv.aux,x = xc,y = ly)$minimum  #Leave-one-out bandwidth selection
      g.hat2=ll(y=ly,xd=xc,x0=xc,h=h1,k=k.epa)                             #local linear estimate
      z.hat2=ly-g.hat2                                                     #residuals
      p.mm2=sqrt(3/(2*mean(z.hat2^2)))                                     #p_hat   (estimate of p)
      f.as2=exp(3/(2*p.mm2)-g.hat2)                                        #f_hat   (estimate of f)
      ase.g2=mean((g.hat2+log(w)-3/(2*p0))^2)                              #Loss function: L(g_hat)
      ase.f2=mean((f.as2-w)^2)                                             #Loss function: L(f_hat)
      
      c(p.mm2,ase.g2,ase.f2,h1)                                            #vector of results
    }
    resultados2[[q]]=result  
    
    if(q==1){p1="dot3"}else{p1=p0}
    write.csv(resultados2[[q]],paste("ase_p",p1,"n",n,".csv",sep = ""))
  }
  final[[v]]=sapply(resultados2,colMeans)                 #mean
  final2[[v]]=sapply(resultados2,colMedians)              #median
  final3[[v]]=sapply(resultados2,colVars)                 #variance
  final4[[v]]=sapply(resultados2,colQuantiles,probs=0.05) #quantiles
  final5[[v]]=sapply(resultados2,colQuantiles,probs=0.95)
}
stopCluster(cl)

#### Constructing tables #####
medias=do.call(rbind,lapply(final,t))
varia=do.call(rbind,lapply(final3,t))  
medianas=do.call(rbind,lapply(final2,t)) 
q5=do.call(rbind,lapply(final4,t)) 
q5=q5[,1]
q95=do.call(rbind,lapply(final5,t)) 
q95=q95[,1]
tabela=cbind(medias,varia[,1],medianas[,1],q5,q95)
ordem=c(3,2,1,5,6,7,8,4)
tabela=tabela[,ordem]

##### Save results ####
write.csv(tabela,"tabela.csv")           

##### Show table for Latex #############
print(xtable(tabela,digits=2),include.rownames = FALSE)

