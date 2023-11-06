library(zipfR)
library(foreach);library(doParallel)
library(matrixStats)
library(xtable)

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
            
#Classical backfitting estimator
cbe.explicit.bivar=function(X, y, h1, k1 = k.epa,details=FALSE,f.sm=ll.sm){
  #obs: nessa função, não normalizei X_1 e X_2 em [0,1]
  n <- length(y)
  a0 <- mean(y)
  x1=X[,1]
  x2=X[,2]
  if(length(h1)==1) h1=rep(h1,2)
  S1=f.sm(xd=x1,zz=x1,h=h1[1],k=k1)
  S2=f.sm(xd=x2,zz=x2,h=h1[2],k=k1)
  S1=scale(S1,scale=FALSE)   #centering by column (I-11'/n)S_1
  S2=scale(S2,scale=FALSE)  #apply(S2, 2, function(y) y - mean(y))
  #if(all.equal(f.sm,ll.sm.loo)){
  #  diag(S1)=0           #se usar leave-one-out, forço a diagonal a ser nula
  #  diag(S2)=0           #ao multiplicar S%*%y, a i-ésima linha de pesos S[i,]
                         #relativo ao x[i], calculada a partir de x[-i], irá anular y[i]
  #}
  W1=diag(n)-solve(diag(n)-S1%*%S2)%*%(diag(n)-S1)
  W2=diag(n)-solve(diag(n)-S2%*%S1)%*%(diag(n)-S2)
  if(details){
    est1=W1%*%y
    est2=W2%*%y
    est=a0+(W1+W2)%*%y
    result=list("comp1"=est1,"comp2"=est2,"estimate"=est)
  }else{
    est=a0+(W1+W2)%*%y
    result=list("estimate"=est)
  }
  return(result)
}

#Loss function to be optimized for selecting the bandwidth
cbe.cv=function(X, y, h1, k1 = k.epa){
  ests=cbe.explicit.bivar(X=X, y=y, h1=h1, k1 = k1,details=FALSE,f.sm=ll.sm.loo)[[1]]
  res = sum((y - ests)^2)
  if(h1[1]<exp(-10)|h1[2]<exp(-10)){res=(abs(h1[1])+abs(h1[2])+1)*exp(10)} #penalizador para o optim não selecionar h<0
  if(h1[1]>1|h1[2]>1){res=(h1[1]+h1[2])*exp(10)}                           #penalizador para o optim não selecionar h>1
  return(res)
}

m1 <- function(x) {
  res <- -3/2*x^2+3*x-1
  return(res)
}                                    #Component 1 (X1): m_1(X_1)
m2 <- function(x) {
  res <- -log(x)/x+log(2)^2/2
  return(res)
}                                    #Component 2 (X2): m_2(X_2)


cbe.graph <- function(x, est_comp, true_comp,tam) {
  nx <- ncol(x)
  for (j in 1:nx) {
    x1 <- x[, j]
    m_est <- est_comp[, j]
    m_true <- true_comp[, j]
    
    plot(
      x = sort(x1),
      y = m_true[order(x1)],
      type = "l",
      ylab = "",lwd=2.5,
      font=6,cex.lab=tam+0.1,cex.axis=tam,font.lab=6,cex=tam,
      ylim = c(min(m_est, m_true), max(m_est, m_true)),
      xlab = paste(TeX("$X_$"), j, sep = "")
    )
    lines(sort(x1), m_est[order(x1)], col = 2,lwd=1.7)
    rug(sort(x1))
    if(j==1){
      legend("topleft",legend=c("true","estimated"),
             col=c(1,2),lty=1,bty = 'n',text.font = 6,
             cex=tam,lwd=c(2.5,1.7))
    }
    
  }
}     #Plot generator

destino="D:R_codes"         #### setting! ####
setwd(dir=destino)

############## Parallelization ###################
noc=5                                                             #n. cores
cl=makeCluster(noc)
registerDoParallel(cl)
clusterEvalQ(cl,c(library(zipfR),library(matrixStats)))           #load packages       
clusterExport(cl,c("k.epa","F.mv.i","rmv","ll.x.sm",
                   "ll.sm","ll.sm.loo","cbe.explicit.bivar","cbe.cv",
                   "m1","m2"))                                    #load custom functions
########### Simulação Monte Carlo ##########################

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
      X <- matrix(runif(2*n,1,2), nrow = n, ncol = 2)  #design matrix
      x1 <- X[, 1]                                     #Input 1
      x2 <- X[, 2]                                     #Input 2
      e <- rmv(n=n,p=p0)                               #Efficiency R~M(p)
      w=exp(-m1(x1) - m2(x2))                          #Production frontier f()
      y = w*e                                          #Multiplicative Model: Y=f_1(X_1)*f_2(X_2)*R
      ly=-log(y)                                       #Log-scale
      
      #bandwidth selection (leave-one-out)
      h1=optim(par=c(0.25,0.25),fn=cbe.cv,method = "Nelder-Mead",X=X,y=ly,k1=k.epa)$par 
      
      #Backfitting using local linear smoothers, the Epanechnikov kernel and the leave-one-out bandwidth
      g.hat2 = cbe.explicit.bivar(X=X, y=ly, h1 = h1, 
                                    k1=k.epa,f.sm = ll.sm,details=TRUE)   

      z.hat2=ly-g.hat2$estimate                                            #residuals
      p.mm2=sqrt(3/(2*mean(z.hat2^2)))                                     #p_hat   (estimate of p)
      f.as2=exp(3/(2*p.mm2)-g.hat2$estimate)                               #f_hat   (estimate of f)
      ase.g=mean((g.hat2$estimate-m1(x1)-m2(x2)-3/(2*p0))^2)               #Loss function: L(g_hat)
      ase.g1=mean((g.hat2$comp1-m1(x1))^2)                                 #Loss function: L(g_hat1)
      ase.g2=mean((g.hat2$comp2-m2(x2))^2)                                 #Loss function: L(g_hat2)
      ase.f=mean((f.as2-w)^2)                                              #Loss function: L(f_hat)
      
      c(p.mm2,ase.g1,ase.g2,ase.g,ase.f,h1)                                #vector of results
    }
    resultados2[[q]]=result  
    
    if(q==1){p1="dot3"}else{p1=p0}
    write.csv(resultados2[[q]],paste("ase_p_bivar",p1,"n",n,".csv",sep = ""))
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
tabela=cbind(medias[,1:5],varia[,1],q5,q95)
ordem=c(5,4,2,3,1,6,7,8)
tabela=tabela[,ordem]

##### Save results ####
write.csv(tabela,"tabela_bivar.csv")    

##### Show table for Latex #############
print(xtable(tabela,digits=2),include.rownames = FALSE)
