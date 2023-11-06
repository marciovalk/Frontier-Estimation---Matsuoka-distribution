library(Benchmarking)
library(akima)

destino="D:/Pesquisa/Pesquisa10_fronteiras/R_codes"         #### settings! ####
setwd(dir=destino)

###### Functions ############
########################### 


k.epa = function(x) 0.75*(1-x^2)*(abs(x) <= 1)          #kernel Epanechnikov

#Classical backfitting estimator
cbe.explicit.bivar=function(X, y, h1, k1 = k.epa,details=FALSE,f.sm=ll.sm){
  
  n <- length(y)
  a0 <- mean(y)
  x1=X[,1]
  x2=X[,2]
  if(length(h1)==1) h1=rep(h1,2)
  S1=f.sm(xd=x1,zz=x1,h=h1[1],k=k1)
  S2=f.sm(xd=x2,zz=x2,h=h1[2],k=k1)
  S1=scale(S1,scale=FALSE)   #centering by column (I-11'/n)S_1
  S2=scale(S2,scale=FALSE)   
  
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
  if(h1[1]<exp(-10)|h1[2]<exp(-10)){res=(abs(h1[1])+abs(h1[2])+1)*exp(10)} #penalty term to avoid optim to select h<0
  if(h1[1]>1|h1[2]>1){res=(h1[1]+h1[2])*exp(10)}                           #penalty term to avoid optim to select h>1
  return(res)
}    

ll.x=function(y,xd,x0,h,k){
  n=length(y)
  u=xd-x0
  s=apply(cbind(1,u,u^2)*k(u/h),2,sum)/(n*h)
  m_hat=sum((s[3]-s[2]*u)*k(u/h)*y)/(n*h*(s[3]*s[1]-s[2]^2))
  return(m_hat)
}                         #auxiliar function for local linear

ll=function(y,xd,x0,h,k){
  m=sapply(x0,function(x) ll.x(y=y,xd=xd,x0=x,h=h,k=k))
  return(m)
}                           #local linear estimator


######## Dados #############
data(milkProd)
N=nrow(milkProd)                #N=108
nc=milkProd$cows                #n cows
xc1=milkProd$vet/nc             #vet. expenses per cow (X)
x1=(xc1-min(xc1))/(max(xc1)-min(xc1))   #normalization 
xc2=milkProd$energy/nc          #energy expenses per cow (X)
x2=(xc2-min(xc2))/(max(xc2)-min(xc2))   #normalization
X=cbind(x1,x2)

y=milkProd$milk/nc/10^3         #milk production per cow (Y)
ly=-log(y)                      #Log-transformation (Z)

#bandwidth selection (leave-one-out)
h0=optim(par=c(0.5,0.5),fn=cbe.cv,method = "Nelder-Mead",
        X=X,y=ly,k1=k.epa)$par                                

#Backfitting using local linear smoothers, the Epanechnikov kernel and the leave-one-out bandwidth
g.hat=cbe.explicit.bivar(X=X, y=ly, h1 = h0, 
                         k1=k.epa,f.sm = ll.sm,details=TRUE)  
z.hat=ly-g.hat$estimate                             #residuals
p.mm=sqrt(3/(2*mean(z.hat^2)))                      #p_hat   (estimate of p)
f.as=exp(3/(2*p.mm)-g.hat$estimate)                 #f_hat   (estimate of f)

####################### FIGURES ######################################
######################################################################

###### 1st step: first step estimation of component functions ########
par(mfrow = c(1, 2),mar=c(4,2,1,1),mgp=c(2.5,.7,0))
x=X
nx <- ncol(x)
tam=1
ymin=min(c(g.hat[[1]],g.hat[[2]]))
ymax=max(c(g.hat[[1]],g.hat[[2]]))
for (j in 1:nx) {
  m_est=g.hat[[j]]
  x1 <- x[, j]
  if(j==1){
    plot(
      x = sort(x1),
      y = m_est[order(x1)],
      type = "l",
      ylab = "",lwd=2.5,
      font=6,cex.lab=tam+0.1,cex.axis=tam,font.lab=6,cex=tam,
      ylim = c(ymin, ymax),
      xlab = "Veterinary/cow"
    )
    legend("topright",legend=expression(paste("Estimated ",g[1](X[1]))),
           col=1,lty=1,bty = 'n',text.font = 6,lwd=2.5,
           cex=tam*0.8)
  }else{
    plot(
      x = sort(x1),
      y = m_est[order(x1)],
      type = "l",
      ylab = "",lwd=2.5,
      font=6,cex.lab=tam+0.1,cex.axis=tam,font.lab=6,cex=tam,
      ylim = c(ymin, ymax),
      xlab = "Energy/cow"
    )
    legend("topright",legend=expression(paste("Estimated ",g[2](X[2]))),
           col=1,lty=1,bty = 'n',text.font = 6,lwd=2.5,
           cex=tam*0.8)
  }
  rug(sort(x1))
  
}


###################### Countour lines of the estimated production frontier ########

x1=X[,1]
x2=X[,2]
resu=list(x1,x2,as.numeric(f.as))

par(mfrow = c(1, 1),mar=c(3.5,5,1,1),mgp=c(2,.7,0))
akima.smooth <- with(resu, interp(x=resu[[1]],y=resu[[2]],z=resu[[3]], nx=500, ny=500,
                                  linear = TRUE,extrap = FALSE))
si.zmin <- min(akima.smooth$z,na.rm=TRUE)
si.zmax <- max(akima.smooth$z,na.rm=TRUE)
breaks <- pretty(c(si.zmin,si.zmax),8)
colors <- gray.colors(length(breaks)-1,start = 1,end=0)
filled.contour(akima.smooth, col=colors,levels=breaks, xlab="vet/cow", ylab="energy/cow",
               plot.axes = { axis(1); axis(2);grid();
                 contour(akima.smooth,add = TRUE,levels=breaks,labcex = .9,
                         vfont =c("sans serif","bold"));
                 points(x=resu[[1]],y=resu[[2]], pch = 19, col= 1)})

