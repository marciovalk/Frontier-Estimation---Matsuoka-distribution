## Generating samples from Matsuoka distribution

###############################################
# Brutforce version of 
# Matsuoka's cdf
# Matsuoka's qdf
# Matsuoka's rdf
###################################
mdf<-function(x,p){
  2*sqrt(-p^3*log(x)/pi)*x^(p-1)
}
mcdf<-function(x,p) integrate(mdf,p=p,0,x)$value
mqdf<-function(x,p) optimize(function(z)(mcdf(z,p=p)-x)^2,c(0,1))$minimum
mrdf<-function(n,p) sapply(runif(n),p=p,mqdf)