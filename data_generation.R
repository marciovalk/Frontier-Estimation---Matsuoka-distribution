## Code used to simulate data from Matsuoka's distribution

##################################
# density plot Example
require(ggplot2)
p=1
x<-mrdf(2000,p=p)
dataf<- data.frame(x=x)

ggplot(dataf, aes(x = x)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  geom_density(lwd = 1, colour = 4,
               fill = 4, alpha = 0.25)