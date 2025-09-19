
lambda <- function(t){
  t
}

q <- rexp

Q <- function(t){
  t^2/2
}

m <- exp(3/8)

ar <- function(){
  while(T){
    s <- q(1)
    integrated_rate <- Q(s)
    if (runif(1)<exp(-integrated_rate+s)/(m)){
      return(integrated_rate)
    }
  }
}

invert_method <- function(ss){
  s <- rexp(ss)
  return(sqrt(2*s))
}

library(ggplot2)
ss <- 1e5
ea <- invert_method(ss)
aro <- replicate(ss,ar())

p <- ggplot()+
  geom_density(aes(x=ea))+
  geom_density(aes(x=aro),col='red')

print(p)


