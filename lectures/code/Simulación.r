

uis <- read.csv(file.choose())

# Sin censura
uis.sc    <- uis[uis$censor==1,]
alpha.hat <- 1/mean(uis.sc$time^2)
alpha.hat

Ssc <- function(t) exp(-alpha.hat*t^2) 
li  <- function(t) Ssc(t) - 1.96*t^2*Ssc(t)*alpha.hat^2/dim(uis.sc)[1]
lf  <- function(t) Ssc(t) + 1.96*t^2*Ssc(t)*alpha.hat^2/dim(uis.sc)[1]

x   <- seq(0,700,by=1)
plot(x,Ssc(x),type="l")
lines(x,li(x),type="l",lty=2)
lines(x,lf(x),type="l",lty=2)

# Censura

alpha.hat.f <- mean(uis$censor)/mean(uis$time^2)
Sc          <- function(t) exp(-alpha.hat.f*t^2) 
lic  <- function(t) Sc(t) - 1.96*t^2*Sc(t)*alpha.hat.f^2/sum(uis.sc$censor)
lfc  <- function(t) Sc(t) + 1.96*t^2*Sc(t)*alpha.hat.f^2/sum(uis.sc$censor)

x   <- seq(0,700,by=1)
plot(x,Ssc(x),type="l")
lines(x,li(x),type="l",lty=2)
lines(x,lf(x),type="l",lty=2)
#
#
lines(x,Sc(x),type="l",col=2)
lines(x,lic(x),type="l",lty=2,col=2)
lines(x,lfc(x),type="l",lty=2,col=2)


#
# Simular data
#

library(glmnet)

N <- 1000
res <- matrix(0,N,10)
for(i in 1:1000)
{
X <- cbind(rbinom(N,1,0.5),rnorm(N,1.7,0.1),runif(N),
           rbinom(N,1,0.5),rbinom(N,1,0.5),rbinom(N,1,0.5),
           rbinom(N,1,0.5),rbinom(N,1,0.5),rbinom(N,1,0.5),
           rbinom(N,1,0.5))
beta   <- c(0.1,-0.5,0.4,rep(0,7))
tim    <- rexp(N,rate=exp(X%*%as.matrix(beta)))
y      <- runif(N,0,5)
censor <- as.numeric(tim<=y)
out    <- cbind(time=tim,status=censor)
cv.model1 <- cv.glmnet(X,out,family="cox",alpha=1)
res[i,] <- as.matrix(coef(cv.model1, s = "lambda.min"))
}








