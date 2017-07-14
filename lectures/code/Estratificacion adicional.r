
# Site 0
li0 <- function(x) exp(-0.034095*x - 1.96*sqrt(c(0,0,0,0,x,0,0)%*%model3$var%*%c(0,0,0,0,x,0,0)))
ls0 <- function(x) exp(-0.034095*x + 1.96*sqrt(c(0,0,0,0,x,0,0)%*%model3$var%*%c(0,0,0,0,x,0,0)))

# Site 1
li1 <- function(x) exp(-0.004*x - 1.96*sqrt(c(0,0,0,0,x,0,x)%*%model3$var%*%c(0,0,0,0,x,0,x)))
ls1 <- function(x) exp(-0.004*x + 1.96*sqrt(c(0,0,0,0,x,0,x)%*%model3$var%*%c(0,0,0,0,x,0,x)))

x <- seq(0,11,1)
vli0 <- apply(as.array(x),1,function(i) li0(i))
vls0 <- apply(as.array(x),1,function(i) ls0(i))
vli1 <- apply(as.array(x),1,function(i) li1(i))
vls1 <- apply(as.array(x),1,function(i) ls1(i))


plot(c(0,11),c(0.2,1.4),type="n",xlab="Diferencia en edad",ylab="95%IC del HR")
points(x,vli0,col=1,type="l")
points(x,vls0,col=1,type="l")
points(x,vli1,col=2,type="l",lty=2)
points(x,vls1,col=2,type="l",lty=2)
abline(h=1)
