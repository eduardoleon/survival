
# 
# Análisis de Supervivencia
# Capítulo 1: Giancarlo Sal y Rosas
# 02/16/17
# * Significa que esa parte puede ser omitida

rm(list=ls(all=TRUE))
library(splines)
library(survival)
library(foreign)
library(numDeriv)

############################################################
#
# HPTN 039
#

setwd("C:/Users/Giancarlo/Dropbox/Courses/Temas de Bioestadistica - Survival analysis/Lectures/Lecture 1 - Introduction")
hptn <- read.csv("hptn.csv")

pdf("fig11.pdf")
plot(survfit(Surv(time=t1,time2=t2,event=event)~arm,data=hptn),conf.int=F,fun="event",
     ylim=c(0,0.3),mark.time=F,col=c(1,2),xlab="Tiempo (dias)",ylab="Probabilidad acumulada de VIH")
legend(50,0.2,legend=c("Placebo","Aciclovir"),col=1:2,lty=c(1,1))
dev.off()

############################################################
#
# Partners notification study
#

pns <- read.table("ptrdatafinal.txt",sep=",")

# variables:
#gender:         1=male, 2=female
#randmass:       5=standard, 6=expedited (recoded so #0= standard; 1 = expedited)
#ptid:           study id
#ct_gc:          initial diagnosis: C=chlamydia, #G=gonorrhea, B=both
#ct_res:         ct result at followup
#gc_res:         gc result at followup
#corrctgc:       positive at followup for a disease #patient had at baseline.
#                (this was the primary outcome: 0=no, #1=yes)
#age:            age in years
#folltime:       days between initial and followup #visits.
#follprim:       followup visit within primary #analysis window? (0=no,1=yes)

names(pns) = c("gender","arm","ptid","ct_gc","ct_res","gc_res","age","corrctgc","folltime","follprim")
pns        = pns[!is.na(pns$folltime),]
pns        = pns[!is.na(pns$corrctgc),]
pns$gender = pns$gender - 1
pns$arm    = pns$arm    - 5
pns        = pns[order(pns$folltime,decreasing=F),]
head(pns)

# 
# Derivada del diagrama de sumas acumuladas de (G,V)*
#

auxiliar = function(i,n,V,G)
{
  if(i==0) aux <- V[(i+1):n]*(G[(i+1):n])^(-1)
  else     aux <- (V[(i+1):n]-V[i])*(G[(i+1):n]-G[i])^(-1)
  aux
}

left = function(G,V)
{
  n <- length(V)
  x <- rep(0,n)
  i <- 0
  while(i < n) 
  {
    aux.1          <- auxiliar(i,n,V,G)
    aux.2          <- min(aux.1)
    aux.3          <- max(c((i+1):n)[aux.1 <= aux.2])
    x[(i+1):aux.3] <- aux.2
    i              <- aux.3
  }  
  x
}

pns  <- pns[order(pns$folltime,decreasing=F),]
resi <- left(G=c(1:dim(pns[pns$arm==1,])[1]),V=cumsum(pns$corrctgc[pns$arm==1]))
resc <- left(G=c(1:dim(pns[pns$arm==0,])[1]),V=cumsum(pns$corrctgc[pns$arm==0]))

pdf("fig12.pdf")
plot(c(0,pns$folltime[pns$arm==1]),c(0,resi),type="s",xlab="Tiempo (dias)"
     ,ylab="Probabilidad acumulada de reinfección",ylim=c(0,1))
points(c(0,pns$folltime[pns$arm==0]),c(0,resc),type="s",col=2)
legend(15,0.8,legend=c("Intervención","Control"),col=1:2,lty=c(1,1))
dev.off()

####################################################

# 
# Desercion Universitaria
#

deser <- read.csv("tesis_database.csv")
N     <- 10000
index <- sample(1:dim(deser)[1],N,replace=FALSE)
deser <- deser[index,]
head(deser)

pdf("fig13.pdf")
mod1 <- survfit(Surv(time,censored)~gender,data=deser)
plot(mod1,ylab="Probabilidad acumulada",xlab="Tiempo (ciclos)",
     mark.time=F,ylim=c(0,1),xlim=c(1,20),col=1:2,fun="event")
legend(2,0.8,legend=c("Mujer","Hombre"),col=1:2,lty=c(1,1))
dev.off()


############################
#
# Hazard function
#
############################

pdf("fig13b.pdf")
x1 = par(mfrow=c(1,3))
t <- seq(0,5,by=0.01)
plot(t,dlnorm(t,0.5,1)/(1-plnorm(t,0.5,1)),type="l",lty=1,yaxt="n",xaxt="n",
     xlab="Tiempo",ylab=expression(lambda(t)),
     main="Pacientes con tuberculosis")
plot(t,(0.5/1)*(t/1)^(3-1),type="l",,ylim=c(0,10),yaxt="n",xaxt="n",
     xlab="Tiempo",ylab=expression(lambda(t)),
     main="Pacientes con leucemia")
plot(t,(0.1/1)*(t/2)^(0.1-1),type="l",,ylim=c(0,0.5),yaxt="n",xaxt="n",
     xlab="Tiempo",ylab=expression(lambda(t)),
     main="Pacientes post cirugia")
par(x1)
dev.off()


####################################################
#
# Exponetial model
#

pdf("fig14.pdf")
x1 <- par(mfrow=c(1,3))
t <- seq(0,5,by=0.01)
plot(t,0.5*exp(-0.5*t),xlab="Tiempo",ylab="Función de Densidad",type="l",
     lty=1,ylim=c(0,0.8))
lines(t,1*exp(-1*t),type="l",lty=2)
lines(t,2*exp(-2*t),type="l",lty=3)
#
plot(t,1-exp(-0.5*t),xlab="Tiempo",ylab="Función de Distribución",type="l",
     lty=1,ylim=c(0,1))
lines(t,1-exp(-t),type="l",lty=2)
lines(t,1-exp(-2*t),type="l",lty=3)
legend(2,0.25,legend=expression(paste(eta,"= 0.5"),paste(eta,"= 1"),
      paste(eta,"= 2")),lty=1:3)
#
plot(t,0.5*t/t,type="l",lty=1,ylim=c(0,2.5),xlab="Tiempo",
     ylab="Riesgo instantáneo")
lines(t,1*t/t,type="l",lty=2)
lines(t,2*t/t,type="l",lty=3)
par(x1)
dev.off()

####################################################
#
# Weibull model
#

pdf("fig15.pdf")
x1 <- par(mfrow=c(1,3))
t <- seq(0,5,by=0.01)
#
plot(t,dweibull(t,1,0.5),xlab="Tiempo",ylab="Función de Densidad",type="l",ylim=c(0,1),col=1)
lines(t,dweibull(t,1,1),type="l",col=2)
lines(t,dweibull(t,1,2),type="l",col=3)
lines(t,dweibull(t,0.5,1),type="l",col=4)
lines(t,dweibull(t,0.1,1),type="l",col=5)
#
plot(t,pweibull(t,1,0.5),xlab="Tiempo",ylab="Función de Distribución",type="l",ylim=c(0,1),col=1)
lines(t,pweibull(t,1,1),type="l",col=2)
lines(t,pweibull(t,1,2),type="l",col=3)
lines(t,pweibull(t,0.5,1),type="l",col=4)
lines(t,pweibull(t,0.1,1),type="l",col=5)
legend(1.2,0.3,
       legend=c("(a,b)=(1,0.5)","(a,b)=(1,1)","(a,b)=(1,2)","(a,b)=(0.5,1)","(a,b)=(0.1,1)"),lty=rep(1,5),col=1:5)
#
plot(t,(1/0.5)*(t/1)^(1-1),type="l",lty=1,ylim=c(0,2.5),xlab="Tiempo",ylab="Riesgo instantáneo")
lines(t,(1/1)*(t/1)^(1-1),type="l",col=2)
lines(t,(1/2)*(t/2)^(1-1),type="l",col=3)
lines(t,(0.5/1)*(t/1)^(0.5-1),type="l",col=4)
lines(t,(0.1/1)*(t/2)^(0.1-1),type="l",col=5)
par(x1)
dev.off()


####################################################
#
# Estimation without censoring
#

t <- c(14.2,7.3,21.6,0.5,13.5,3.0,9.7,7.4,3.1,1.9)

# Exponential model
lexp  <- function(t,eta){
  n   <- length(t)  
  res <- n*log(eta) - eta*sum(t,na.rm=T)
}
mexp     <- optimize(lexp,lower=0.05,upper=4,t=t,maximum=TRUE)


# Weibull model
lweibull <- function(x,t)
{
  n   <- length(t)
  a   <- x[1]
  b   <- x[2]
  res <- n*log(a) - n*log(b) + (a-1)*sum(log(t)) - sum(t^a)/(b^a) - (n*a-n)*log(b)
  return(-res)
}  
mwei     <- nlminb(c(1,1),lweibull,t=t)
mwei
ee <- solve(hessian(lweibull,x=mwei$par,t=t))
ee

aes <- mwei$par[1]
bes <- mwei$par[2]

pdf("fig16.pdf")
plot(seq(0,max(t),by=0.1),1-pexp(seq(0,max(t),0.1),mexp$maximum),xlab="Tiempo (meses)",ylab=expression(S[n]),main=" ",lty=1,type="l",col=1)
segments(x0=5,x1=5,y0=0,y1=1-pexp(5,mexp$maximum),lty=2)
segments(x0=-1,x1=5,y0=1-pexp(5,mexp$maximum),y1=1-pexp(5,mexp$maximum),lty=2)
segments(x0=10,x1=10,y0=0,y1=1-pexp(10,mexp$maximum),lty=2)
segments(x0=-1,x1=10,y0=1-pexp(10,mexp$maximum),y1=1-pexp(10,mexp$maximum),lty=2)
dev.off()

pdf("fig17.pdf")
plot(seq(0,max(t),by=0.1),exp(-(seq(0,max(t),0.1)/bes)^aes),xlab="Tiempo (meses)",ylab=expression(S[n]),main=" ",lty=1,type="l",col=1)
segments(x0=5,x1=5,y0=0,y1=exp(-(5/bes)^aes),lty=2)
segments(x0=-1,x1=5,y0=exp(-(5/bes)^aes),y1=exp(-(5/bes)^aes),lty=2)
segments(x0=10,x1=10,y0=0,y1=exp(-(10/bes)^aes),lty=2)
segments(x0=-1,x1=10,y0=exp(-(10/bes)^aes),y1=exp(-(10/bes)^aes),lty=2)
dev.off()

# Nonparametric model

Survival <- function(t)
{
 t <- t[order(t,decreasing=F)]
 n <- length(t)
 S <- apply(array(1:n),1,function(i) 1-mean(t<=t[i]))
 return(list(t=c(0,t),Su=c(1,S)))
}

pdf("fig18.pdf")
aux = Survival(t)
plot(aux$t,aux$Su,type="n",
     xlab="Tiempo (meses)",ylab=expression(S[n]),main=" ")
for(i in 1:(length(aux$t)+1)){
  segments(x0=aux$t[i],y0=aux$Su[i],
           x1=aux$t[i+1],y1=aux$Su[i])
  points(aux$t[i],aux$Su[i],lwd=2)
}
dev.off()


# Nonparametric vs. parametric

pdf("fig19.pdf")
aux = Survival(t)
plot(aux$t,aux$Su,type="n",
     xlab="Tiempo (meses)",ylab=expression(S[n]),main=" ")
for(i in 1:(length(aux$t)+1)){
  segments(x0=aux$t[i],y0=aux$Su[i],
           x1=aux$t[i+1],y1=aux$Su[i])
  points(aux$t[i],aux$Su[i],lwd=2)
}
points(seq(0,max(t),by=0.1),1-pexp(seq(0,max(t),0.1),mexp$maximum),type="l",col=2)
points(seq(0,max(t),by=0.1),exp(-(seq(0,max(t),0.1)/bes)^aes),type="l",col=3)
legend(14,0.7,legend=c("Non parametrica","Exponencial","Weibull"),col=1:3,lty=rep(1,3))
dev.off()

####################################################
# 
# Censoring graphic
#


pdf("fig110.pdf")
x1 <- par(mfrow=c(1,2))
plot(c(0,1),c(0,1),type="n",xlab=" ",ylab=" ",axes=F)
segments(x0=0.1,y0=0.12,x1=0.1,y1=1)
abline(h=0.15)
segments(x0=0.9,y0=0.12,x1=0.9,y1=1)
text(0.11,0.1,"Inicio del")
text(0.1,0.05,"estudio")
text(0.9,0.1,"Fin del")
text(0.9,0.05,"estudio")
text(0.5,0.1,"Tiempo")
text(0.5,0.05,"calendario")
segments(x0=0.3,y0=0.3,x1=0.7,y1=0.3)
points(x=0.3,y=0.3,pch=4)
points(x=0.7,y=0.3,pch=1)
segments(x0=0.1,y0=0.6,x1=0.6,y1=0.6)
points(x=0.1,y=0.6,pch=4)
points(x=0.6,y=0.6,pch=4)
segments(x0=0.2,y0=0.8,x1=0.9,y1=0.8)
segments(x0=0.9,y0=0.8,x1=0.98,y1=0.8,lty=3)
points(x=0.2,y=0.8,pch=4)
points(x=0.98,y=0.8,pch=4)
#
plot(c(0,1),c(0,1),type="n",xlab=" ",ylab=" ",axes=F)
segments(x0=0.1,y0=0.12,x1=0.1,y1=1)
abline(h=0.15)
text(0.1,0.1,"0")
text(0.5,0.1,"Tiempo del")
text(0.5,0.05,"paciente")
segments(x0=0.1,y0=0.3,x1=0.5,y1=0.3)
points(x=0.5,y=0.3,pch=1)
segments(x0=0.1,y0=0.6,x1=0.6,y1=0.6)
points(x=0.6,y=0.6,pch=4)
segments(x0=0.1,y0=0.8,x1=0.88,y1=0.8)
points(x=0.88,y=0.8,pch=1)
par(x1)
dev.off()


#
# Estructura de los datos
#

head(lung)
Surv(time=lung$time,event=lung$status)[1:10]


hptn <- hptn[!is.na(hptn$t1),]
head(hptn)
Surv(time = hptn$t1,time2 = hptn$t2,event = hptn$event)[1:10]

