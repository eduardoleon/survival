
#
# Lecture 4: Parametric model
# Giancarlo Sal y Rosas
# 04/06/15 
#


rm(list=ls(all=TRUE))

library(splines)
library(survival)


pdf("fig41.pdf")
t  <- seq(0,2,0.01)
plot(t,exp(t)-1,type="l",
     xlab=expression(beta[1]),ylab=expression(exp(beta[1])-1))
lines(t,t,col=2)
dev.off()



pdf("fig42.pdf")
curve(exp(x-exp(x)),-3,3,ylab=expression(paste("f(",epsilon,")")),xlab=expression(epsilon))
dev.off()

#############################################
#
# Exponential regression
#
#############################################


#
# Cervical cancer
#

cancer <- read.csv(file.choose())

# Vamos a reasignar los estadios en estadios I y II                   
cancer$stage <- ifelse(cancer$Stage=="IIA" | cancer$Stage=="IIB","II","I")

model1 <- survreg(Surv(time,dead)~factor(Group),dist="exp",data=cancer)
summary(model1)

pdf("fig43.pdf")
l1  <- exp(-coef(model1)[1])
l2  <- exp(-coef(model1)[1]-coef(model1)[2])
t   <- seq(0,max(cancer$time),0.01)
plot(survfit(Surv(time,dead)~factor(Group),data=cancer),
     ylab="Función de supervivencia",xlab="Tiempo(años)",col=1:2)
lines(t,exp(-l1*t),type="l")
lines(t,exp(-l2*t),type="l",lty=2,col=2)
legend("topright",legend=c("Radioterapia","Histerectomia"),lty=rep(1,2),col=1:2)
dev.off()



model2 <- survreg(Surv(time,dead)~Age,dist="exp",data=cancer)
summary(model2)

# Graphic of hazard functión 
pdf("fig44.pdf")
l1  <- exp(-coef(model2)[1]-coef(model2)[2]*quantile(cancer$Age,0.5))
l2  <- exp(-coef(model2)[1]-coef(model2)[2]*quantile(cancer$Age,0.25))
l3  <- exp(-coef(model2)[1]-coef(model2)[2]*quantile(cancer$Age,0.75))
t   <- seq(0,max(cancer$time),0.01)
plot(t,exp(-l1*t),ylim=c(0,1),type="l",ylab="Función de supervivencia",xlab="Tiempo(años)")
points(t,exp(-l2*t),type="l",col=2,lty=2)
points(t,exp(-l3*t),type="l",col=3,lty=3)
edad <- quantile(cancer$Age,c(0.5,0.25,0.75))
legend("topright",legend=edad,lty=1:3,title="Edad",col=1:3)
dev.off()

model3 <- survreg(Surv(time,dead)~factor(Group)+Age,dist="exp",data=cancer)
summary(model3)

# Graphic of hazard functión 
l1  <- exp(-coef(model4)[1]-coef(model4)[3]*mean(cancer$Age))
l2  <- exp(-coef(model4)[1]-coef(model4)[3]*quantile(cancer$Age,0.25))
l3  <- exp(-coef(model4)[1]-coef(model4)[3]*quantile(cancer$Age,0.75))
t   <- seq(0,4,0.01)

pdf("fig45.pdf")
plot(t,exp(-l1*t),ylim=c(0,1),type="l",ylab="Función de supervivencia",xlab="Tiempo(años)")
points(t,exp(-l2*t),type="l",lty=2)
points(t,exp(-l3*t),type="l",lty=3)
edad <- round(c(mean(cancer$Age),quantile(cancer$Age,c(0.25,0.75))),0)
legend(3,0.9,legend=edad,lty=1:3,title="Edad")
dev.off()

#
# HIV: HPTN 039
#

hiv <- read.csv(paste(dirdata,"hiv.csv",sep=""))
n <- dim(hiv)[1]
npartner <- rowSums(hiv[,9:11],na.rm=T)

model1 <- survreg(Surv(time,cen)~arm+age+circum+factor(edu)+factor(country)+npartner,dist="exp",data=hiv)
summary(model1)

model1 <- survreg(Surv(time,cen)~factor(edu),dist="exp",data=hiv)
summary(model1)

# Graphic of hazard functión 
l1  <- exp(-coef(model1)[1]-coef(model1)[3]*mean(hiv$age))
l2  <- exp(-coef(model1)[1]-coef(model1)[3]*quantile(hiv$age,0.25))
l3  <- exp(-coef(model1)[1]-coef(model1)[3]*quantile(hiv$age,0.75))
t   <- seq(0,500,0.01)
pdf("fig46.pdf")
plot(t,exp(-l1*t),ylim=c(0.6,1),type="l",ylab="Función de supervivencia",xlab="Tiempo(años)")
points(t,exp(-l2*t),type="l",lty=2)
points(t,exp(-l3*t),type="l",lty=3)
edad <- round(c(mean(hiv$age),quantile(hiv$age,c(0.25,0.75))),0)
legend(3,0.9,legend=edad,lty=1:3,title="Edad (años)")
dev.off()

#############################################
#
# Weibull regression
#
#############################################

fun <- function(t,s) s^(-1)*exp(t/s-exp(t/s))
t   <- seq(-5,5,0.01)

pdf("fig46.pdf")
plot(t,fun(t,s=0.5),ylab="Función de densidad",type="l",xlab=expression(epsilon),lty=1,lwd=2)
points(t,fun(t,s=1.5),type="l",lty=2,lwd=2)
points(t,fun(t,s=1),type="l",lty=3,lwd=2)
legend(1.5,0.6,legend=c(expression(paste(sigma,"= 0.5"),paste(sigma,"= 1.5"),
                               paste(sigma,"= 1.0"))),lty=c(1,2,3),lwd=rep(2,3))
dev.off()


#
model1 <- survfit(Surv(time,status)~ sex,data=cancer)
pdf("fig47.pdf")
x1 <- par(mfrow=c(1,2))
plot(model1,xlab="Tiempo (dias)",ylab="Función de supervivencia",
     col=1:2)
legend("topright",legend=c("Hombres","Mujeres"),col=1:2,lty=rep(1,2))
plot(model1,fun="cloglog",xlab="Tiempo (dias)",
     ylab="log(-log(S(t)))",col=1:2)
legend("topleft",legend=c("Hombres","Mujeres"),col=1:2,lty=rep(1,2))
par(x1)
dev.off()
#

model.weibull.1 <- survreg(Surv(time,status)~ sex,dist="w",data=cancer)
summary(model.weibull.1)

pdf("fig48.pdf")
t <- seq(0,3,0.01)
a <- 1/model.weibull.1$scale
plot(t,a*exp(-a*coef(model.weibull.1)[1])*t^(0.32),type="l",col=1,
     xlab="Tiempo (dias)",ylab="Riesgo instantaneo")
lines(t,a*exp(-a*sum(coef(model.weibull.1)))*t^(0.32),
      type="l",col=2)
legend("topleft",legend=c("Hombres","Mujeres"),col=1:2,lty=rep(1,2))
dev.off()


model.weibull.2 <- survreg(Surv(time,status)~ sex +age,dist="w",data=cancer)
summary(model.weibull.2)

#############################################
#
# Log normal regression
#
#############################################

fun <- function(t,s) dnorm(log(t)/s)*(1-pnorm(log(t)/s))^(-1)*(s*t)^(-1)
t   <- seq(0,5,0.01)

pdf(paste(dirfile,"lognormal.pdf",sep=""))
plot(t,fun(t,s=0.5),ylab="Función de riesgo",type="l",lty=1,lwd=2)
points(t,fun(t,s=1.5),type="l",lty=2,lwd=2)
points(t,fun(t,s=1),type="l",lty=3,lwd=2)
legend(3,1,legend=c(expression(paste(sigma,"= 0.5"),paste(sigma,"= 1.5"),
                               paste(sigma,"= 1.0"))),lty=c(1,2,3),lwd=rep(2,3))
dev.off()

model5 <- survreg(Surv(time,dead)~factor(Group)+Age,dist="weibull",data=cancer)
summary(model5)

#############################################
#
# Log normal regression
#
#############################################

fun <- function(t,s) dnorm(log(t)/s)*(1-pnorm(log(t)/s))^(-1)*(s*t)^(-1)
t   <- seq(0,5,0.01)

pdf(paste(dirfile,"lognormal.pdf",sep=""))
plot(t,fun(t,s=0.5),ylab="Función de riesgo",type="l",lty=1,lwd=2)
points(t,fun(t,s=1.5),type="l",lty=2,lwd=2)
points(t,fun(t,s=1),type="l",lty=3,lwd=2)
legend(3,1,legend=c(expression(paste(sigma,"= 0.5"),paste(sigma,"= 1.5"),
                               paste(sigma,"= 1.0"))),lty=c(1,2,3),lwd=rep(2,3))
dev.off()

#############################################
#
# Log logistic regression
#
#############################################


fun <- function(t) exp(t)*(1+exp(t))^(-2)
t   <- seq(-5,5,0.01)

pdf(paste(dirfile,"loglogit.pdf",sep=""))
plot(t,fun(t),ylab="Función de distribución",type="l",lty=1,lwd=2,ylim=c(0,0.6))
points(t,dnorm(t),type="l",lty=2,lwd=1)
legend(2,0.5,legend=c("Log logistico","Normal"),lty=c(1,2),lwd=rep(2,1))
dev.off()
