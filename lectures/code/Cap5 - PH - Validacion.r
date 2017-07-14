
############################################
#
# Lecture: Model validation
# Date: 052417
#
############################################

rm(list=ls(all=TRUE))

library(survival)
library(foreign)

uis       <- read.csv(file.choose())
uis$ok    <- complete.cases(uis$age, uis$beck, uis$ndrugtx, uis$nivhx, uis$site, uis$race, uis$treat)
uis       <- uis[uis$ok,]
uis$nivhx <- ifelse(is.na(uis$ivhx),NA,ifelse(uis$ivhx < 3,0,1))

model1    <- coxph(Surv(time, censor)~ treat + age+beck+ndrugtx+age*site + race*site,data=uis,method="breslow")
summary(model1)

#########################################################################
#
# Residuos
#
#########################################################################

rM    <- residuals(model1,type="martingale")

##############################
#
# Martingalas
# 
#############################


pdf("fig71.pdf")
plot(uis$age,rM,ylab="Residuos de Martingala",xlab="Edad (Age)")
lines(smooth.spline(uis$age,rM),col="red",lwd=2)
lines(lowess(uis$age,rM),col="blue",lwd=2)
legend("bottomright",legend=c("Splines","Lowess"),col=c("red","blue")
       ,lty=rep(1,1),lwd=c(2,2))
dev.off()

pdf("fig72.pdf")
plot(density(as.numeric(rM)),main="",xlab="Residuos martingalas",ylab="Densidad")
dev.off()

pdf("fig73.pdf")
x1 <- par(mfrow=c(2,2))
plot(uis$age,rM,ylab="Residuos de Martingala",xlab="Edad (Age)")
lines(smooth.spline(uis$age,rM),col="red",lwd=2)
lines(lowess(uis$age,rM),col="blue",lwd=2)
legend("bottomright",legend=c("Splines","Lowess"),col=c("red","blue"),
       lty=rep(1,1),lwd=c(2,2))
#
plot(uis$beck,rM,ylab="Residuos de Martingala",xlab="Depresión (Beck)")
lines(smooth.spline(uis$beck,rM),col="red",lwd=2)
lines(lowess(uis$beck,rM),col="blue",lwd=2)
legend("bottomright",legend=c("Splines","Lowess"),col=c("red","blue"),
       lty=rep(1,1),lwd=c(2,2))
#
plot(uis$ndrugtx,rM,ylab="Residuos de Martingala",
     xlab="# tratamientos previos")
lines(smooth.spline(uis$ndrugtx,rM),col="red",lwd=2)
lines(lowess(uis$ndrugtx,rM),col="blue",lwd=2)
legend("bottomright",legend=c("Splines","Lowess"),col=c("red","blue"),
       lty=rep(1,1),lwd=c(2,2))
par(x1)
dev.off()

########################
#
# Cox-Snell
#
########################

pdf("fig74.pdf")
r_CS    <- uis$censor - residuals(model1,type="martingale")
fitres  <- survfit(coxph(Surv(r_CS,censor)~1,method='breslow',data=uis),type='aalen')
plot(fitres$time,-log(fitres$surv),type='s',xlab="Residuos de Cox-Snell", 
     ylab=expression(paste(Lambda,"(residuos)")))
abline(0,1,col='red',lty=2)
dev.off()

#######################
#
# Deviance residuals
#
########################


pdf("fig75.pdf")
x1 <- par(mfrow=c(1,2))
rM    <- residuals(model1,type="deviance")
lp    <- predict(model1,type="lp")
plot(uis$age,rM,ylab="Residuos de Desviación",xlab="Edad",ylim=c(-3.1,3.5))
abline(h=c(-3,3),col="red")
#
plot(lp,rM,ylab="Residuos de Desviación",xlab="Predictor lineal",ylim=c(-3.1,3.5))
abline(h=c(-3,3),col="red")
par(x1)
dev.off()


########################################################
#
# Score residuals
#
########################################################

rM    <- residuals(model1,type="score")

pdf("fig76.pdf")
plot(uis$age,as.numeric(rM[,1]),xlab="Edad",ylab="Residuos de Score")
abline(h=0,col="red")
dev.off()

pdf("fig77.pdf")
x1 <- par(mfrow=c(2,2))
plot(uis$age,as.numeric(rM[,1]),xlab="Edad",ylab="Residuos de Score")
abline(h=0,col="red")
#
plot(uis$age,as.numeric(rM[,2]),xlab="Depresión",ylab="Residuos de Score")
abline(h=0,col="red")
#
plot(uis$age,as.numeric(rM[,3]),xlab="Número de tratamientos 1",ylab="Residuos de Score")
abline(h=0,col="red")
#
plot(uis$age,as.numeric(rM[,4]),xlab="Número de tratamientos 2",ylab="Residuos de Score")
abline(h=0,col="red")
par(x1)
dev.off()


#
# Delta-beta
#

rdb           <- residuals(model1,type="dfbeta")

# Variación %
per_rdb       <- t(apply(rdb,1,function(x) 100*x/model1$coef))


pdf("fig78.pdf")
boxplot(per_rdb[,1],xlab="Edad",ylab="Diferencia porcentual")
dev.off()

pdf("fig79.pdf")
boxplot(per_rdb,xlab="Variables",ylab="Diferencia porcentual")
abline(h=-10,col=2)
abline(h=10,col=2)
dev.off()

#
# Mirar obsevacione especificas
#

# Observaciones que hacer cambiar el coeficiente 10% o mas
round(per_rdb[rowSums(abs(per_rdb) > 10)>0,],3)
uis[rowSums(abs(per_rdb) > 10)>0,]

pdf("fig710.pdf")
a <- boxplot(uis$beck,xlab="Beck")
text(x=c(1,1),y=a$out,labels=a$out,col=rep(2,2),pos=4)
dev.off()



#
# Estudio de la  proporcionalidad
#

ph <- cox.zph(model1)
ph

# Depresion
plot(ph[2])
abline(h=0, lty=3)

# Tratamiento
plot(ph[7])
abline(h=0, lty=3)


phl <- cox.zph(model1,transform="log")
phl

# Depresion
plot(phl[2])
abline(h=0, lty=3)

# Tratamiento
plot(phl[7])
abline(h=0, lty=3)


