
#
# Lecture 5: Modelo de riesgos proporcionales
# Giancarlo Sal y Rosas
# 04/26/17 
#


rm(list=ls(all=TRUE))

install.packages("glmnet")
install.packages("rms")
install.packages("pec")

library(survival)
library(glmnet)
library(rms)
library(pec)

uis <- read.csv(file.choose())

#
# Construcción de un modelo
#


model1 <- coxph( Surv(time, censor)~treat, data=uis)
pdf("fig61.pdf")
plot(model1,col=1:2,xlab="Tiempo (dias)")
legend(700,0.6,legend=c("Tratamiento Corto","Tratamiento Largo"),col=1:2,lty=c(1,1))
dev.off()

summary(model1)


#
# Analisis Bivariado 
#

coxph(Surv(time, censor)~age,data=uis)
coxph(Surv(time, censor)~beck,data=uis)
coxph(Surv(time, censor)~factor(hercoc),data=uis)
coxph(Surv(time, censor)~factor(ivhx),data=uis)
coxph(Surv(time, censor)~ndrugtx,data=uis)
coxph(Surv(time, censor)~race,data=uis)
coxph(Surv(time, censor)~treat,data=uis)
coxph(Surv(time, censor)~site,data=uis)
coxph(Surv(time, censor)~los,data=uis)



# Multivariate analysis

model1 <- coxph(Surv(time,censor)~age + beck + factor(hercoc) + factor(ivhx) + ndrugtx 
          + race + treat + site,data=uis)
summary(model1)

# Removemos la variable hercoc

model2     <- coxph(Surv(time,censor)~age + beck + factor(ivhx) + ndrugtx + race + treat + 
          site,data=uis)
pvalue <- 1 - pchisq(2*(model1$loglik[2] - model2$loglik[2]),3)
pvalue

summary(model2)

#
# ivhx
#

uis$nivhx <- ifelse(is.na(uis$ivhx),NA,ifelse(uis$ivhx < 3,0,1))
model3    <- coxph(Surv(time,censor)~age + beck + nivhx + ndrugtx + race + treat + site,data=uis)
pvalue    <- 1 - pchisq(2*(model2$loglik[2] - model3$loglik[2]),1)
pvalue

summary(model3)

#
# Escala continua 
#

uis$cage     <- cut(uis$age,breaks=quantile(uis$age,probs=c(0,0.25,0.5,0.75,1),na.rm=TRUE),include.lowest=TRUE)
uis$cbeck    <- cut(uis$beck,breaks=quantile(uis$beck,probs=c(0,0.25,0.5,0.75,1),na.rm=TRUE),include.lowest=TRUE)
uis$cndrugtx <- cut(uis$ndrugtx,breaks=quantile(uis$ndrugtx,probs=c(0,0.25,0.5,0.75,1),na.rm=TRUE),include.lowest=TRUE)

model4 <- coxph( Surv(time, censor)~factor(cage)+beck+ndrugtx+nivhx+race+treat+site, data=uis)
model5 <- coxph( Surv(time, censor)~age+factor(cbeck)+ndrugtx+nivhx+race+treat+site, data=uis)
model6 <- coxph( Surv(time, censor)~age+beck+factor(cndrugtx)+nivhx+race+treat+site, data=uis)

cbind(c(24,30.5,35.5,47.5),c(0,as.numeric(summary(model4)$coef[1:3,1])))

pdf("fig62.pdf")
x1 <- par(mfrow=c(2,2))
plot(c(24,30.5,35.5,47.5),c(0,as.numeric(summary(model4)$coef[1:3,1])),type="l",xlab="age",ylab="Logaritmo del riesgo")
points(c(24,30.5,35.5,47.5),c(0,as.numeric(summary(model4)$coef[1:3,1])),col=2)
#
plot(c(5,12.5,20,40),c(0,as.numeric(summary(model5)$coef[2:4,1])),type="l",xlab="beck",ylab="Logaritmo del riesgo")
points(c(5,12.5,20,40),c(0,as.numeric(summary(model5)$coef[2:4,1])),col=2)
#
plot(c(0.5,2.5,5.0,23.5),c(0,as.numeric(summary(model6)$coef[3:5,1])),type="l",xlab="ndrugtx",ylab="Logaritmo del riesgo")
points(c(0.5,2.5,5.0,23.5),c(0,as.numeric(summary(model6)$coef[3:5,1])),col=2)
par(x1)
dev.off()
#
# Fractional polynomials
#

# Vamos a tratar de modelar mejor beck y ndrugtx

library(mfp)
# Tenemos que transformar las variables para que no tomen valores 0
uis$fdrugtx  <- 10/(uis$ndrugtx+1) # Esta es una de las transformaciones mas usadas
m8 <- mfp( Surv(time, censor)~ fp(fdrugtx,df=2) + beck+ age +nivhx+race+treat+site, data=uis,family=cox)
summary(m8)
m9 <- mfp( Surv(time, censor)~ fp(fdrugtx,df=4) + beck+ age +nivhx+race+treat+site, data=uis,family=cox)
summary(m9)

# Comparar modelos

pvalue <- 1 - pchisq(2*(m9$loglik[2] - m9$loglik[1]),4)
pvalue

pvalue <- 1 - pchisq(2*(m9$loglik[2] - m3$loglik[1]),3)
pvalue

pvalue <- 1 - pchisq(2*(m9$loglik[2] - m8$loglik[1]),3)
pvalue

#
# Grafico en base a el nuevo modelo de polinomios fraccionales
#

y      <- m9$x[,6]*m9$coef[6] + m9$x[,7]*m9$coef[7]
# Delete missing values
uis$ok <-complete.cases(uis$age, uis$beck, uis$ndrugtx, uis$nivhx, uis$site, uis$race, uis$treat)
to.uis <-uis[uis$ok,]
f51.d  <- cbind(x=to.uis$ndrugtx,y)
f51.d  <- f51.d[order(f51.d[,1]),]
par(cex=.8)
plot(f51.d[,1], f51.d[,2], type="l", xlab="NDRUGTX", ylab="Log riesgo instántaneo")


#
# Interacciones
#

uis$ndrugfp1 <- 1/((uis$ndrugtx+1)/10)
uis$ndrugfp2 <- (1/((uis$ndrugtx+1)/10))*log((uis$ndrugtx+1)/10)

uis$agesite      <- uis$age*uis$site
uis$becksite     <- uis$beck*uis$site
uis$ndrugfp1site <- uis$ndrugfp1*uis$site
uis$ndrugfp2site <- uis$ndrugfp2*uis$site
uis$nivhxsite    <- uis$nivhx*uis$site
uis$racesite    <- uis$race*uis$site

uis$beckage      <- uis$beck*uis$age
uis$ndrugfp1age  <- uis$ndrugfp1*uis$age
uis$ndrugfp2age  <- uis$ndrugfp2*uis$age
uis$raceage      <- uis$race*uis$age

uis$ndrugfp1race  <- uis$ndrugfp1*uis$race
uis$ndrugfp2race  <- uis$ndrugfp2*uis$race
uis$beckrace  <- uis$ndrugfp1*uis$race

m10 <- coxph( Surv(time, censor)~age+beck+ndrugfp1+ndrugfp2+nivhx+race+treat+site+
                agesite + becksite + ndrugfp1site + ndrugfp2site + nivhxsite + racesite + beckage + ndrugfp1age +
                ndrugfp2age + raceage + ndrugfp1race + ndrugfp2race, data=uis,method="breslow", na.action=na.exclude)


m11 <- coxph( Surv(time, censor)~age+beck+ndrugfp1+ndrugfp2+nivhx+race+treat+site+
                agesite +  racesite , data=uis,method="breslow", na.action=na.exclude)
xtable(m11)


#
# Lasso regression
#

uis   <- uis[complete.cases(uis),]
y     <- cbind(time=uis$time,status=uis$censor)
x     <- as.matrix(uis[,-c(1,10:16)])
fit1  <- glmnet(x, y, family = "cox")

pdf("fig64.pdf")
plot(fit1,label=T)
dev.off()

# Parametro adicional
cv.fit <- cv.glmnet(x, y, family="cox", alpha=1)

pdf("fig65.pdf")
plot(cv.fit)
dev.off()

cv.fit$lambda.min
cv.fit$lambda.1se

coef(cv.fit, s = "lambda.min")


ff       <- age+beck + ndrug + factor(ivhx) +race + treat + site + age*site +  race*site
x1       <- model.matrix(~age + beck + ndrugtx + factor(ivhx) +race 
            + treat + site + age*site +  race*site,uis)
cv.fit1  <- cv.glmnet(x1, y, family="cox", alpha=1)

plot(cv.fit1)

pdf("fig64.pdf")
plot(fit1,label=T)
dev.off()


#
# Back and forward in R
#

library(pec)
fit2 <- selectCox(Surv(time,censor)~age+beck+hercoc+nivhx+
                  ndrugtx+race+treat+site+site*race+age*site,data=uis,rule="aic")
fit2$fit

fit3 <- selectCox(Surv(time,censor)~age+beck+hercoc+nivhx+
                    ndrugtx+race+treat+site+site*race+age*site,data=uis,rule="p")
fit3$fit

fit4 <- selectCox(Surv(time,censor)~age+beck+hercoc+nivhx+
                  ndrugtx+race+treat+site+site*race+age*site + 
                  beck*site+nivhx*site+ndrugtx*site,data=uis,rule="p")
fit4$fit

#
# pbc dataset
#



