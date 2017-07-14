
#
# Lecture 5: Modelo de riesgos proporcionales
# Giancarlo Sal y Rosas
# 04/26/17 
#


rm(list=ls(all=TRUE))

library(survival)

#
# Interpretation
#

#
# Binary response
#

model1 <- coxph(Surv(time,status)~sex,data=cancer)
summary(model1)

# Intervalo de confianza
confint(model1)
exp(confint(model1))

# Politomica

cancer$ph.ecog <- ifelse(cancer$ph.ecog < 2,cancer$ph.ecog,2)
model2         <- coxph(Surv(time,status)~factor(ph.ecog),data=cancer)
summary(model2)

confint(model2)
exp(confint(model2))


model2$var
c(0.56 - 1.96*sqrt(0.06),0.56 + 1.96*sqrt(0.06))
exp(c(0.56 - 1.96*sqrt(0.06),0.56 + 1.96*sqrt(0.06)))

# Continua
model3 <- coxph(Surv(time,status)~age,data=cancer)
summary(model3)

# Intervalo de confianza
c(0.02*10 - 1.96*10*sqrt(model3$var),0.02 + 1.96*10*sqrt(model3$var))
exp(c(0.02*10 - 1.96*10*sqrt(model3$var),0.02 + 1.96*10*sqrt(model3$var)))

# Multivariado

model4 <- coxph(Surv(time,status)~age + sex + factor(ph.ecog),data=cancer)
summary(model4)

# Comparación de modelos
model5 <- coxph(Surv(time,status)~ sex + factor(ph.ecog),data=cancer)
summary(model5)

anova(model5,model4)


#############################################
#
# Estimation of the baseline hazard function
#
#############################################



# Función de supervivencia

pdf("fig51.pdf")
plot(survfit(model1,newdata=data.frame(sex=1),conf.int=FALSE),
     xlab="Tiempo (dias)",ylab=" Función de supervivencia ",col=1)
lines(survfit(model1,newdata=data.frame(sex=2),conf.int=FALSE),col=2)
legend("topright",legend=c("Hombre","Mujeres"),title="Sexo del paciente",
       lty=c(1,1),col=1:2)
dev.off()


#
# Función basehaz
# 

# Base on the basehaz function
cancer$sex <- as.numeric(cancer$sex) - 1
model1     <- coxph(Surv(time,status)~sex,data=cancer)
haz        <- basehaz(model1,centered=FALSE)
pdf("fig52.pdf")
plot(haz$time,haz$hazard,xlab="Tiempo (dias)",ylab="Riesgo acumulado basal",
     col=1,type="s")
aux <- survfit(model1,newdata=data.frame(sex=0),conf.int=FALSE)
lines(aux$time,-log(aux$surv),col=2)
legend("topleft",legend=c("basehaz","survfit"),title="Método",lty=c(1,1),col=1:2)
dev.off()

#
# Adjusting by age
#

model3 <- coxph(Surv(time,status)~age+sex,data=cancer)
summary(model3)
aux1 <- survfit(model3,newdata=data.frame(sex=1,age=62.5),conf.int=FALSE)
aux2 <- survfit(model3,newdata=data.frame(sex=2,age=62.5),conf.int=FALSE)

pdf("fig53.pdf")
plot(aux1,xlab="Tiempo (dias)",ylab="Función de supervivencia")
lines(aux2$time,aux2$surv,col=2)
legend("topright",legend=c("Hombres","Mujeres"),lty=c(1,1),col=1:2)
dev.off()


#
# Risk score
#

uis <- read.csv(file.choose())

model4 <- coxph(Surv(time,censor)~age+beck+race+treat+site,data=uis)
summary(model4)

lp <- predict(model4,newdata=uis,type="lp")
summary(lp)

pdf("fig54.pdf")
hist(lp,prob=TRUE,xlab="Score",main="")
dev.off()


s25 <- survfit(model4)$surv^exp(quantile(model4$linear.predictors,0.25))
s50 <- survfit(model4)$surv^exp(quantile(model4$linear.predictors,0.5))
s75 <- survfit(model4)$surv^exp(quantile(model4$linear.predictors,0.75))

pdf("fig55.pdf")
plot(survfit(model4)$time,s25,xlab="Tiempo (dias)",
     ylab=" Función de supervivencia ",col=1,type="l",ylim=c(0,1))
lines(survfit(model4)$time,s50,col=2,type="l")
lines(survfit(model4)$time,s75,col=3,type="l")
legend(500,0.8,legend=c("Percentil 25","Percentil 50","Percentil 75"),
       title="Cuantiles del score de riesgo",lty=c(3,3),col=c(1,2,3))
dev.off()

pdf("fig56.pdf")
sa <- survfit(model4)$surv^exp(quantile(model4$linear.predictors,0.5))
sb <- survfit(model4)$surv^exp(quantile(model4$linear.predictors-coef(model4)[4],0.5))
plot(survfit(model4)$time,sa,xlab="Tiempo (dias)",
     ylab=" Función de supervivencia ",col=1,type="l",ylim=c(0,1))
lines(survfit(model4)$time,sb,col=2,type="l")
legend(700,0.8,legend=c("Largo","Corto"),
       title="Tratamiento",lty=c(1,1),col=1:2)
dev.off()

