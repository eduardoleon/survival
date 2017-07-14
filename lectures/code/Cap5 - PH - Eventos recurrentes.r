
#
# Eventos recurrentes
# 06/02/17
#


rm(list=ls(all=TRUE))
library(survival)


# Estudio: Tumor de vejiga

head(bladder2)

#
# Proceso de conteo
#

model1 <- coxph(Surv(start,stop,event) ~ factor(rx) + number + size ,
               data=bladder2)
summary(model1)

model2 <- coxph(Surv(start,stop,event) ~ factor(rx) + number + size ,
                data=bladder2,robust = TRUE)
summary(model2)

#
# Proceso condicional
#

bladder[bladder$id==14,]
#
bladderB       <- bladder2
bladderB$start <- 0
bladderB$stop <- bladder2$stop - bladder2$start
bladderB[bladderB$id==14,]

#
bladder2[bladder2$id==14,]

model3 <- coxph(Surv(stop,event) ~ factor(rx) + number + size + strata(enum),
                data=bladder)
summary(model3)

H0 <- basehaz(model3,centered=FALSE)
pdf("fig91.pdf")
plot(H0$time[as.numeric(H0$strata)==1],H0$hazard[as.numeric(H0$strata)==1],
     type="s",ylab="Función de riesgo basal",xlab="Tiempo (dias)")
for(i in 2:4) lines(H0$time[as.numeric(H0$strata)==i],H0$hazard[as.numeric(H0$strata)==i],
                    type="s",col=i)
legend("topleft",legend=c("1","2","3","4"),col=1:4,lty=rep(1,4),title="# de recurrencias")
dev.off()

model4 <- coxph(Surv(start,stop,event) ~ factor(rx) + number + size + strata(enum),
                data=bladder2)
summary(model4)

pdf("fig92.pdf")
H0b <- basehaz(model4,centered=FALSE)
plot(H0b$time[as.numeric(H0b$strata)==1],H0b$hazard[as.numeric(H0b$strata)==1],
     type="s",ylab="Función de riesgo basal",xlab="Tiempo (dias)",ylim=c(0,6))
for(i in 2:4) lines(H0b$time[as.numeric(H0b$strata)==i],H0b$hazard[as.numeric(H0b$strata)==i],
                    type="s",col=i)
legend("topleft",legend=c("1","2","3","4"),col=1:4,lty=rep(1,4),title="# de recurrencias")
dev.off()

#
# Condicional B
#

model5 <- coxph(Surv(stop,event) ~ factor(rx) + number + size + strata(enum),
                data=bladderB,robust = TRUE)
summary(model5)


