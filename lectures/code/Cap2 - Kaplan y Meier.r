
###########################
#
# Lecture 2: Estimación de función de supervivencia
# Giancarlo Sal  y Rosas
# Date: 03/11/15
#
##########################

rm(list=ls(all=TRUE))
library(splines)
library(survival)
library(xtable)

#
# Right censoring
#

# Example: Ovarian cancer
head(ovarian)

dat <- ovarian[ovarian$fustat==1,]
dat <- dat[order(dat$futime,decreasing=F),]
dat$futime
mean(dat$futime>353)

ovarian <- ovarian[order(ovarian$futime,decreasing=F),]
ovarian$futime
mean(ovarian$futime>353)


model0 <- summary(survfit(Surv(futime,fustat)~1,data=ovarian))
cols   <- lapply(c(2:4) , function(x) model0[x])
tbl    <- do.call(data.frame, cols)

Surv(ovarian$futime,ovarian$fustat)
tbl

cols   <- lapply(c(2:4,6) , function(x) model0[x])
tbl    <- do.call(data.frame, cols)
tbl

cols   <- lapply(c(2:6,8) , function(x) model0[x])
tbl    <- do.call(data.frame, cols)
tbl

cols   <- lapply(c(2:6,8:10) , function(x) model0[x])
tbl    <- do.call(data.frame, cols)
tbl


km <- survfit(Surv(futime,fustat)~1,data=ovarian)
pdf("fig21.pdf")
plot(km,conf.int = F,mark.time=TRUE,xlab="Tiempo (dias)",
     ylab="Función de supervivencia")
dev.off()

# Varianza
model0 <- summary(survfit(Surv(futime,fustat)~1,data=ovarian))
cols   <- lapply(c(2:4,6,8) , function(x) model0[x])
tbl    <- do.call(data.frame, cols)
tbl

# Intervalo de confianza
model0 <- summary(survfit(Surv(futime,fustat)~1,data=ovarian))
model0

# Diferent type of confindence intervals

summary(survfit(Surv(futime,fustat)~1,data=ovarian,conf.type="plain"))

summary(survfit(Surv(futime,fustat)~1,data=ovarian,conf.type="log"))

summary(survfit(Surv(futime,fustat)~1,data=ovarian,conf.type="log-log"))


#
# km with CI
#

pdf("fig22.pdf")
plot(survfit(Surv(futime,fustat)~1,data=ovarian),xlab="Tiempo (dias)",
     ylab="Función de supervivencia",mark.time=TRUE)
dev.off()

# lung

head(lung)
lung <- lung[order(lung$time,decreasing=F),]
Surv(lung$time,lung$status)[1:60]
summary(survfit(Surv(time,status)~1,data=lung))

pdf("fig23.pdf")
km <- survfit(Surv(time,status)~1,data=lung)
plot(km,mark.time=TRUE,xlab="Tiempo (dias)",
     ylab="Función de supervivencia")
dev.off()


summary(survfit(Surv(futime,fustat)~1,data=ovarian))


#
# Función para calcular mediana e intervalo de confianza
#

per_surv <- function(mod,p=0.5)
{
  qp  <- ifelse(min(mod$surv)  > 1-p,Inf,min(mod$time[mod$surv <= 1-p]))
  lqp <- ifelse(min(mod$lower) > 1-p,Inf,min(mod$time[mod$lower <= 1-p]))
  uqp <- ifelse(min(mod$upper) > 1-p,Inf,min(mod$time[mod$upper <= 1-p]))
  return(list(qp=qp,ic=c(lqp,uqp)))
}

km <- survfit(Surv(futime,fustat)~1,data=ovarian)
per_surv(km,p=0.5)
per_surv(km,p=0.25)


# Median and confidence interval

pdf("fig24.pdf")
km <- survfit(Surv(futime,fustat)~1,data=ovarian)
q  <- per_surv(km,p=0.5)
plot(km,xlab="Tiempo (dias)",
     ylab="Función de supervivencia",xaxt="n")
abline(h=0.5,lty=4,col=2)
segments(x0=q$qp,x1=q$qp,y0=0,y1=0.5,lty=4,col=2)
segments(x0=q$ic[1],x1=q$ic[1],y0=0,y1=0.5,lty=4,col=2)
axis(1,c(q$ic[1],q$qp))
dev.off()


km1 <- survfit(Surv(time,status)~1,data=lung)
per_surv(km1,p=0.5)
per_surv(km1,p=0.25)

# Median and confidence interval: Cáncer de pulmon

pdf("fig25.pdf")
km <- survfit(Surv(time,status)~1,data=lung)
q  <- per_surv(km,p=0.5)
plot(km,xlab="Tiempo (dias)",
     ylab="Función de supervivencia",xaxt="n")
abline(h=0.5,lty=4,col=2)
segments(x0=q$qp,x1=q$qp,y0=0,y1=0.5,lty=4,col=2)
segments(x0=q$ic[1],x1=q$ic[1],y0=0,y1=0.5,lty=4,col=2)
segments(x0=q$ic[2],x1=q$ic[2],y0=0,y1=0.5,lty=4,col=2)
axis(1,c(q$ic[1],q$qp,q$ic[2]))
dev.off()

#
# Media
#

smean <- function(mod,tim){
  Sur <- rep(mod$surv,times=mod$n.risk,mod$n.event)
  res <- sum(Sur*diff(c(0,tim)))
  res
}

km <- survfit(Surv(time,status)~1,data=lung)
smean(km,tim=lung$time)

print(km,rmean="common")
print(km,rmean="individual")
