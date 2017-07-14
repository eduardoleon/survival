
###########################
#
# Lecture 3: Prueba de hipótesis
# Giancarlo Sal  y Rosas
# Date: 03/11/15
#
##########################

rm(list=ls(all=TRUE))

library(splines)
library(survival)


#############################
#
# HPTN 039: Analyses
#
#############################

hiv <- read.csv(file.choose())
hiv <- hiv[!is.na(hiv$t1),]
hiv <- hiv[!is.na(hiv$t2),]

pdf("fig31.pdf")
plot(survfit(Surv(time=t1,time2=t2,event=event)~arm,data=hiv),conf.int=F,fun="event",
     ylim=c(0,0.3),mark.time=F,col=c(1,2),xlab="Tiempo (dias)",ylab="Probabilidad acumulada de VIH")
legend(50,0.2,legend=c("Placebo","Aciclovir"),col=1:2,lty=c(1,1))
dev.off()


#############################################
#
# Cervical Cancer
#
#############################################

cancer <- read.csv(file.choose())

pdf("fig32.pdf")
plot(survfit(Surv(time=time,event=rec)~Group,data=cancer),conf.int=F,
     ylim=c(0,1),mark.time=F,col=c(1,2),xlab="Tiempo (años)",
     ylab="Probabilidad de supervivencia")
legend(0.5,0.25,legend=c("Abandono","Completo"),col=1:2,lty=c(1,1))
dev.off()

summary(survfit(Surv(time=time,event=dead)~Group,data=cancer))

library(survival)
survdiff(Surv(time=time,event=dead)~factor(Stage),data=cancer)

plot(survfit(Surv(time=time,event=dead)~factor(Stage),data=cancer))


###########################################
#
# Analysis
#
###########################################


mod1 <- summary(survfit(Surv(futime,fustat)~rx,data=ovarian))
mod1

# Cancer de ovario

mod1 <- survdiff(Surv(futime,fustat)~rx,data=ovarian)
mod1
1-pchisq(mod1$chisq,1)

# Recurrente event: Cancer
survdiff(Surv(time=time,event=rec)~Group,data=cancer)

# Time to dead: Cancer
survdiff(Surv(time=time,event=dead)~Group,data=cancer)


#
# Time to HIV infection: HPTN dataset
#

survdiff(Surv(time=time,event=cen)~arm,data=hiv)

#
# Generalized log rank test 
#

hiv$nage <- cut(hiv$age,breaks=c(18-1,24,30,38,max(hiv$age)))
survdiff(Surv(time=time,event=cen)~factor(nage),data=hiv)

pdf(paste(dirfile,"hptn039_age.pdf",sep=""))
plot(survfit(Surv(time=time,event=cen)~factor(nage),data=hiv),conf.int=F,fun="event",
     ylim=c(0,0.25),mark.time=F,lty=1:4,xlab="Tiempo (dias)",ylab="Probabilidad acumulada de VIH")
legend(50,0.2,legend=c("18-24","25-30","31-38","39-79"),title="Edad",lty=1:4)
dev.off()


#
# Stratified log rank test
#

# Categorization of the age variable
cancer$nage <- cut(cancer$Age,breaks=c(20,40,60,max(cancer$Age)))
survdiff(Surv(time=time,event=rec)~Group+strata(nage),data=cancer)

# Effect of Age stratified by country
survdiff(Surv(time=time,event=cen)~factor(nage)+strata(country),data=hiv)

# Ponderado
survdiff(Surv(time=time,event=rec)~Group,rho=0,data=cancer)
survdiff(Surv(time=time,event=rec)~Group,rho=1,data=cancer)

# Grafico de riesgos

pdf("fig34.pdf")

par(mfrow=c(1,3))

t<-seq(0,6,0.1)

h1<-1/4*t+6
h2<-(-1/4)*t+4
plot(t,h1,ylim=c(2,8),type="l",xlab="Tiempo",ylab="Función de riesgo",col="red")
points(t,h2,type="l")

h3<-1/4*t+3
h4<-(-1/4)*t+7
plot(t,h3,ylim=c(2,8),type="l",xlab="Tiempo",ylab="Función de riesgo",col="red")
points(t,h4,type="l")

h5<-1/2*t+3.5
h6<-(-1/2)*t+6.5
plot(t,h5,ylim=c(2,8),type="l",xlab="Tiempo",ylab="Función de riesgo",col="red")
points(t,h6,type="l")

dev.off()


#
#  
#

uis <- read.csv(file.choose())


