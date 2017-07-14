
#####################################################
#
# Time dependent covariates
#
######################################################

rm(list=ls(all=TRUE))

library(survival)
library(foreign)

#
# HPTN 039
#

hptn <- read.dta(file.choose())
hptn  <- subset(hptn,select=c("visit","ptid","t1","t2","age","circum","qsmp1isx","qsmp1rsx",
        "qsmp2isx","qsmp2rsx","qsmp3isx","qsmp3rsx","site","qsmp1new","qsmp2new","qsmp3new",
        "pos","arm"))


# Any new  partner
hptn$anyn <- ifelse((hptn$qsmp1new ==1 & !is.na(hptn$qsmp1new)) | (hptn$qsmp2new==1 & !is.na(hptn$qsmp2new)) | (hptn$qsmp3new==1 & !is.na(hptn$qsmp3new)),1,0)


hptn          <- hptn[!is.na(hptn$t1),]
hptn$ site    <- ifelse(hptn$site==203,"NY",ifelse(hptn$site==204,"SF",ifelse(hptn$site==205
                  ,"Seattle",ifelse(hptn$site==627,"Pucallpa",ifelse(hptn$site==682,"Iquitos","Lima")))))

hptn$qsmp1isx <- ifelse(is.na(hptn$qsmp1isx),0,hptn$qsmp1isx)
hptn$qsmp2isx <- ifelse(is.na(hptn$qsmp2isx),0,hptn$qsmp2isx)
hptn$qsmp3isx <- ifelse(is.na(hptn$qsmp3isx),0,hptn$qsmp3isx)

hptn$qsmp1rsx <- ifelse(is.na(hptn$qsmp1rsx),0,hptn$qsmp1rsx)
hptn$qsmp2rsx <- ifelse(is.na(hptn$qsmp2rsx),0,hptn$qsmp2rsx)
hptn$qsmp3rsx <- ifelse(is.na(hptn$qsmp3rsx),0,hptn$qsmp3rsx)

hptn$tii <-  hptn$qsmp1isx +  hptn$qsmp2isx +  hptn$qsmp3isx   			 
hptn$tri <-  hptn$qsmp1rsx +  hptn$qsmp2rsx +  hptn$qsmp3rsx 				 	 

hptn$ta  <- hptn$tii + hptn$tri
hptn$pia <- ifelse(hptn$ta>0,hptn$tii/(hptn$tii + hptn$tri),0)

# Any partner in the last three months
hptn$anyp3 <- ifelse(hptn$ta > 0,1,0)

hptn[1:5,c(2:3,5:13)]
table(hptn$site)

pdf("fig81.pdf")
plot(survfit(Surv(t1,t2,pos) ~ factor(arm), data=hptn,
     conf.int = FALSE),ylim=c(0,0.5),col=1:2,fun="event",
     xlab="Tiempo (Dias)",ylab="Función acumulada de probabilidad")
legend("topleft",legend=c("Placebo","Aciclovir"),col=1:2,lty=rep(1,2))
dev.off()


model1 <- coxph(Surv(t1,t2,pos)~ arm,data=hptn)
summary(model1)

model2 <- coxph(Surv(t1,t2,pos)~ arm + strata(site),data=hptn)
summary(model2)


#
# Graphic of the baseline functions
#
haz     <- basehaz(model2,centered=FALSE)
haz$s   <- as.numeric(haz$strata)

pdf("fig82.pdf")
plot(haz$time,haz$hazard,type="n",xlab="Tiempo (dias)",ylab="Función acumulada de riesgo")
for(i in 1:5) points(haz$time[haz$s==i],haz$hazard[haz$s==i],type="l",col=i)
legend("topleft",legend=names(table(haz$strata)),lty=rep(1,5),col=1:5)
dev.off()

###################################################
#
# Graphic of the hazard ratio
#
###################################################

hptn[1:25,c("ptid","t1","t2","tii","ta","pia","circum")]

model3 <- coxph(Surv(t1,t2,pos)~ circum*pia + strata(site),data=hptn)
summary(model3)


pdf("fig83.pdf")
lp <- exp(coef(model3)[1] + coef(model3)[3]*seq(0,1,by=0.01))
plot(seq(0,1,0.01),lp,xlab="Proporción de actos insertivos",ylab="Cociente de riesgos",type="l")
abline(h=1,lty=2)
dev.off()
