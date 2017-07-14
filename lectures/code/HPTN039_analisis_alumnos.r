
#####################################################
#
# Análisis: HPTN 039
#
######################################################

rm(list=ls(all=TRUE))

library(survival)
library(foreign)

hptn <- read.dta(file.choose())

aux <- (hptn$ptid == 625030074 & hptn$t2 == 352 & !is.na(hptn$t2))
hptn <- hptn[!(hptn$ptid == 625030074 & hptn$t2 == 352 & !is.na(hptn$t2)),]
hptn <- hptn[!(hptn$ptid == 205001122 & hptn$t2 == 56  & !is.na(hptn$t2)),]
hptn <- hptn[!(hptn$ptid == 629002856 & hptn$t2 == 274 & !is.na(hptn$t2)),]
hptn <- hptn[!is.na(hptn$circum),]


#
# Data descrita en cada visita trimestral
#

# Ciudad
hptn$site <- ifelse(hptn$site==203,"NY",ifelse(hptn$site==204,"SF",
             ifelse(hptn$site==205,"Seattle",ifelse(hptn$site==627,"Pucallpa",
             ifelse(hptn$site==628,"Iquitos","Lima")))))

###############################################################
#
# Descripción de las tres ultimas parejas
#

# Suma de parejas nuevas
hptn$tnew <- rowSums(cbind(hptn$qsmp1new==1,hptn$qsmp2new==1,hptn$qsmp3new==1),na.rm=T)

# Total de actos insertivos y receptivos

hptn$tii  <-  rowSums(cbind(hptn$qsmp1isx,hptn$qsmp2isx,hptn$qsmp3isx),na.rm=T)   			 
hptn$tri  <-  rowSums(cbind(hptn$qsmp1rsx,hptn$qsmp2rsx,hptn$qsmp3rsx),na.rm=T) 				 	 
hptn$ta   <-  hptn$tii + hptn$tri
hptn$pia  <-  ifelse(hptn$ta>0,hptn$tii/(hptn$tii + hptn$tri),0)

########################################################################
#
# Análisis
#

# Regresión simple
plot(survfit(Surv(t1,t2,pos) ~ factor(circum), data=hptn,
     conf.int = FALSE),col=1:2,fun="event",ylim=c(0,0.4),
     xlab="Tiempo (Dias)",ylab="Función acumulada de probabilidad")
legend(10,0.34,legend=c("No","Si"),title="Circuncisión",col=1:2,lty=rep(1,2))
#

# Efecto de circuncisión
model1 <- coxph(Surv(t1,t2,pos)~ factor(circum)*pia + 
          factor(circum)*anyp3 + strata(site),data=hptn)
summary(model1)


