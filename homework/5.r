library("survival")
source("code/general.r")

# 1.a, 1.c
cm = coxph(Surv(time, status) ~ sex + age + ph.ecog, data = cancer)
summary(cm)

# 1.b
bhaz = basehaz(cm, FALSE)
resp1 = exp(c(1,cm$means[2],cm$means[3]) %*% cm$coef)[1]
resp2 = exp(c(2,cm$means[2],cm$means[3]) %*% cm$coef)[1]
jpeg("dump/5.1.jpg")
new.plot(xlim = c(0,1022), xlab = "Tiempo de vida",
         ylim = c(0,1), ylab = "Función de supervivencia",
         leg = c("Hombres", "Mujeres"), col = c("blue", "red"))
lines(bhaz$time, exp(-bhaz$hazard * resp1), type = "s", col = "blue")
lines(bhaz$time, exp(-bhaz$hazard * resp2), type = "s", col = "red")
dev.off()

# 1.d
bhaz51 = bhaz$hazard[findInterval(51, bhaz$time)]
bhaz51 * exp(c(1,63,cm$means[3]) %*% cm$coef)[1]
bhaz51 * exp(c(2,63,em$means[3]) %*% cm$coef)[1]

# 2.a
uis = read.csv("data/uis.csv")
uis = with(uis, uis[complete.cases(race, age, treat, ndrugtx),])
cm = coxph(Surv(time, censor) ~ race + age + treat + ndrugtx, data = uis)
summary(cm)

# 2.b
bhaz = basehaz(cm, FALSE)
resp0 = exp(c(cm$means[1],cm$means[2],0,cm$means[4]) %*% cm$coef)[1]
resp1 = exp(c(cm$means[1],cm$means[2],1,cm$means[4]) %*% cm$coef)[1]
jpeg("dump/5.2.jpg")
new.plot(xlim = c(0,1172), xlab = "Tiempo a recaída en drogas",
         ylim = c(0,1), ylab = "Función de supervivencia",
         leg = c("Tratam. corto", "Tratam. largo"), col = c("red", "blue"))
lines(bhaz$time, exp(-bhaz$hazard * resp0), type = "s", col = "red")
lines(bhaz$time, exp(-bhaz$hazard * resp1), type = "s", col = "blue")
dev.off()

# 2.c
d = c(0,5,1,0)
m = d %*% cm$coef
s = sqrt(d %*% cm$var %*% d)
exp(m)
qlnorm(interval(0.95), m, s)
plnorm(1, m, s)

# 2.d
dm = coxph(Surv(time, censor) ~ age + treat + ndrugtx, data = uis)
anova(cm,dm)
