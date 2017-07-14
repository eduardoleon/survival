library("survival")
source("code/general.r")

# Pregunta 1
km = survfit(Surv(time, status) ~ sex, data = cancer)
jpeg("dump/3.1.jpg")
plot(km, col = c("blue", "red"), xlab = "Tiempo de vida", ylab = "Función de supervivencia")
legend("topright", c("Hombres", "Mujeres"), col = c("blue", "red"), lty = 1)
dev.off()
survdiff(Surv(time, status) ~ sex, data = cancer)

# Pregunta 2
survdiff(Surv(time, status) ~ quantize(age), data = cancer)

# Pregunta 3
survdiff(Surv(time, status) ~ ph.ecog, data = cancer)
