library("FAdist")
library("haven")
library("survival")
source("code/general.r")

# Leer la base de datos
salud = read_sav("data/ensusalud.sav")
llegada = as.numeric(salud$C1P8_H) * 60 + as.numeric(salud$C1P8_M)
atencion = salud$C1P14
institucion = salud$INSTITUCION
idioma = salud$C1P80

# Pregunta 3
em = survreg(Surv(llegada) ~ 1, dist = "exp")
wm = survreg(Surv(llegada) ~ 1, dist = "wei")
km = survfit(Surv(llegada) ~ 1)
exp(em$coef)
exp(wm$coef)
1 / wm$scale
sink("dump/1.3.txt")
summary(km)
sink()
jpeg("dump/1.3.jpg")
plot(km, xmax = 200, xlab = "Tiempo de llegada", ylab = "Función de supervivencia")
dev.off()

# Pregunta 4
em = survreg(Surv(atencion) ~ 1, dist = "exp")
wm = survreg(Surv(atencion) ~ 1, dist = "wei")
km = survfit(Surv(atencion) ~ 1)
exp(em$coef)
exp(wm$coef)
1 / wm$scale
sink("dump/1.4.txt")
summary(km)
sink()
jpeg("dump/1.4.jpg")
plot(km, xlab = "Tiempo de atención", ylab = "Función de supervivencia")
dev.off()

# Pregunta 5
lm = survreg(Surv(atencion) ~ 1, dist = "loglog")
exp(lm$coef)
1 / lm$scale
pllog(10, lm$scale, lm$coef, FALSE)
pllog(10, lm$scale, lm$coef, FALSE) - pllog(13, lm$scale, lm$coef, FALSE)

# Pregunta 6
quality(em)
quality(wm)
quality(lm)

# Pregunta 7
wm = survreg(Surv(atencion) ~ factor(institucion), dist = "wei")
shape = 1 / wm$scale
scale1 = exp(c(1,0,0,0) %*% wm$coef)
scale2 = exp(c(1,1,0,0) %*% wm$coef)
scale3 = exp(c(1,0,1,0) %*% wm$coef)
scale4 = exp(c(1,0,0,1) %*% wm$coef)
surv = function(x, scale) pweibull(x, shape, scale, FALSE)
jpeg("dump/1.7.jpg")
new.plot(xlim = c(0,90), xlab = "Tiempo de atención",
         ylim = c(0,1), ylab = "Función de supervivencia",
         leg = c("MINSA", "ESSALUD", "FF.AA. - PNP", "Clínicas"),
         col = c("red", "green", "blue", "black"))
curve(surv(x, scale1), 0, 90, col = "red", add = TRUE)
curve(surv(x, scale2), 0, 90, col = "green", add = TRUE, lty = 2)
curve(surv(x, scale3), 0, 90, col = "blue", add = TRUE)
curve(surv(x, scale4), 0, 90, col = "black", add = TRUE)
dev.off()

# Pregunta 8
lm = survreg(Surv(atencion) ~ factor(idioma), dist = "loglog")
scale1 = c(1,0,0) %*% lm$coef
scale2 = c(1,1,0) %*% lm$coef
scale3 = c(1,0,1) %*% lm$coef
haz = function(x, scale) dllog(x, lm$scale, scale) / pllog(x, lm$scale, scale, FALSE)
jpeg("dump/1.8.jpg")
new.plot(xlim = c(0,90), xlab = "Tiempo de atención",
         ylim = c(0,0.3), ylab = "Función de riesgo",
         leg = c("Español", "Quechua", "Otros"),
         col = c("red", "green", "blue"))
curve(haz(x, scale1), 0, 90, col = "red", add = TRUE)
curve(haz(x, scale2), 0, 90, col = "green", lty = 2, add = TRUE)
curve(haz(x, scale3), 0, 90, col = "blue", add = TRUE)
dev.off()
