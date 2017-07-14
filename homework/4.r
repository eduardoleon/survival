library("survival")
source("code/general.r")


################################################################################
# Preguntas 1, 2, 3

block = function(s,a) {
    n = length(a)
    o = rep(1,n)
    s = rep(s,n)
    matrix(c(o,s,a), n, 3)
}

# 1.a, 1.c
em = survreg(Surv(time, status) ~ sex + age, data = cancer, dist = "exp")
summary(em)

# 1.b
jpeg("dump/4.1.b.jpg")
new.plot(xlim = c(0,1022), xlab = "Tiempo de vida",
         ylim = c(0,1), ylab = "Función de supervivencia",
         leg = c("Hombres", "Mujeres"), col = c("blue", "red"))
rate1 = exp(-c(1,1,em$means[3]) %*% em$coef)
rate2 = exp(-c(1,2,em$means[3]) %*% em$coef)
curve(pexp(x, rate1, FALSE), 0, 1022, col = "blue", add = TRUE)
curve(pexp(x, rate2, FALSE), 0, 1022, col = "red", add = TRUE)
dev.off()

# 1.d
jpeg("dump/4.1.d.jpg")
new.plot(xlim = c(39,82), xlab = "Edad",
         ylim = c(0,1022), ylab = "Tiempo de vida",
         leg = c("Hombres", "Mujeres"), col = c("blue", "red"))
curve(exp(block(1,x) %*% em$coef), 39, 82, col = "blue", add = TRUE)
curve(exp(block(2,x) %*% em$coef), 39, 82, col = "red", add = TRUE)
dev.off()

# 1.e
51 / exp(c(1,1,63) %*% em$coef)
51 / exp(c(1,2,63) %*% em$coef)

# 2.a, 2.c
wm = survreg(Surv(time, status) ~ sex + age, data = cancer, dist = "wei")
summary(wm)

# 2.b
jpeg("dump/4.2.b.jpg")
new.plot(xlim = c(0,1022), xlab = "Tiempo de vida",
         ylim = c(0,1), ylab = "Función de supervivencia",
         leg = c("Hombres", "Mujeres"), col = c("blue", "red"))
shape = 1 / wm$scale
scale1 = exp(c(1,1,wm$means[3]) %*% wm$coef)
scale2 = exp(c(1,2,wm$means[3]) %*% wm$coef)
surv = function(x, scale) pweibull(x, shape, scale, FALSE)
curve(surv(x, scale1), 0, 1022, col = "blue", add = TRUE)
curve(surv(x, scale2), 0, 1022, col = "red", add = TRUE)
dev.off()

# 2.d
jpeg("dump/4.2.d.jpg")
new.plot(xlim = c(39,82), xlab = "Edad",
         ylim = c(0,1022), ylab = "Tiempo de vida",
         leg = c("Hombres", "Mujeres"), col = c("blue", "red"))
curve(exp(block(1,x) %*% wm$coef), 39, 82, col = "blue", add = TRUE)
curve(exp(block(2,x) %*% wm$coef), 39, 82, col = "red", add = TRUE)
dev.off()

# 2.e
(51 / exp(c(1,1,63) %*% wm$coef)) ^ shape
(51 / exp(c(1,2,63) %*% wm$coef)) ^ shape

# 3.a, 3.c
lm = survreg(Surv(time, status) ~ sex + age, data = cancer, dist = "loglog")
summary(lm)

# 3.b
jpeg("dump/4.3.b.jpg")
new.plot(xlim = c(0,1022), xlab = "Tiempo de vida",
         ylim = c(0,1), ylab = "Función de supervivencia",
         leg = c("Hombres", "Mujeres"), col = c("blue", "red"))
scale1 = c(1,1,lm$means[3]) %*% lm$coef
scale2 = c(1,2,lm$means[3]) %*% lm$coef
surv = function(x, scale) pllog(x, lm$scale, scale, FALSE)
curve(surv(x, scale1), 0, 1022, col = "blue", add = TRUE)
curve(surv(x, scale2), 0, 1022, col = "red", add = TRUE)
dev.off()

# 3.d
jpeg("dump/4.3.d.jpg")
new.plot(xlim = c(39,82), xlab = "Edad",
         ylim = c(0,1022), ylab = "Tiempo de vida",
         leg = c("Hombres", "Mujeres"), col = c("blue", "red"))
curve(exp(block(1,x) %*% lm$coef), 39, 82, col = "blue", add = TRUE)
curve(exp(block(2,x) %*% lm$coef), 39, 82, col = "red", add = TRUE)
dev.off()

# 3.e
-log(pllog(51, lm$scale, c(1,1,63) %*% lm$coef, FALSE))
-log(pllog(51, lm$scale, c(1,2,63) %*% lm$coef, FALSE))


################################################################################
# Pregunta 4

# 4.a
quality(em)
quality(wm)
quality(lm)


################################################################################
# Pregunta 5

uis = read.csv("data/uis.csv")
em = survreg(Surv(time, censor) ~ treat + race + age, data = uis, dist = "exp")
wm = survreg(Surv(time, censor) ~ treat + race + age, data = uis, dist = "wei")

# 5.a, 5.b, 5.d, 5.e
summary(wm)
exp(c(0,1,0,0) %*% wm$coef)

# 5.c
shape = 1 / wm$scale
scale0 = exp(c(1,0,wm$means[3],wm$means[4]) %*% wm$coef)
scale1 = exp(c(1,1,wm$means[3],wm$means[4]) %*% wm$coef)
surv = function(x, scale) pweibull(x, shape, scale, FALSE)
jpeg("dump/4.5.jpg")
new.plot(xlim = c(0,1172), xlab = "Tiempo a recaída en drogas",
         ylim = c(0,1), ylab = "Función de supervivencia",
         leg = c("Tratam. corto", "Tratam. largo"), col = c("red", "blue"))
curve(surv(x, scale0), 0, 1172, col = "red", add = TRUE)
curve(surv(x, scale1), 0, 1172, col = "blue", add = TRUE)
dev.off()

# 5.f
anova(em,wm)


################################################################################
# Pregunta 6

pns = read.csv("data/pns.csv")
left = ifelse(pns$censor, NA, pns$time)
right = ifelse(pns$censor, pns$time, NA)
em = survreg(Surv(left, right, type = "interval2") ~ pns$z, dist = "exp")
summary(em)
exp(-c(1,0) %*% em$coef)
exp(-c(1,1) %*% em$coef)
