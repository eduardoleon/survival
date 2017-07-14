library("survival")
source("code/empirical.r")
source("code/general.r")

# Pregunta 1
km = survfit(Surv(futime, fustat) ~ 1, data = jasa)
survival(km, 500)
hazard(km, 51)
quantile(km)
km
quantile(survfit(Surv(futime, fustat) ~ transplant, data = jasa))

# Pregunta 2
km = survfit(Surv(time, status) ~ 1, data = cancer)
survival(km, 500)
hazard(km, 51)
quantile(km)
km
qnorm(c(0.025, 0.975), 376.3, 19.7)
quantile(survfit(Surv(time, status) ~ sex, data = cancer))

# Pregunta 4
km = survfit(Surv(futime, fustat) ~ transplant, data = jasa)
jpeg("dump/2.4.jpg")
new.plot(xlim = c(0,1500), xlab = "Tiempo de vida",
         ylim = c(0,700), ylab = "Media residual",
         leg = c("Sin transplante", "Con transplante"),
         col = c("red", "blue"))
lines(rmean(km, 1), col = "red")
lines(rmean(km, 2), col = "blue")
dev.off()

# Pregunta 5
t = rweibull(300, 0.5, 0.5)
y = runif(300, 0, 5)
x = pmin(t, y)
d = t <= y
km = survfit(Surv(x, d) ~ 1)
ep = sum(d) / sum(x)
wp = list(shape = 0.5, scale = 0.5)
jpeg("dump/2.5.jpg")
new.plot(xlim = c(0,5), xlab = "Tiempo de vida",
         ylim = c(0,1), ylab = "Función de supervivencia",
         leg = c("Kaplan-Meier", "Exponencial", "Weibull"),
         col = c("red", "green", "blue"))
lines(km, col = "red", conf.int = FALSE)
curve(pexp(x, ep, FALSE), col = "green", add = TRUE)
curve(pweibull(x, 0.5, 0.5, FALSE), col = "blue", add = TRUE)
dev.off()

# Pregunta 6
d = sum(cancer$status - 1)
x = sum(cancer$time)
b = interval(0.95)
# 6.a
l = d/x
l
qlnorm(b, log(l), sqrt(1/d))
# 6.c
(1+d) / (1+x)
qgamma(b, 1 + d, 1 + x)

# Pregunta 8
pns = read.csv("data/pns.csv")
left = with(data = pns, ifelse(censor, NA, time))
right = with(data = pns, ifelse(censor, time, NA))
em = survreg(Surv(left, right, type = "interval2") ~ 1, dist = "exp")
summary(em)
exp(-em$coef)
