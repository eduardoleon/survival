library("survival")
library("mfp")
source("code/general.r")

# Example: Worcester Heart Attack Study
whas = read.table("hosmer/2ed/whas500.dat")


################################################################################
# Step 1: Select individually significant covariates with p = 0.2

# Discrete covariate selection: V10, V12 discarded due to insufficient data
survdiff(Surv(V21, V22) ~ V3, data = whas)    # Gender
survdiff(Surv(V21, V22) ~ V8, data = whas)    # History of Cardiovascular Disease
survdiff(Surv(V21, V22) ~ V9, data = whas)    # Atrial Fibrillation
survdiff(Surv(V21, V22) ~ V10, data = whas)   # Cardiogenic Shock
survdiff(Surv(V21, V22) ~ V11, data = whas)   # Congestive Heart Complications
survdiff(Surv(V21, V22) ~ V12, data = whas)   # Complete Heart Block
survdiff(Surv(V21, V22) ~ V13, data = whas)   # MI Order
survdiff(Surv(V21, V22) ~ V14, data = whas)   # MI Type

# Continuous covariate selection: No covariates discarded
coxph(Surv(V21, V22) ~ V2, data = whas)       # Age
coxph(Surv(V21, V22) ~ V4, data = whas)       # Initial Heart Rate
coxph(Surv(V21, V22) ~ V5, data = whas)       # Initial Systolic Blood Pressure
coxph(Surv(V21, V22) ~ V6, data = whas)       # Initial Dyastolic Blood Pressure
coxph(Surv(V21, V22) ~ V7, data = whas)       # Body Mass Index


################################################################################
# Step 2:  Build an initial combined linaer model

m1 = coxph(Surv(V21, V22) ~ V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V11 + V13 + V14, data = whas)
summary(m1)


################################################################################
# Step 3:  Discard non-significant covariantes in the combined model

changes = function(m1, m2) {
    m1 = m1$coef
    m2 = m2$coef
    sapply(names(m2), function(n) (m2[n] - m1[n]) / m1[n])
}

# V5, V8 discarded, p > 0.9, changes under 20%
m2 = coxph(Surv(V21, V22) ~ V2 + V3 + V4 + V6 + V7 + V9 + V11 + V13 + V14, data = whas)
summary(m2)
changes(m1, m2)

# V13 discarded, p > 0.7, changes under 20%
m3 = coxph(Surv(V21, V22) ~ V2 + V3 + V4 + V6 + V7 + V9 + V11 + V14, data = whas)
summary(m3)
changes(m2, m3)

# V9, V14 discarded, p > 0.3, changes under 20%
m4 = coxph(Surv(V21, V22) ~ V2 + V3 + V4 + V6 + V7 + V11, data = whas)
summary(m4)
changes(m3, m4)

# V3 is of clinical importance, can't be discarded


################################################################################
# Step 4:  Confirm that discarded variables aren't statistically significant

# omitted


################################################################################
# Step 5:  Examine the scale of the continuous covariates

# Generate quartile midpoints
mids = c(1,3,5,7) / 8
x1 = quantile(whas$V2, mids)
x2 = quantile(whas$V4, mids)
x3 = quantile(whas$V6, mids)
x4 = quantile(whas$V7, mids)

# Rerun the model after quantizing each continuous covariate
fq = function(x) factor(quantize(x))
y1 = c(0, coxph(Surv(V21, V22) ~ fq(V2) + V3 + V4 + V6 + V7 + V11, data = whas)$coef[1:3])
y2 = c(0, coxph(Surv(V21, V22) ~ V2 + V3 + fq(V4) + V6 + V7 + V11, data = whas)$coef[3:5])
y3 = c(0, coxph(Surv(V21, V22) ~ V2 + V3 + V4 + fq(V6) + V7 + V11, data = whas)$coef[4:6])
y4 = c(0, coxph(Surv(V21, V22) ~ V2 + V3 + V4 + V6 + fq(V7) + V11, data = whas)$coef[5:7])

jpeg("dump/p.5.jpg")
p = par(mfrow = c(2,2))
plot(x = x1, y = y1, type = "l", ylab = "Log-hazard", xlab = "Age")
plot(x = x2, y = y2, type = "l", ylab = "Log-hazard", xlab = "Initial Heart Rate")
plot(x = x3, y = y3, type = "l", ylab = "Log-hazard", xlab = "Diastolic Blood Pressure")
plot(x = x4, y = y4, type = "l", ylab = "Log-hazard", xlab = "Body Mass Index")
par(p)
dev.off()

# Only V7 has nonlinear respnose
summary(mfp(Surv(V21, V22) ~ fp(V2) + V3 + fp(V4) + fp(V6) + fp(V7) + V11, data = whas, family = cox))

# Main effects model
V7a = (whas$V7 / 10) ^ 2
V7b = (whas$V7 / 10) ^ 3
m5 = coxph(Surv(V21, V22) ~ V2 + V3 + V4 + V6 + V7a + V7b + V11, data = whas)
summary(m5)


################################################################################
# Step 6:  Select interaction terms

# Evaluate age and gender interactions
anova(m5, coxph(Surv(V21, V22) ~ V2 + V3 + V4 + V6 + V7a + V7b + V11 + V2 * V3, data = whas))
anova(m5, coxph(Surv(V21, V22) ~ V2 + V3 + V4 + V6 + V7a + V7b + V11 + V2 * V4, data = whas))
anova(m5, coxph(Surv(V21, V22) ~ V2 + V3 + V4 + V6 + V7a + V7b + V11 + V2 * V6, data = whas))
anova(m5, coxph(Surv(V21, V22) ~ V2 + V3 + V4 + V6 + V7a + V7b + V11 + V2 * V7a + V2 * V7b, data = whas))
anova(m5, coxph(Surv(V21, V22) ~ V2 + V3 + V4 + V6 + V7a + V7b + V11 + V2 * V11, data = whas))
anova(m5, coxph(Surv(V21, V22) ~ V2 + V3 + V4 + V6 + V7a + V7b + V11 + V3 * V4, data = whas))
anova(m5, coxph(Surv(V21, V22) ~ V2 + V3 + V4 + V6 + V7a + V7b + V11 + V3 * V6, data = whas))
anova(m5, coxph(Surv(V21, V22) ~ V2 + V3 + V4 + V6 + V7a + V7b + V11 + V3 * V7a + V3 * V7b, data = whas))
anova(m5, coxph(Surv(V21, V22) ~ V2 + V3 + V4 + V6 + V7a + V7b + V11 + V3 * V11, data = whas))

# Preliminary model
m6 = coxph(Surv(V21, V22) ~ V2 + V3 + V4 + V6 + V7a + V7b + V11 + V2 * V3 + V2 * V6 + V2 * V11 + V3 * V6, data = whas)
summary(m6)
anova(m5, m6)


################################################################################
# Step 7:  Validate the model

# Compute martingale residuals
mr = residuals(m6, type = "martingale")
martingale = function(z, xlab) {
    plot(z, mr, xlab = xlab, ylab = "Martingale residuals")
    lines(smooth.spline(z, mr), col = "red")
    lines(lowess(z, mr), col = "blue")
    legend("bottomright", legend = c("Splines", "Lowess"), col = c("red", "blue"), lty = 1)
}

# Plot martingale residuals
jpeg("dump/p.7.mr.jpg")
p = par(mfrow = c(2,2))
martingale(whas$V2, xlab = "Age")
martingale(whas$V4, xlab = "Initial Heart Rate")
martingale(whas$V6, xlab = "Diastolic Blood Pressure")
martingale(whas$V7, xlab = "Body Mass Index")
par(p)
dev.off()

# Cox-Snell residuals
jpeg("dump/p.7.cr.jpg")
rm = survfit(coxph(Surv(whas$V22 - mr, whas$V22) ~ 1, method = "breslow"), type = "aalen")
plot(rm, fun = "cumhaz", conf.int = FALSE)
abline(0, 1, col = "red")
dev.off()

# Compute deviance residuals
dr = residuals(m6, type = "deviance")
deviance = function(z, xlab) {
    plot(z, dr, xlab = xlab, ylab = "Deviance residuals")
    abline(h = c(-3, 3), col = "red")
}

# Plot deviance residuals
jpeg("dump/p.7.dr.jpg")
p = par(mfrow = c(2,3))
deviance(whas$V2, xlab = "Age")
deviance(whas$V4, xlab = "Initial Heart Rate")
deviance(whas$V6, xlab = "Diastolic Blood Pressure")
deviance(whas$V7, xlab = "Body Mass Index")
deviance(predict(m6, type = "lp"), xlab = "Linear Predictor")
par(p)
dev.off()

# Compute score residuals
sr = residuals(m6, type = "score")
score = function(z, i, xlab) {
    plot(z, as.numeric(sr[,i]), xlab = xlab, ylab = "Score residuals")
    abline(h = 0, col = "red")
}

# Plot score residuals
jpeg("dump/p.7.sr.jpg")
p = par(mfrow = c(2,2))
score(whas$V2, "V2", xlab = "Age")
score(whas$V4, "V4", xlab = "Initial Heart Rate")
score(whas$V6, "V6", xlab = "Diastolic Blood Pressure")
par(p)
dev.off()

# Compute Schoenfeld residuals
ph = cox.zph(m6)
schoenfeld = function(i) {
    plot(ph[i])
    abline(h = 0, col = "red")
}

# Plot Schoenfeld residuals
jpeg("dump/p.7.ph.jpg")
p = par(mfrow = c(2,2))
schoenfeld(1)
schoenfeld(3)
schoenfeld(4)
par(p)
dev.off()


################################################################################
# Exercises

m7 = coxph(Surv(V21, V22) ~ V2 + V3 + V4 + V6 + V7a + V7b + V11 + V2 * V3 + V2 * V6 + V2 * V11 + V3 * V6 + strata(V15), data = whas)
changes(m6, m7)
