library("survival")
library("glmnet")
source("code/general.r")

# Example: Worcester Heart Attack Study
whas = read.table("hosmer/2ed/whas500.dat")

x = with(whas, model.matrix(~ V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V11 + V13 + V14))
y = with(whas, Surv(V21, V22))
model = cv.glmnet(x, y, family = "cox")
coef(model, s = "lambda.min")
