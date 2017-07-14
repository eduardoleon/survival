library("survival")
library("pec")
source("code/general.r")

# Example: Worcester Heart Attack Study
whas = read.table("hosmer/2ed/whas500.dat")

# This function does all the work
selectCox(Surv(V21, V22) ~ V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V11 + V13 + V14, data = whas, rule = "aic")
selectCox(Surv(V21, V22) ~ V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V11 + V13 + V14, data = whas, rule = "p")
