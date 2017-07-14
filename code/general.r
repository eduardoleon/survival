interval = function(p) c(1-p, 1, 1+p) / 2

quantize = function(x, ...) cut(x, breaks = quantile(x, na.rm = TRUE, ...), include.lowest = TRUE)

new.plot = function(xlim, xlab, ylim, ylab, leg, col) {
    plot(0, xlim = xlim, xlab = xlab, ylim = ylim, ylab = ylab, type = "n")
    legend("topright", legend = leg, col = col, lty = 1)
}

quality = function(m) {
    aic = extractAIC(m)
    par = c(aic[1], m$loglik[2], aic[2])
    names(par) = c("df", "loglik", "aic")
    par
}

changes = function(m1, m2) {
    m1 = m1$coef
    m2 = m2$coef
    sapply(names(m2), function(n) (m2[n] - m1[n]) / m1[n])
}

count.valid = function(df, ...) {
    p = rep(0, nrow(df))
    for (col in c(...)) p = p + !is.na(df[col])
    p
}
