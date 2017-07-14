survival = function(m, t) {
    n = findInterval(t, m$time)
    c(lower = m$lower[n], surv = m$surv[n], upper = m$upper[n])
}

hazard = function(m, t) {
    n = findInterval(t, m$time)
    -log(c(lower = m$upper[n], hazard = m$surv[n], upper = m$lower[n]))
}

rmean = function(m, g = NULL) {
    # Work with the whole data set
    if (is.null(g)) {
        t = m$time
        s = m$surv
        n = length(t)
    }
    
    # Select the specified stratum
    else {
        n = m$strata[g]
        d = sum(m$strata[1:(g-1)])
        r = d + 1:n
        t = m$time[r]
        s = m$surv[r]
    }
    
    # Generate consecutive points of the curve
    w = diff(c(1,s))
    q = rev(cumsum(rev(w*t)) / cumsum(rev(w)))
    m = list()
    m$x = c(0, rep(t, each = 2, len = 2*n - 1))
    m$y = rep(q, each = 2) - m$x
    class(m) = "rmean"
    m
}

plot.rmean = function(m, ...) plot(m$x, m$y, type = "l", ...)
