perturb = function(df, n1, n2, nc, dx) {
    s1 = deparse(substitute(n1))
    s2 = deparse(substitute(n2))
    sc = deparse(substitute(nc))
    d = 0
    c = df[1,sc]
    for (i in 1:nrow(df)) {
        if (c == df[i,sc]) {
            df[i,s1] = df[i,s1] + d
        } else {
            d = 0
            c = df[i,sc]
        }
        
        if (df[i,s1] < df[i,s2]) {
            d = 0
        } else {
            d = d + dx
            df[i,s2] = df[i,s2] + d
        }
    }
    df
}
