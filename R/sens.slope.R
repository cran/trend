sens.slope <-
function(x, level=0.95){
    na.fail(x)
    res.mk <- mk.test(x)
    varS <- res.mk$Varianz
    n <- length(x)
    d <- NULL
    retval <- list(b.sen = NULL, b.sen.up = NULL, b.sen.lo = NULL,
                   intercept = NULL, nobs = n, method="SSLP", D=NULL,
                   conf.int = level * 100, varS = NULL)
    k <- 0
    for (i in 1:(n-1)) {
        for (j in (i+1):n){
            k <- k + 1
            d[k] <- (x[j] - x[i]) / ( j - i)
        }
    }
    b.sen <- median(d, na.rm=TRUE)
    C <- qnorm(1 - (1 - level)/2) * sqrt(varS)
    rank.up <- round((k + C) / 2 + 1)
    rank.lo <- round((k - C) / 2)
    rank.d <- sort(d)
    retval$b.sen.lo <- rank.d[rank.lo]
    retval$b.sen.up <- rank.d[rank.up]
    intercepts <- x - b.sen * (1:n)
    retval$intercept <- median(intercepts)
    retval$b.sen <- b.sen
    retval$D <- d
    retval$varS <- varS
    class(retval) <- "trend.test"
    return(retval)
}

