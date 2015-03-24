csmk.test <-
function(x){
# Correlated Seasonal Mann Kendall Test
    if (!is.ts(x))
        stop("error: input must be ts object")
    fr <- frequency(x)
    if(fr < 2)
        stop("error: minimum of two seasons required.")
    ts.tsp <- tsp(x)
    n <- length(abs(ts.tsp[1]):abs(ts.tsp[2]))
    y <- matrix(data=NA, nrow=n, ncol=fr)
    for (i in 1:fr) {
        y[,i] <- x[cycle(x) == i]
    }
    res <- multivar.MK.test(y, method="CSMK")
    class(res) <- "trend.test"
    res
}

