sea.sens.slope <-
    function(x){
        na.fail(x)
        if (!is.ts(x)){
            stop("error: input must be ts object")
        }
        fr <- frequency(x)
        if (fr < 2){
            stop("error: minimum of two seasons required.")
        }

        n <- length(x)
        retval <- list(b.sen = NULL,
                   intercept = NULL, nobs = n, method="SeaSLP")
        y <- NULL
        d <- NULL
        varS <- 0
        for (i in 1:fr) {
            y <- x[cycle(x) == i]
            res <- sens.slope(ts(y))
            D.m <- res$D
            d <- c(d, D.m)
            varS <- varS + varS
        }

        b.sen <- median(d, na.rm=TRUE)
        intercepts <- as.vector(x) - b.sen * (1:n)
        retval$intercept <- median(intercepts)
        retval$b.sen <- b.sen
        class(retval) <- "trend.test"
        return(retval)

    }
