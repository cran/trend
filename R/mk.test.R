mk.test <-
function (x)
{
    if (!is.ts(x))
        stop("error: input must be ts object")
    y <- data.frame(x)
    res <- multivar.MK.test(as.matrix(y), method = "SMK")
    class(res) <- "trend.test"
    res$method <- "MK"
    res
}
