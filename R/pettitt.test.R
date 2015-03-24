pettitt.test <- function(x){
    na.fail(x)
    n <- length(x)
   
    DNAME <- deparse(substitute(x))
    retval <- list(nobs = n, 
                   statistic = NULL, estimate=NULL,
                   p.value = NULL,
                   data.name=DNAME,
                   alternative="true change point is present in the series",
                   method = "Pettitt's test for single change-point detection")
    res <- .Fortran(F_pettitt, n=as.integer(n),
                x=as.single(x), pval=as.single(0),
                tau=as.integer(0), Kt=as.integer(0))              
    K <- res$Kt
    tau <- res$tau
    pval <- res$pval
    names(K) <- "K"
    retval$statistic <- K
    names(tau) <- "probable change point at tau"
    retval$estimate <- tau
    retval$p.value <-  pval
    class(retval) <- "htest"
    return(retval)
}
