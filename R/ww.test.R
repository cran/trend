ww.test <- function(x) {
##    Copyright (C) 2016  Thorsten Pohlert
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
##    This function computes the Wald-Wolfowitz test
##    for independence and stationarity.
##
##    Wald, A. and J. Wolfowitz, 1943: An exact test for 
##    randomness in the non-parametric case based on 
##    serial correlation. Annual Mathematical Statistics, 14: 378--388.
##
    DNAME <- deparse(substitute(x))
    stopifnot(is.numeric(x))
    x <- na.omit(x)

    n <- length(x)
    if (n < 2) {
        stop("sample size must be greater than 1")
    }
    
# The test statistic r
    r <- sum(x[1:(n-1)] * x[2:n]) + x[1] * x[n]

    s1 <- sum(x^1)
    s2 <- sum(x^2)
    s3 <- sum(x^3)
    s4 <- sum(x^4)

# expected value of r
    er <- (s1^2 - s2) / (n - 1)

# expected variance of r
    vr <- (s2^2 - s4) / (n - 1) - er^2 +
          (s1^4 - 4 * s1^2 * s2 + 4 * s1 * s3 + s2^2 - 2 * s4) /
          ((n - 1) * (n - 2))

# Standardised quantile (approaches normal distribution)	  
    z <- (r - er) / sqrt(vr)
    names(z) <- "z-value"
    
# p-value, two-sided case
    if (z >= 0) {
        pval <- 2 * pnorm(z, lower.tail = FALSE)
    }
    else {
        pval <- 2 * pnorm(z, lower.tail = TRUE)
    }

    out <- list(statistic = z, p.value= pval,
                alternative = "The series is significantly different from  \n independence and stationarity",
                method = "Wald-Wolfowitz test for independence and stationarity",
                data.name = DNAME, nobs = n)
			
    class(out) <- "htest"
    return(out)
}
