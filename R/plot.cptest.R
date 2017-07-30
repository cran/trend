##    file plot.cptest.R part of package trend
##
##    Copyright (C) 2017 Thorsten Pohlert
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
##   
#     Changes
#     2017-07-01
#     - initial writing
#
#' @title Plotting cptest-objects
#' 
#' @description
#' Plotting method for objects inheriting from class "cptest"
#'
#' @param x an object of class "cptest"
#' @param \ldots further arguments, currently ignored
#'
#' @examples
#' data(Nile)
#' (out <- br.test(Nile))
#' par(mfrow=c(2,1))
#' plot(Nile) ; plot(out)
#'
#' @export
plot.cptest <- function (x, ...)
{
    OK <- inherits(x, c("cptest"))
    if (!OK)
        stop ("'x' is not an object of class 'cptest'")

    ylab <- attr(x$data, 'nm')
    inp <- list(x = x$data,
                ylab = ylab,
                main = x$method,
                type = "l")

    if (is.ts(x$data)){
        do.call("plot.ts", inp)
    } else {
        do.call("plot.default", inp)
    }

    ## 95 \\% Confidence 
    ##if(grepl("Buishand U", x$method)){
    ##    n <- length(x$data)
    ##    xx <- c(10, 20, 30, 40, 50, 100, Inf)
    ##    uu <- c(0.416, 0.440, 0.447, 0.451, 0.453, 0.457, 0.461)
    ##    u <- approx(xx, uu, xout=n)$y
    ##    u <- u * n
    ##    abline(h = u, lty=2)     
    ##    legend("topleft", legend = c("U", "0.95-CI"), lty=c(1,2))
    ##}
}
