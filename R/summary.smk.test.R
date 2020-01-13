##   file summary.smk.test.R part of package trend
##   Copyright (C) 2015-2018 Thorsten Pohlert
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
#' @title Object summaries
#' @rdname summary.smktest
#' @description
#' Generic function "summary" for objects of class \code{smktest}. 
#' @param object an object of class \code{smktest}
#' @param \ldots further arguments, currently ignored
#' @importFrom stats symnum
#' @export
summary.smktest <- function(object, ...)
{ 
    if(!inherits(object, "smktest")){
        stop("'object' must be of class 'smktest'")
    }
    cat(paste0("\n\t", object$method, "\n\n"))
    cat(paste0("data: ", object$data.name, "\n"))
    cat(paste0("alternative hypothesis: ", object$alternative, "\n\n"))
    cat(paste0("Statistics for individual seasons\n\n"))
    
    if (object$alternative == "two.sided"){
        pform <- "Pr(>|z|)"
        H0 <- "S = 0"
    } else if (object$alternative == "greater"){
        pform <- "Pr(>z)"
        H0 <- "S <= 0"
    } else {
        pform <- "Pr(<z)"
        H0 <- "S >= 0"
    }
    
    ##   Symbols
    symp <- symnum(object$pvalg, corr=FALSE,
                   cutpoints = c(0,  .001,.01,.05, .1, 1),
                   symbols = c("***","**","*","."," "))

    da <- data.frame(S = object$Sg,
                     varS = round(object$varSg, 1),
                     tau = round(object$taug, 3),
                     z = round(object$Zg, 3),
                     p = format.pval(object$pvalg),
                     symp)
    names(da) <- c("S", "varS", "tau", "z", pform, "")
    
    rownames(da) <- sapply(1:length(object$Zg),
                           function(i) paste0("Season ", i, ":   ",H0))
    cat("H0\n")
    print(da)
    cat("---\n")
    cat("Signif. codes: ", attr(symp, 'legend'),"\n\n")
    invisible(object)
}

