\name{partial.cor.trend.test}
\Rdversion{1.1}
\alias{partial.cor.trend.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
 Partial correlation trend test
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Performs a partial correlation trend test with either the Pearson's or
  the Spearman's correlation coefficients (\eqn{r(tx.z)}).
}
\usage{
partial.cor.trend.test(x, z, method = c("pearson", "spearman"))
}
\arguments{
  \item{x}{
    a "vector" or "ts" object that contains the variable, which is tested
    for trend (i.e. correlated with time)
  }
  \item{z}{
    a "vector" or "ts" object that contains the variable, which effect
    on "x" is partialled out
  }

  \item{method}{
    a character string indicating which correlation coefficient
    is to be computed. One of "pearson" (default) or "spearman",
    can be abbreviated.
  }
}
\details{
  This function conducts a partial correlation trend test using either
  the "pearson" correlation coefficient, or the "spearman" rank
  correlation coefficient (Hipel and McLoed (2005), p. 882).
  The partial correlation coefficient
  for the response variable "x" with time "t",
  when the effect of the explanatory variable "z" is partialled out,
  is defined as:
  \deqn{r_{tx.z} = \frac{r_{tx} - r_{tz}~r_{xz}}{\sqrt{1 - r_{tz}^2} ~ \sqrt{1-r_{xz}^2}}}{%
    r_{tx.z} = (r_{tx} - r_{tz}~r_{xz}) / (sqrt(1 - r_{tz}^2) sqrt(1-r_{xz}^2))
  }

  The H0: \eqn{r_{tx.z} = 0}{r(tx.z) = 0} (i.e. no trend for "x", when
  effect of "z" is partialled out) is tested against the
  alternate Hypothesis, that there is a trend for "x", when the effect of
  "z" is partialled out.

  The partial correlation coefficient is tested for significance with
  the student t distribution on \eqn{df = n -2} degree of freedom.

}

\value{
  An object of class "htest"

  \item{method}{
    a character string indicating the chosen test
  }
  \item{data.name}{
    a character string giving the name(s) of the data
  }
  \item{cor}{
    the correlation-matrix with t, x and z
  }
  \item{statistic}{
    	the value of the test statistic
  }
  \item{estimate}{
    the partial correlation coefficient \eqn{r(tx.z)}
  }
  \item{parameter}{
    the degrees of freedom of the test statistic in the case that it follows a t distribution
  }
  \item{alternative}{
    a character string describing the alternative hypothesis
  }
  \item{p.value}{
    the p-value of the test
  }
}


\references{
  Hipel, K.W. and McLeod, A.I., (2005).
  Time Series Modelling of Water Resources and Environmental Systems.
  \url{http://www.stats.uwo.ca/faculty/aim/1994Book/}.

 Bahrenberg, G., Giese, E. and Nipper, J., (1992): Statistische Methoden
 in der Geographie, Band 2 Multivariate Statistik, Teubner, Stuttgart.
}
\author{

  T. Pohlert
}

\section{Warning}{
  Current Version is for complete observations only.
}
\seealso{
  \code{\link{cor}},
  \code{\link{cor.test}},
  \code{\link[psych]{partial.r}},
  \code{\link{partial.mk.test}},
}
\examples{
data(maxau)
a <- tsp(maxau) ; tt <- a[1]:a[2]
s <- maxau[,"s"] ; Q <- maxau[,"Q"]
maxau.df <- data.frame(Year = tt, s =s, Q = Q)
plot(maxau.df)

partial.cor.trend.test(s,Q, method="pearson")
partial.cor.trend.test(s,Q, method="spearman")
}

\keyword{ts}
\keyword{nonparametric}
\keyword{multivariate}

