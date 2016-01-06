\name{wm.test}
\Rdversion{1.1}
\alias{wm.test}

\title{
 Wallis and Moore phase-frequency test

}
\description{
  Performes the non-parametric Wallis and Moore phase-frequency
  test for testing the H0-hypothesis, whether the series
  comprises random data, against the HA-Hypothesis, that the series
  is significantly different from randomness (two-sided test).
}
\usage{
wm.test(x)
}

\arguments{
  \item{x}{
    a vector or a time series object of class "ts"
  }
}

\details{
  The test statistic of the phase-frequency test for \eqn{n > 30}
  is calculated as:
  
 \deqn{z = \frac{\| h - \frac{2 n - 7}{3} \|}{\sqrt{\frac{16 n - 29}{90}}}}{%
   z = abs(h - (2 * n - 7) / 3) / sqrt((16 * n - 29)/ 90)}
 

 where \eqn{h} denotes the number of phases, whereas the first and the
 last phase is not accounted. The \eqn{z}-statistic is
 normally distributed. For \eqn{n \le 30} a continuity correction of
 \eqn{-0.5} is included in the denominator.
}

\value{
  An object of class "htest"

  \item{method}{
    a character string indicating the chosen test
  }
  \item{data.name}{
    a character string giving the name(s) of the data
  }
  \item{statistic}{
    the Wallis and Moore z-value
  }
  \item{alternative}{
    a character string describing the alternative hypothesis
  }
  \item{p.value}{
    the p-value for the test
  }
}

\section{Note}{
  NA values are omitted. Many ties in the series will lead to reject H0
  in the present test.
}

\references{
  L. Sachs (1997), \emph{Angewandte Statistik}. Berlin: Springer.
  
  C.-D. Schoenwiese (1992), \emph{Praktische Statistik}. Berlin:
  Gebr. Borntraeger.

  W. A. Wallis and G. H. Moore (1941): A significance test for time series
  and other ordered observations. Tech. Rep. 1. National Bureau of
  Economic Research. New York.
}
\author{
  T. Pohlert
}

\seealso{
  \code{\link{mk.test}}
}
\examples{
## Example from Schoenwiese (1992, p. 113)
## Number of frost days in April at Munich from 1957 to 1968
## z = -0.124, Accept H0
frost <- ts(data=c(9,12,4,3,0,4,2,1,4,2,9,7), start=1957)
wm.test(frost)

## Example from Sachs (1997, p. 486)
## z = 2.56, Reject H0 on a level of p < 0.05
x <- c(5,6,2,3,5,6,4,3,7,8,9,7,5,3,4,7,3,5,6,7,8,9)
wm.test(x)

data(nottem)
wm.test(nottem)

data(maxau)
wm.test(maxau[, "s"])
}

\keyword{htest}
\keyword{nonparametric}
