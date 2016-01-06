\name{cs.test}
\Rdversion{1.1}
\alias{cs.test}

\title{
 Cox and Stuart trend test

}
\description{
  Performes the non-parametric Cox and Stuart trend test
 (two-sided test).
}
\usage{
cs.test(x)
}

\arguments{
  \item{x}{
    a vector or a time series object of class "ts"
  }
}

\details{
  First, the series is devided by three. It is compared, whether the
  data of the first third of the series are larger or smaller than the
  data of the last third of the series.
  The test statistic of the Cox-Stuart trend test for \eqn{n > 30}
  is calculated as:
  
 \deqn{z = \frac{\| S - \frac{n}{6} \|}{\sqrt{\frac{n}{12}}}}{%
   z = abs(S - n / 6) / sqrt(n / 12)}
 
 where \eqn{S} denotes the maximum of the number of signs, i.e. \eqn{+}
 or \eqn{-}, respectively. The \eqn{z}-statistic is
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
    the Cox-Stuart z-value
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

  D. R. Cox and A. Stuart (1955), Quick sign tests for trend in location
  and dispersion. \emph{Biometrika} 42, 80-95.
}
\author{
  T. Pohlert
}

\seealso{
  \code{\link{mk.test}}
}
\examples{
## Example from Schoenwiese (1992, p. 114)
## Number of frost days in April at Munich from 1957 to 1968
## z = -0.5, Accept H0
frost <- ts(data=c(9,12,4,3,0,4,2,1,4,2,9,7), start=1957)
cs.test(frost)

## Example from Sachs (1997, p. 486-487)
## z ~ 2.1, Reject H0 on a level of p = 0.0357
x <- c(5,6,2,3,5,6,4,3,7,8,9,7,5,3,4,7,3,5,6,7,8,9)
cs.test(x)

data(Nile)
cs.test(Nile)
}

\keyword{htest}
\keyword{nonparametric}
