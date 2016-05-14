\name{ww.test}
\Rdversion{1.1}
\alias{ww.test}

\title{
 Wald-Wolfowitz test for independence and stationarity

}
\description{
  Performes the non-parametric Wald-Wolfowitz test for independence and stationarity
 (two-sided test).
}
\usage{
ww.test(x)
}

\arguments{
  \item{x}{
    a vector or a time series object of class "ts"
  }
}

\details{
  Let \eqn{x_1, x_2, ..., x_n} denote the sampled data, then the test
  statistic of the Wald-Wolfowitz test is calculated as:

  \deqn{R = \sum_{i=1}^{n-1} x_i x_{i+1} + x_1 x_n}{%
    R = sum(x[1:(n-1)] * x[2:n]) + x[1] * x[n]}

  The expected value of R is:

  \deqn{E(R) = \frac{s_1^2 - s_2}{n - 1}}{%
    E(R) = (s1^2 - s2) / (n - 1)}

  The expected variance is:

  \deqn{V(R) = \frac{s_2^2 - s_4}{n - 1} - E(R)^2 + \frac{s_1^4 -
      4 s_1^2 s_2 + 4 s_1 s_3 + s_2^2 - 2 s_4}{(n - 1) (n - 2)}}{%
    V(R) = (s2^2 - s4) / (n - 1) - er^2 +
    (s1^4 - 4 * s1^2 * s2 + 4 * s1 * s3 + s2^2 - 2 * s4) /
    ((n - 1) * (n - 2))}

  with:
	
  \deqn{s_t = \sum_{i=1}^{n} x_i^t, ~~ t = 1, 2, 3, 4}{%
    st = sum(x^t), ~~ t = 1, 2, 3, 4}
  
  For \eqn{n > 10} the test statistic is normally distributed, with:

  \deqn{z = \frac{R - E(R)}{\sqrt{V(R)}}}{%
    z = (R - E(R)) / V(R)^0.5}

  ww.test calculates p-values from the standard normal distribution for
  the two-sided case.
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
    the Wald-Wolfowitz z-value
  }
  \item{alternative}{
    a character string describing the alternative hypothesis
  }
  \item{p.value}{
    the p-value for the test
  }
}

\section{Note}{
  NA values are omitted.
}

\references{
  R. K. Rai, A. Upadhyay, C. S. P. Ojha and L. M. Lye (2013),
  Statistical analysis of hydro-climatic variables. In:
  R. Y. Surampalli, T. C. Zhang, C. S. P. Ojha, B. R. Gurjar,
  R. D. Tyagi and C. M. Kao (ed. 2013), \emph{Climate change modelling,
    mitigation, and adaptation}. Reston, VA: ASCE. doi =
  10.1061/9780784412718.
  
  A. Wald and J. Wolfowitz (1943), An exact test for randomness in the
  non-parametric case based on serial correlation.
  \emph{Annual Mathematical Statistics} 14, 378--388.

  WMO (2009), \emph{Guide to Hydrological Practices}. Volume II,
  Management of Water Resources and Application of Hydrological
  Practices, WMO-No. 168.
}
\author{
  T. Pohlert
}

\examples{
data(nottem)
ww.test(nottem)

data(Nile)
ww.test(Nile)

set.seed(200)
x <- rnorm(100)
ww.test(x)
}

\keyword{htest}
\keyword{nonparametric}
