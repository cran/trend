\name{partial.mk.test}
\Rdversion{1.1}
\alias{partial.mk.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Partial Mann-Kendall Test
%%  ~~function to do ... ~~
}
\description{
  Performs a partial Mann-Kendall test
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
partial.mk.test(x, z)
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
}

\details{
  According to Libiseller and Grimvall (2002), the expected value for
  the Mann-Kendall Score of variable "a" that is correlated with the covariate
  "b" is:
  \eqn{E = S_b * \sigma_{ab} / \sigma_{bb}}{%
    E = S(b) * sigma(ab) / sigma(bb)
  }

  Furthermore, the expected variance of "a" is:
  \eqn{Var = \sigma_{aa} - \sigma_{ab} / \sigma_{bb} * \sigma_{ba}}{%
    Var = sigma(aa) - sigma(ab) / sigma(bb) * sigma(ba)
  }

  Finally, the Z-Quantile is defined as:
  \eqn{Z = (S_a - E) / \sqrt(Var)}{%
    Z = (S(a) - E) / sqrt(Var)
  }
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
    	the value of the test statistic
  }
  \item{estimate}{
    the Mann-Kendall score S and the variance varS
  }

  \item{alternative}{
    a character string describing the alternative hypothesis
  }
  \item{p.value}{
    the p-value of the test
  }
}

\references{
Libiseller, C. and Grimvall, A., (2002).
Performance of partial Mann-Kendall tests for trend detection in the
presence of covariates. Environmetrics 13, 71-84, \url{http://dx.doi.org/10.1002/env.507}.
}
\author{
  T. Pohlert
}

\section{Warning}{
  Current Version is for complete observations only.
  The "method=='Partial'" is in testing mode.
}
\seealso{
  \code{\link{partial.cor.trend.test}},
  \code{\link{summary.trend.test}}
}
\examples{
data(maxau)
s <- maxau[,"s"]; Q <- maxau[,"Q"]
partial.mk.test(s,Q)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
\keyword{nonparametric}% __ONLY ONE__ keyword per line
\keyword{univar}
