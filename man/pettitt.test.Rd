\name{pettitt.test}
\Rdversion{1.1}
\alias{pettitt.test}
%\alias{pettitt.f}

\title{
 Pettitt's test for change-point-detection

}
\description{
 Performes a nonparametric U-test after Pettitt in order to test for a
 shift in the central tendcency of a time series. The
 H0-hypothesis is: no change is tested against the HA-Hypothesis: change.
}
\usage{
pettitt.test(x)
}

\arguments{
  \item{x}{
    a vector or a time series object of class "ts"

}
}

\details{
  The test is conducted in accordance to Pettitt's proposed
  method. First a vector that comprises the U-Statistics is calculated:

 \deqn{U_{i} = \sum_{k=1}^i{~\sum_{j=i+1}^n{ \textnormal{sgn}(x_k -
       x_j)}}}{%
   U_{i} = \sum_{k=1}^i{~\sum_{j=i+1}^n{ \textnormal{sgn}(x_k -
       x_j)}}
 }

 for \eqn{(1 \le i \le n )}. For a two-sided test, the maximum of the
 absolute obtained U-statistics is taken as the test statistic K:

 \deqn{K = \max \left| U_i \right|}{%
K = max(abs(U(i))).
 }

 Finally, the approximate probability for a two-sided test is calculated
 according to

 \deqn{2 \exp^{-6K^2 / (T^3 + T^2)}}{%
   2 exp[-6K^2 / (T^3 + T^2)].
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
    the maximum K of the absolute values of the U-statistics
  }
  \item{estimate}{
    probable change point at time step tau
  }
  \item{nobs}{
    number of total observations
  }
  \item{alternative}{
    a character string describing the alternative hypothesis
  }
  \item{p.value}{
    the p-value for the test
  }
}

\section{Warning}{
  The current function is for complete observations only.
  The approximate probability is good for \eqn{p \le 0.5}.
}

\references{
CHR (ed., 2010): Das Abflussregime des Rheins und seiner Nebenfluesse im
20. Jahrhundert, Report no I-22 of the CHR, p. 172.

Pettitt, A. N., (1979). A non-parametric approach to the change point
problem. Journal of the Royal Statistical Society Series C, Applied
Statistics 28, 126-135.

%% ~put references to the literature/web site here ~
}
\author{
  %%  ~~who you are~~
  T. Pohlert
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
%%\section{Warning}{
%%  Current Version is for testing purposes!!
%%  Probability calculation needs to be checked with original publication
%%  of Pettitt 1979. Current Version is for complete observations only.
%%}
\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
  \code{\link[strucchange]{efp}}
  \code{\link[strucchange]{sctest.efp}}
}
\examples{
data(maxau) ; plot(maxau[,"s"])
s.res <- pettitt.test(maxau[,"s"])
n <- s.res$nobs
i <- s.res$estimate
s.1 <- mean(maxau[1:i,"s"])
s.2 <- mean(maxau[(i+1):n,"s"])
s <- ts(c(rep(s.1,i), rep(s.2,(n-i))))
tsp(s) <- tsp(maxau[,"s"])
lines(s, lty=2)
print(s.res)


data(PagesData) ; pettitt.test(PagesData)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
\keyword{nonparametric}% __ONLY ONE__ keyword per line
\keyword{univar}
