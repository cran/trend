\name{sens.slope}
\Rdversion{1.1}
\alias{sens.slope}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Sen's slope and intercept
%%  ~~function to do ... ~~
}
\description{
  Computes Sen's slope  for linear rate of change and the corresponding intercept
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sens.slope(x, level = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    an time series object of class "ts"
%%     ~~Describe \code{x} here~~
}
\item{level}{
  the level of significance
%%     ~~Describe \code{level} here~~
}
}

\details{
  This test computes both the slope (i.e. linear rate of change) and
  intercept according to Sen's method. First, a set of linear slopes is
  calculated as follows:
  \deqn{d_{k} = \frac{x_j - x_i}{j - i}}{%
    d(k) = (x(j) - x(i)) / (j - i)}

  for \eqn{\left(1 \le i < j \le n \right)}{(1 <= i < j <= n)}, where d
  is the slope, x denotes the variable, n is the number of data, and i,
  j are indices.

  Sen's slope is then calculated as the median from all slopes:
  \eqn{b_{Sen} = \textnormal{median}(d_k)}{b = Median(d(k))}. The
  intercepts are computed for each timestep i as given by

  \deqn{a_i = x_i - b_{Sen} \times i}{%
    a(i) = x(i) - b * i}

  and the corresponding intercept is as well the median of all
  intercepts.

  This function also computes the upper and lower confidence limits for
  sens slope.

}


\value{
  A list of class "trend.test".
  
  \item{b.sen}{
    The value of Sen's slope
  }

  \item{b.sen.up, b.sen.lo}{
    The corresponding value of Sen's slope for
    the upper (lower) confidence limit
  }
  \item{intercept}{
    The intercept
  }
  \item{nobs}{
    Number of records
  }
  \item{D}{
    A vector of individual slopes
  }
  \item{varS}{
    The variance of the Mann-Kendall score
  }
  \item{conf.int}{
    The confidence intervall for the upper and lower bounds
    }
  
}

\references{
    Hipel, K.W. and McLeod, A.I., (2005). Time Series Modelling of Water
    Resources and Environmental Systems. \url{http://www.stats.uwo.ca/faculty/aim/1994Book/}.

    Sen, P.K., (1968): Estimates of the regression coefficient based on
    Kendall's tau, Journal of the American Statistical Association, 63,
    1379--1389.
}
\author{
  %%  ~~who you are~~
  T. Pohlert
}
%% \note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning}{
  Current Version is for complete observations only.
}
\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
  \code{\link{mk.test}},
  \code{\link[zyp]{zyp.sen}},
  \code{\link{lm}}
}
\examples{
data(maxau) ; plot(maxau[,"s"])
sen.res <- sens.slope(maxau[,"s"])
t <- (1:(length(maxau[,"s"])))
s.pred <- sen.res$intercept + sen.res$b.sen * t
s.pred.ts <- ts(s.pred)
tsp(s.pred.ts) <- tsp(maxau[,"s"])
lines(s.pred.ts, lty=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
\keyword{nonparametric}% __ONLY ONE__ keyword per line
\keyword{univar}
