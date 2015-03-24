\name{sea.sens.slope}
\Rdversion{1.1}
\alias{sea.sens.slope}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Seasonal Sen's slope and intercept
%%  ~~function to do ... ~~
}
\description{
  Computes seasonal Sen's slope  for linear rate of change and the corresponding intercept
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sea.sens.slope(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    an time series object of class "ts"
  }
}

\value{
  A list of class "trend.test".

  \item{b.sen}{
    The value of Sen's slope
  }


  \item{intercept}{
    The intercept
  }
  \item{nobs}{
    Number of records
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
sea.sens.slope(nottem)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
\keyword{nonparametric}% __ONLY ONE__ keyword per line
\keyword{univar}
