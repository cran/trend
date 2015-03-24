\name{print.trend.test}
\Rdversion{1.1}
\alias{print.trend.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Print Method for Class 'trend.test'
%%  ~~function to do ... ~~
}
\description{
  Print Method for Class 'trend.test'
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
%\usage{
%print.trend.test(object, ...)
%}
\usage{
 \method{print}{trend.test}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An object of class "trend.test"
%%     ~~Describe \code{object} here~~
}
\item{\dots}{
  Further arguments. (Currently not in use!).
%%     ~~Describe \code{\dots} here~~
}
}

\references{
  %% ~put references to the literature/web site here ~
  Hipel, K.W. and McLeod, A.I., (2005). Time Series Modelling of Water
  Resources and Environmental Systems. \url{http://www.stats.uwo.ca/faculty/aim/1994Book/}.

  Libiseller, C. and Grimvall, A., (2002). Performance of partial
  Mann-Kendall tests for trend detection in the presence of covariates.
  Environmetrics 13, 71-84.
}
\author{
  %%  ~~who you are~~
  T. Pohlert
}

\section{Warning}{
  Current Version is for complete observations only.
}
\seealso{
  %% ~~objects to See Also as \code{\link{help}},
  \code{\link{csmk.test}},
  \code{\link{mk.test}},
  \code{\link{partial.mk.test}},
  \code{\link{smk.test}},
  \code{\link[base]{summary}}
}
\examples{

### ** Examples
data(Nile)
plot(Nile)
res <- mk.test(Nile)
print.trend.test(res)

data(maxau)
Q <- maxau[,"Q"]
plot(Q)
res <- mk.test(Q)
print.trend.test(res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
\keyword{nonparametric}% __ONLY ONE__ keyword per line
\keyword{univar}
\keyword{multivariate}
