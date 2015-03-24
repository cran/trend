\name{mk.test}
\Rdversion{1.1}
\alias{mk.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  Mann-Kendall Test
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Performs a univariate Mann-Kendall test
}
\usage{
mk.test(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    A time series object
%%     ~~Describe \code{x} here~~
}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
  The Mann-Kendall statistics is computed
  taking ties into account.
}
\value{
  %%  ~Describe the value returned
%%  If it is a LIST, use
\item{method}{
  The chosen method (i.e. MK)
}
\item{Sg, Stot}{
  Mann-Kendall score
}
\item{varSg, Varianz}{
  variance of Mann-Kendall score
}
\item{Zg, Z}{
  Z-value (quantile of the normal distribution)
}
\item{pvalg, pvalue}{
  p-value
}
\item{Covar}{
  Variance - Covariance matrix
}
\item{Correl}{
  Correlation matrix
}

 Generic function summary.trend.test is provided.
}
\references{
  %% ~put references to the literature/web site here ~
  Hipel, K.W. and McLeod, A.I., (2005).
Time Series Modelling of Water Resources and Environmental Systems.
Electronic reprint of our book orginally published in 1994.
\url{http://www.stats.uwo.ca/faculty/aim/1994Book/}.

Libiseller, C. and Grimvall, A., (2002).
Performance of partial Mann-Kendall tests for trend detection in the
presence of covariates. Environmetrics 13, 71-84, \url{http://dx.doi.org/10.1002/env.507}.
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

  \code{\link{cor.test}},
  \code{\link[Kendall]{MannKendall}},
  \code{\link{partial.mk.test}},
  \code{\link{sens.slope}}
  \code{\link{summary.trend.test}}
}
\examples{
data(Nile)
plot(Nile)
res <- mk.test(Nile)
summary.trend.test(res)

## Compare with cor()
n <- length(Nile)
cor.test(x=(1:n),y=Nile, meth="kendall")

data(maxau)
Q <- maxau[,"Q"]
plot(Q)
res <- mk.test(Q)
summary.trend.test(res)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
\keyword{nonparametric}% __ONLY ONE__ keyword per line
\keyword{univar}
