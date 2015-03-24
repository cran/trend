\name{csmk.test}
\Rdversion{1.1}
\alias{csmk.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  Correlated Seasonal Mann-Kendall Test
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Performs a Seasonal Mann-Kendall test under the presence of correlated
  seasons.
}
\usage{
csmk.test(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    A time series object that comprises \eqn{>= 2} seasons
}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
  The Mann-Kendall statistics are first computed for each season seperately.
  The variance - covariance matrix is computed
  according to Libiseller and Grimvall (2002) and ties are taken into
  account. Finally the corrected Z-statistics for the entire series
  is calculated as follows:

  \deqn{
    Z = \mathbf{1}^T \mathbf{S} / \mathbf{1}^T \mathbf{X}~\mathbf{1}}{%
    Z = 1^T S / 1^T X 1}

  Z denotes the quantile of the normal distribution, 1 indicates a vector
  with all elements equal to one, S is the vector of Mann-Kendall scores
  for each season and X denotes the variance - covariance matrix.
}
\value{
\item{method}{
  The chosen method (i.e. CSMK)
}
\item{Sg}{
  Vector of Mann-Kendall scores for each season
}
\item{varSg}{
  Vector of the variance of Mann-Kendall scores for each season
}
\item{Zg}{
  Vector of Z-values for each season
}
\item{pvalg}{
  Vector of p-values for each season
}
\item{Covar}{
  Variance - Covariance matrix for the seasonal Mann-Kendall statistics
}
\item{Correl}{
  Correlation matrix for the seasonal Mann-Kendall statistics
}
\item{Stot}{
  Mann-Kendall score for the entire series
}
\item{Z}{
  Z-value for the entire series
}
\item{Varianz}{
  variance for the entire series
}
\item{pvalue}{
  p-value for the entire series
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
  \code{\link{cor}},
  \code{\link{cor.test}},
  \code{\link{mk.test}},
  \code{\link{smk.test}},
  \code{\link[Kendall]{SeasonalMannKendall}},
  \code{\link{summary.trend.test}}
}
\examples{
data(co2)
plot(co2)
res <- csmk.test(co2)
summary.trend.test(res)

data(nottem)
plot(nottem)
res <- csmk.test(nottem)
summary.trend.test(res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
\keyword{nonparametric}% __ONLY ONE__ keyword per line
\keyword{multivariate}
