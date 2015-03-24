\name{trend-package}
\Rdversion{1.1}
\alias{trend-package}
\alias{trend}
\docType{package}
\title{
 Non-Parametric Trend Tests and Change-Point Detection
}
\description{The analysis of environmental data often requires
	     the detection of trends and change-points. 
	     This package provides the Mann-Kendall Trend Test,
             seasonal Mann-Kendall Test,
             correlated seasonal Mann-Kendall Test,
             partial Mann-Kendall Trend test,
	     (Seasonal) Sen's slope, partial correlation trend test and
	     change-point test after Pettitt.
}
\details{
\tabular{ll}{
Package: \tab trend\cr
Type: \tab Package\cr
Version: \tab 0.0.1\cr
Date: \tab 2015-03-23\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
}

\tabular{ll}{
maxau  \tab                 Suspended sediment concentration (s) and\cr
       \tab              discharge (Q) for the River Rhine at Maxau,\cr
                  \tab      annual averages\cr
csmk.test           \tab    Correlated Seasonal Mann-Kendall Test\cr
mk.test             \tab    Mann-Kendall Trend Test\cr
PagesData           \tab    Simulated data of Page \cr
partial.cor.trend.test    \tab  Partial correlation trend test\cr
partial.mk.test     \tab    Partial Mann-Kendall Test\cr
pettitt.test        \tab     Pettitt's test for change-point-detection\cr
sea.sens.slope   \tab     Seasonal Sen's slope and intercept \cr
sens.slope          \tab    Sen's slope and intercept\cr
smk.test            \tab    Seasonal Mann-Kendall Test\cr
summary.trend.test  \tab    Summary Method for Class 'trend'\cr
}
}
\author{
Thorsten Pohlert

Maintainer: Thorsten Pohlert <pohlert@gmx.de>
}
\references{
  Bahrenberg, G., Giese, E. and Nipper, J., (1992): Statistische Methoden
  in der Geographie, Band 2 Multivariate Statistik, Teubner, Stuttgart.

  CHR (ed., 2010): Das Abflussregime des Rheins und seiner Nebenfluesse
  im 20. Jahrhundert, Report no I-22 of the CHR, p. 172.

  Hipel, K.W. and McLeod, A.I., (2005).
  Time Series Modelling of Water Resources and Environmental Systems.
  \url{http://www.stats.uwo.ca/faculty/aim/1994Book/}.

  Hirsch, R., J. Slack and R. Smith, (1982): Techniques of Trend Analysis
  for Monthly Water Quality Data. Water Resour. Res., 18, 107-121.

  Libiseller, C. and Grimvall, A., (2002).
  Performance of partial Mann-Kendall tests for trend detection in the
  presence of covariates. Environmetrics 13,
  71-84, \url{http://dx.doi.org/10.1002/env.507}.

  Pettitt, A. N., (1979). A non-parametric approach to the change point
  problem. Journal of the Royal Statistical Society Series C, Applied
  Statistics 28, 126-135.

  Sen, P.K., (1968): Estimates of the regression coefficient based on
  Kendall's tau, Journal of the American Statistical Association, 63,
  1379--1389. 
}

\keyword{ package }
\keyword{ts}
\keyword{nonparametric}% __ONLY ONE__ keyword per line
\keyword{univar}
\keyword{multivariate}

\seealso{
\code{\link{cor}},
\code{\link{cor.test}},
\code{\link{csmk.test}},
\code{\link{mk.test}},
\code{\link{PagesData}},
\code{\link{partial.mk.test}},
\code{\link{partial.cor.trend.test}},
\code{\link{pettitt.test}},
\code{\link{print.trend.test}},
\code{\link{sea.sens.slope}},
\code{\link{sens.slope}}
\code{\link{smk.test}},
\code{\link{summary.trend.test}}
}
