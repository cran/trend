\name{trend-package}
\Rdversion{0.0.2}
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
Version: \tab 0.2.0\cr
Date: \tab 2016-05-14\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
}

\tabular{ll}{
csmk.test           \tab    Correlated Seasonal Mann-Kendall Test\cr
cs.test            \tab      Cox and Stuart trend test\cr
maxau  \tab                 Suspended sediment concentration (s) and\cr
       \tab              discharge (Q) for the River Rhine at Maxau,\cr
       \tab      annual averages\cr
mk.test             \tab    Mann-Kendall Trend Test\cr
PagesData           \tab    Simulated data of Page \cr
partial.cor.trend.test    \tab  Partial correlation trend test\cr
partial.mk.test     \tab    Partial Mann-Kendall Test\cr
pettitt.test        \tab     Pettitt's test for change-point-detection\cr
sea.sens.slope   \tab     Seasonal Sen's slope and intercept \cr
sens.slope          \tab    Sen's slope and intercept\cr
smk.test            \tab    Seasonal Mann-Kendall Test\cr
summary.trend.test  \tab    Summary Method for Class 'trend'\cr
wm.test            \tab      Wallis and Moore phase-frequency test\cr
ww.test            \tab      Wald-Wolfowitz test for independence and stationarity\cr
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

  D. R. Cox and A. Stuart (1955), Quick sign tests for trend in location
  and dispersion. Biometrika 42, 80-95.
    
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

  Richard O. Gilbert. (1987). Statistical methods for environmental
  pollution monitoring. John Wiley & Sons. pp 230

  R. K. Rai, A. Upadhyay, C. S. P. Ojha and L. M. Lye (2013),
  Statistical analysis of hydro-climatic variables. In:
  R. Y. Surampalli, T. C. Zhang, C. S. P. Ojha, B. R. Gurjar,
  R. D. Tyagi and C. M. Kao (ed. 2013), \emph{Climate change modelling,
    mitigation, and adaptation}. Reston, VA: ASCE. doi =
  10.1061/9780784412718.
  
  L. Sachs (1997), \emph{Angewandte Statistik}. Berlin: Springer.
  
  Sen, P.K., (1968): Estimates of the regression coefficient based on
  Kendall's tau, Journal of the American Statistical Association, 63,
  1379--1389.

  C.-D. Schoenwiese (1992), \emph{Praktische Statistik}. Berlin:
  Gebr. Borntraeger.

  W. A. Wallis and G. H. Moore (1941): A significance test for time series
  and other ordered observations. Tech. Rep. 1. National Bureau of
  Economic Research. New York.

  A. Wald and J. Wolfowitz (1943), An exact test for randomness in the
  non-parametric case based on serial correlation.
  \emph{Annual Mathematical Statistics} 14, 378--388.

  WMO (2009), \emph{Guide to Hydrological Practices}. Volume II,
  Management of Water Resources and Application of Hydrological
  Practices, WMO-No. 168.
}

\keyword{ package }

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
\code{\link{summary.trend.test}},
\code{\link{wm.test}},
\code{\link{ww.test}},
\code{\link{cs.test}}
}
