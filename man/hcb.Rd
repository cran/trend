% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hcb.R
\docType{data}
\name{hcb}
\alias{hcb}
\title{Monthly concentration of particle bound HCB, River Rhine}
\format{
a time series object of class "mts"

\itemize{
\item we first column, series of station Weil (RKM 164.3)
\item ka second column, series of station
Karlsruhe-Iffezheim (RKM 333.9)
\item mz third column, series of station Mainz (RKM 498.5)
\item ko fourth column, series of station Koblenz (RKM 590.3)
\item bh fith column, series of station Bad Honnef(RKM 645.8)
\item bi sixth column, series of station Bimmen (RKM 865.0)
}
}
\source{
International Commission for the Protection of the River Rhine
}
\usage{
data(hcb)
}
\description{
Time series of monthly concentration of particle bound
 Hexachlorobenzene (HCB) in \eqn{\mu}g/kg at six different
 monitoring sites at the River Rhine, 1995.1-2006.12
}
\details{
NO DATA values in the series were filled with estimated values
using linear interpolation (see \code{\link{approx}}.

The Rhine Kilometer (RKM) is in increasing order from source to mouth
of the River Rhine.
}
\examples{
data(hcb)
plot(hcb)
mult.mk.test(hcb)
}
\references{
T. Pohlert, G. Hillebrand, V. Breitung (2011), Trends of persistent
organic pollutants in the suspended matter of the River Rhine,
\emph{Hydrological Processes} 25, 3803--3817.
\doi{10.1002/hyp.8110}
}
\keyword{datasets}
