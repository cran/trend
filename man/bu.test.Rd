% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bu.test.R
\name{bu.test}
\alias{bu.test}
\title{Buishand U Test for Change-Point Detection}
\usage{
bu.test(x, m = 20000)
}
\arguments{
\item{x}{a vector of class "numeric" or a time series object of class "ts"}

\item{m}{numeric, number of Monte-Carlo replicates, defaults to 20000}
}
\value{
A list with class "htest" and "cptest"
\item{data.name}{character string that denotes the input data}
\item{p.value}{the p-value}
\item{statistic}{the test statistic}
\item{null.value}{the null hypothesis}
\item{estimates}{the time of the probable change point}
\item{alternative}{the alternative hypothesis}
\item{method}{character string that denotes the test}
\item{data}{numeric vector of Sk for plotting}
}
\description{
Performes the Buishand U test for change-point detection
of a normal variate.
}
\details{
Let \eqn{X} denote a normal random variate, then the following model
with a single shift (change-point) can be proposed:

\deqn{
  x_i = \left\{
      \begin{array}{lcl}
       \mu + \epsilon_i, & \qquad & i = 1, \ldots, m \\
       \mu + \Delta + \epsilon_i & \qquad & i = m + 1, \ldots, n \\
      \end{array} \right.}{%
  x[i] = \mu + \epsilon[i] for i = 1, ..., m and x[i] = \mu + \delta
  + \epsilon_i for i = m + 1, ..., n}

with \eqn{\epsilon \approx N(0,\sigma)}. The null hypothesis \eqn{\Delta = 0}
is tested against the alternative \eqn{\Delta \ne 0}{\delta != 0}.

In the Buishand U test, the rescaled adjusted partial sums
are calculated as

\deqn{S_k = \sum_{i=1}^k \left(x_i - \bar{x}\right) \qquad (1 \le i \le n)}{%
S[k] = \sum (x[i] - xmean)   (1, <= i <= n)}

The sample standard deviation is
\deqn{
   D_x = \sqrt{n^{-1} \sum_{i=1}^n \left(x_i - \bar{x}\right)}}{%
   Dx = sqrt(n^(-1) \sum(x - \mu))}

The test statistic is calculated as:
\deqn{U = \left[n \left(n + 1 \right) \right]^{-1} \sum_{k=1}^{n-1}
\left(S_k / D_x \right)^2 }{%
U = 1 / [n * (n + 1)] * \sum_{k=1}^{n-1} (S[k] - Dx)^2}. 

The \code{p.value} is estimated with a Monte Carlo simulation
using \code{m} replicates.

Critical values based on \eqn{m = 19999} Monte Carlo simulations
are tabulated for \eqn{U} by Buishand (1982, 1984).
}
\note{
The current function is for complete observations only.
}
\examples{
data(Nile)
(out <- bu.test(Nile))
plot(out)

data(PagesData)
bu.test(PagesData)
 
}
\references{
T. A. Buishand (1982), Some Methods for Testing the Homogeneity
of Rainfall Records, \emph{Journal of Hydrology} 58, 11--27.

T. A. Buishand (1984), Tests for Detecting a Shift in the Mean of
Hydrological Time Series, \emph{Journal of Hydrology} 73, 51--69.
}
\seealso{
\code{\link[strucchange]{efp}}
\code{\link[strucchange]{sctest.efp}}
}
\keyword{htest}
\keyword{ts}
\keyword{univar}
