\name{maxau}
\Rdversion{1.1}
\alias{maxau}
\docType{data}
\title{
  %%   ~~ data name/kind ... ~~
  Annual suspended sediment concentration and flow data, River Rhine}
\description{
  %%  ~~ A concise (1-5 lines) description of the dataset. ~~
  Annual time series of average suspended sediment
  concentration (s) in mg/l and average discharge (Q) in m^3 / s at
  the River Rhine, 1965.1-2009.1
}
\usage{data(maxau)}
\format{
  The format is:
 mts [1:45, 1:2] 37.4 31.3 26.4 28.8 32.9 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:2] "s" "Q"
 - attr(*, "tsp")= num [1:3] 1965 2009 1
 - attr(*, "class")= chr [1:2] "mts" "ts"
}

\source{

BfG, Germany 
}

\examples{
data(maxau)
plot(maxau)
}
\keyword{datasets}
