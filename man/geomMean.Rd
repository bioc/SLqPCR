\name{geomMean}
\alias{geomMean}
\title{ Geometric Mean }
\description{
  Computation of the geometric mean.
}
\usage{
geomMean(x, na.rm = FALSE)
}
\arguments{
  \item{x}{ numeric vector of non-negative Reals }
  \item{na.rm}{ a logical value indicating whether \code{NA} values should be
          stripped before the computation proceeds.  }
}
\details{
  The computation of the geometric mean is done via \code{prod(x)^(1/length(x))}.
}
\value{
  geometric mean
}
%\references{ }
\author{ Dr. Matthias Kohl (SIRS-Lab GmbH) \email{kohl@sirs-lab.com}}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{}
%\examples{
%}
\keyword{data}
