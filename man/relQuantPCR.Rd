\name{relQuantPCR}
\alias{relQuantPCR}
\title{ Compute relative expression values for realtime quantitative RT-PCR data }
\description{
  Compute relative expression values for realtime quantitative RT-PCR data based on 
  Ct or take-off values, respectively. The computations use the PCR efficiency.
 }
\usage{
relQuantPCR(x, E = 2, na.rm = FALSE)
}
\arguments{
  \item{x}{ numeric vector containing raw data }
  \item{E}{ PCR efficiency }
  \item{na.rm}{ a logical value indicating whether \code{NA} values should be
          stripped before the computation proceeds.  }
}
%\details{
%}
\value{
  vector of relative expression values w.r.t. specified PCR efficiency.
}
\references{ 
  Jo Vandesompele, Katleen De Preter, Filip Pattyn et al. (2002). Accurate 
  normalization of real-time quantitative RT-PCR data by geometric averaging 
  of multiple internal control genes. 
  Genome Biology 2002. 3(7):research0034.1-0034.11.
  \url{http://genomebiology.com/2002/3/7/research/0034/}
}
\author{ Dr. Matthias Kohl (SIRS-Lab GmbH) \email{kohl@sirs-lab.com}}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{}
%\examples{
%}
\keyword{data}
