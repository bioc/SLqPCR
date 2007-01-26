\name{normPCR}
\alias{normPCR}
\title{ Normalization of real-time quantitative RT-PCR data }
\description{
  This function can be used to normalize real-time quantitative RT-PCR data.
}
\usage{
normPCR(relData, HKs, method = "Vandesompele", na.rm = FALSE)
}
\arguments{
  \item{relData}{ matrix or data.frame containing relative quantities (genes in columns) }
  \item{HKs}{ integer, column numbers of housekeeping genes }
  \item{method}{ method for the computation }
  \item{na.rm}{ a logical value indicating whether \code{NA} values should be
          stripped before the computation proceeds. }
}
\details{
  This function can be used to normalize real-time quantitative RT-PCR data. 
  The default method \code{"Vandesompele"} was proposed by 
  Vandesompele et al. (2002).

  Currently, only the method by Vandesompele et al. (2002) is implemented.
}
\value{Normalized expression data}
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
\examples{
data(SLqPCRdata)
relData <- apply(SLqPCRdata, 2, relQuantPCR)
geneStabM(relData[,c(3,4)])
exprData <- normPCR(SLqPCRdata, c(3,4))
}
\keyword{data}
