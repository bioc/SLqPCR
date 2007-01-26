\name{geneStabM}
\alias{geneStabM}
\title{ Gene expression stability value M }
\description{
  Computation of the gene expression stability value M for real-time 
  quantitativ RT-PCR data. For more details we refer to 
  Vandesompele et al. (2002). 
}
\usage{
geneStabM(relData, na.rm = FALSE)
}
\arguments{
  \item{relData}{ matrix or data.frame containing real-time quantitative 
  		  RT-PCR data }
  \item{na.rm}{ a logical value indicating whether \code{NA} values should be
          stripped before the computation proceeds. }
}
\details{
  The gene expression stability value M is defined as the average pairwise
  normalization factor; i.e., one needs to specify data from at least two
  genes. For more details see Vandesompele et al. (2002).
}
\value{
  numeric vector with gene expression stability values 
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
