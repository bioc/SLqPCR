\name{selectHKgenes}
\alias{selectHKgenes}
\title{ Selection of reference/housekeeping genes }
\description{
  This function can be used to determine a set of reference/housekeeping (HK) 
  genes for gene expression experiments. 
}
\usage{
selectHKgenes(relData, method = "Vandesompele", minNrHK = 2, geneSymbol, 
	      trace = TRUE, na.rm = FALSE)
}
\arguments{
  \item{relData}{ matrix or data.frame containing relative expression values }
  \item{method}{ method to compute most stable genes }
  \item{minNrHK}{ minimum number of HK genes that should be considered }
  \item{geneSymbol}{ gene symbols }
  \item{trace}{ logical, print additional information }
  \item{na.rm}{ a logical value indicating whether \code{NA} values should be
          stripped before the computation proceeds.  }
}
\details{
  This function can be used to determine a set of reference/housekeeping (HK) genes
  for gene expression experiments. The default method \code{"Vandesompele"}
  was proposed by Vandesompele et al. (2002).

  Currently, only the method by Vandesompele et al. (2002) is implemented.

  Vandesompele et al. (2002) propose a cut-off value of 0.15 for the 
  pairwise variation. Below this value the inclusion of an additional
  housekeeping gene is not required. 
}
\value{
  If \code{method = "Vandesompele"} a list with the following components is
  returnd
  \item{ranking }{ ranking of genes from best to worst where the two most
  		   stable genes cannot be ranked }
  \item{variation }{ pairwise variation during stepwise selection }
  \item{meanM }{ average expression stability M }
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
\examples{
data(vandesompele)
res.BM <- selectHKgenes(vandesompele[1:9,], method = "Vandesompele", geneSymbol = names(vandesompele), minNrHK = 2, trace = TRUE, na.rm = TRUE)
}
\keyword{data}
