\name{SLqPCRdata}
\alias{SLqPCRdata}
\docType{data}
\title{ SIRS-Lab inhouse qPCR data }
\description{
  This data is part of a SIRS-Lab inhouse real-time quantitative PCR experiment.
}
\usage{data(SLqPCRdata)}
\format{
  A data frame with 16 observations on the following 4 variables.
  \describe{
    \item{\code{Gene1}}{a numeric vector, average take-off values of gene 1 }
    \item{\code{Gene2}}{a numeric vector, average take-off values of gene 2 }
    \item{\code{HK1}}{a numeric vector, average take-off values of housekeeper 1 }
    \item{\code{HK2}}{a numeric vector, average take-off values of housekeeper 2 }
  }
}
\details{
  The row names of this data set indicate the probes which were investigated. 
  The take-off values are mean values of three replicates. 
}
\source{
  \url{www.sirs-lab.com}
}
\references{
  \url{www.sirs-lab.com}
}
\examples{
data(SLqPCRdata)
SLqPCRdata
}
\keyword{datasets}
