###############################################################################
## normalization of real-time quantitative RT-PCR data
###############################################################################
## relData: matrix or data.frame containing relative quantities (genes in columns)
## HKs: integer, column numbers of housekeeping genes
## method: method for the computation
## na.rm: remove NA values
normPCR <- function(relData, HKs, method = "Vandesompele", na.rm = FALSE){
  if(!is.matrix(relData) & !is.data.frame(relData))
    stop("'relData' needs to be of class matrix or data.frame")

  if(method == "Vandesompele" & length(HKs) < 2)
    stop("you need data from at least 2 housekeeping genes")
  n <- ncol(relData)
  if(method == "Vandesompele" & n < 3)
    stop("you need data from at least 1 gene of interest")

  if(method == "Vandesompele"){
    NF <- apply(relData[,HKs], 1, geomMean, na.rm = na.rm)
    exprData <- relData[,-HKs]/NF
  }else{
    stop("specified method not yet implemented")
  }
}
