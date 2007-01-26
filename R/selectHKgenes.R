###############################################################################
## selection of housekeeping (HK) genes for real-time quantitative RT-PCR data
###############################################################################
## relData: matrix or data.frame with data from at least 3 HK genes
## method: method for the computation
## minNrHK: minimum number of HK genes which should be considered
## geneSymbol: character, gene symbols for genes
## trace: locical, print information
## na.rm: remove NA values
selectHKgenes <- function(relData, method = "Vandesompele", minNrHK = 2, 
                          geneSymbol, trace = TRUE, na.rm = FALSE){
  if(!is.matrix(relData) & !is.data.frame(relData))
    stop("'relData' needs to be of class matrix or data.frame")
  n <- ncol(relData)
  if(n < 3)
    stop("you need data from at least 3 genes")
  if(minNrHK >= n)
    stop("'minNrHK' must be smaller than 'ncol(relData)'")
  if(minNrHK < 2){
    warning("'minNrHK' is set to 2")
    minNrHK <- 2
  }
  if(length(geneSymbol) != n)
    stop("'geneSymbol' has wrong length")

  if(method == "Vandesompele"){
    V <- numeric(n-minNrHK)
    names(V) <- paste(((n-1):minNrHK), "/", (n:(minNrHK+1)), sep = "")
    meanM <- numeric(n-minNrHK+1)
    names(meanM) <- as.character(n:minNrHK)
    R <- character(n)
    names(R) <- as.character(c(rep(1, minNrHK),(minNrHK+1):length(R)))
    for(i in n:minNrHK){
      M <- geneStabM(relData, na.rm = na.rm)
      names(M) <- geneSymbol
      ind <- which.max(M)
      meanM[n-i+1] <- mean(M)
      if(i == minNrHK)
        R[1:minNrHK] <- geneSymbol
      else
        R[i] <- geneSymbol[ind]
                      
      if(i > 2){
        NF.old <- apply(relData, 1, geomMean, na.rm = na.rm)
        NF.new <- apply(relData[,-ind], 1, geomMean, na.rm = na.rm)
        V[n-i+1] <- sd(log2(NF.new/NF.old), na.rm = TRUE)
      }
      
      if(trace){
        cat("###############################################################\n")
        cat("Step ", n-i+1, ":\n")
        cat("gene expression stability values M:\n")
        print(sort(M))
        cat("average expression stability M:\t", meanM[n-i+1], "\n")
        if(i > 2){
          cat("gene with lowest stability (largest M value):\t", geneSymbol[ind], "\n")
          cat("Pairwise variation, (", i-1, "/", i, "):\t", V[n-i+1], "\n")
        }
      }
      relData <- relData[,-ind]
      geneSymbol <- geneSymbol[-ind]
    }
    return(list(ranking = R, variation = V, meanM = meanM))
  }else{
    stop("specified method not yet implemented")
  }
}
