################################################################################
## Vandesompele et al. (2002): gene-stability measure M
################################################################################
## relData: matrix or data.frame containing real-time quant. RT-PCR data
## na.rm: remove NA values
geneStabM <- function(relData, na.rm = FALSE){
  if(!is.data.frame(relData) & !is.matrix(relData))
    stop("'relData' has to of class matrix or data.frame")

  n <- ncol(relData)
  if(n == 1) stop("you need at least two genes for this computation")

  M <- numeric(n)
  for(j in 1:n){
    A <- log2(relData[,j]/relData[,-j])
    if(n > 2)
      M[j] <- mean(apply(A, 2, sd, na.rm = na.rm))
    else
      M[j] <- sd(A, na.rm = na.rm)
  }
  if(is.data.frame(relData))
    names(M) <- names(relData)
  else
    names(M) <- colnames(relData)

  return(M)
}
