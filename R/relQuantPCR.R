## Compute relative expression level for realtime quantitative RT-PCR data
## x: vector containing RT-PCR data
## E: RT-PCR efficiency
## na.rm: remove NA values
relQuantPCR <- function(x, E = 2, na.rm = FALSE){ 
  if(!is.numeric(x) && !is.complex(x) && !is.logical(x)){
    warning("argument is not numeric or logical: returning NA")
    return(as.numeric(NA))
  }
  if(na.rm) x <- x[!is.na(x)]

  return(E^(min(x, na.rm = TRUE)-x))
}
