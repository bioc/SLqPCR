## Compute geometric mean for numeric vector
## x: numeric vector of non-negative Reals
## na.rm: remove NA values
geomMean <- function(x, na.rm = FALSE){ 
  if(!is.numeric(x) && !is.complex(x) && !is.logical(x)){
    warning("argument is not numeric or logical: returning NA")
    return(as.numeric(NA))
  }
  if(any(x < 0)) 
    stop("'x' contains negative value(s)")
  if(na.rm) x <- x[!is.na(x)]

  return(prod(x)^(1/length(x)))
}
