missby <- function(ruling, bounds, truth){
  if(ruling == TRUE  | is.na(ruling)) return(NA)
  else if( ruling == FALSE & truth$value > bounds[2]){
   miss <- truth$value - bounds[2]
  }else if( ruling == FALSE & truth$value < bounds[2]){
    miss <- bounds[1] - truth$value
  }

  return(miss)
}
