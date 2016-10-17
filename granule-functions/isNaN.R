isNaN = function(x){
  x[is.na(x) == TRUE] = 0
  return(x)
}
