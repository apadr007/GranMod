isZeroFunc = function(x){
  x[which(truefalse == 0)] <- FALSE
  x[which(truefalse == 1)] <- TRUE
  return(x)
}