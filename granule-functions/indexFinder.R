indexFinder <- function(c, z){
  l <- unlist(apply(c, 2, list), recursive=F)
  logic <- mapply(function(x,y)x==y, l, z)
  output <- which(.rowSums(logic, m=nrow(logic), n=ncol(logic)) == ncol(logic))
  return(output)
  
}