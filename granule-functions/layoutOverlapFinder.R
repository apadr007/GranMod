layoutOverlapFinder = function(m, x){
  #m is layout.old matrix
  #x is a 2D position of the node you are querying about moving
  #val = NA
  #if ( isTRUE(is.na(x)) ){ 
  #    FALSE }
  #val <- any(apply(m,1,function(n,x) all(n==x),x=x))
  #if (is.na(val)==TRUE){ val <- FALSE}
  #return(val)
  
  paste(x[1], x[2], sep='&') %in% paste(m[,1], m[,2], sep='&')
}