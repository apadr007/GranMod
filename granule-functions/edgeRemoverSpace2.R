edgeRemoverSpace2 <- function(g, layout.old, near){
  n = sapply(v, neighbors, g)
  
  for (v in 1:nrow(layout.old) ) {
    print(v)
    f = near[[v]]
    b = which( !n[[v]]%in%f )
    if (isEmpty(b)==FALSE){
      for ( i in 1:length(b) ) {
        b.idx = b[i]
        n.rdy = n[ b.idx ]
        g = delete.edges(g, E(g, P=c(v, n.rdy)) )
      }
    } else {  }
  }
  return(g)
}