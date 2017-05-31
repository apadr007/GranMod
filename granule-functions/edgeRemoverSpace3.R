edgeRemoverSpace3 = function(g, layout.old){
  f = myDistOneFast(layout.old)
  n <- lapply(1:nrow(layout.old), neighbors, graph=g)
  
  for (v in 1:nrow(layout.old)){
    b = which( !n[[v]]%in%f[[v]] )
    if (isEmpty(b)==FALSE){
      for ( i in 1:length(b) ) {
        b.idx = b[i]
        n.rdy = n[ b.idx ]
        g = delete.edges(g, E(g, P=c(v, n.rdy)) )
      }
    }
  }
  return(g)
}