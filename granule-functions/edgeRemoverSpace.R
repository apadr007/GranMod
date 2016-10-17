edgeRemoverSpace= function(g, layout.old){
  for (v in 1:nrow(layout.old)){
    n = neighbors(g, v)
    f = findNearestMovers3(layout.old)[[v]]
    b = which( !n%in%f)
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
