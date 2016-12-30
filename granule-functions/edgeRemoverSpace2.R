g2=g

edgeRemoverSpace2 = function(g, layout.old, nearest.mRNA3) {
  v = 1:nrow(layout.old)
  n = sapply(v, neighbors, g=g )  
  b = !n%in%nearest.mRNA3
  
  for ( i in 1:length(b) ) {
    if ( b[i] == TRUE ) {
      b.idx = b[i]
      n.rdy = n[ b.idx ]
      g = delete.edges(g, E(g, P=c(v[i], n.rdy)) )
      }
    }
    return(g)
}

edgeRemoverSpace <- function(g, layout.old){
  near = findNearestMovers4(layout.old)
  for (v in 1:nrow(layout.old) ) {
    n = neighbors(g, v)
    f = near[[v]]
    b = which( !n%in%f )
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

edgeRemoverSpace2 <- function(g, layout.old, near){
  n = sapply(v, neighbors, g=g)
  
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