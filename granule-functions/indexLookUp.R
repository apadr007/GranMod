indexLookUp = function(gran.members, posit.move, ...){
  out = list()
  for (i in 1:length(gran.members)) {
      
    x = layout.old[gran.members[i] ,1]
    y = layout.old[gran.members[i] ,2]  
    
    A = c(x + 1, y)
    B = c(x - 1, y)
    C = c(x, y + 1)
    D = c(x, y - 1)
    E = c(x, y)
    m = as.matrix(rbind(A, B, C, D, E))
    
    out[[i]] <- m[ posit.move,]
  }
  return(do.call(rbind, out) )
}
