indexLookUp2 = function(gran.members, val, layout.old, spaceMax){
  dicFrame <- data.frame(val = c("A", "B", "C", "D", "E"), val2 = c("B", "A", "D", "C", "E"), stringsAsFactors = F)
  
  m = matrix(0, nrow=4, ncol=2)
  m[1,] = c(1, 0)
  m[2,] = c(-1, 0)
  m[3,] = c(0, 1)
  m[4,] = c(0, -1)
  m = rbind(m, c(0,0) )
  rownames(m)=c('A','B','C','D','E')

  val2 = if(any(layout.old[gran.members, ] >= spaceMax | layout.old[gran.members, ] <= 1)) dicFrame[dicFrame$val == val, 'val2'] else val
  out = t(layout.old[gran.members, ]) + m[val2, ]
  out = t(out)
  
  
  return(list(as.matrix(out), val2) )
}