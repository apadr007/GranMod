indexLookUp3 = function(df, val, spaceMax){
  dicFrame <- data.frame(val = c("A", "B", "C", "D", "E"), val2 = c("B", "A", "D", "C", "E"), stringsAsFactors = F)
  
  m = matrix(0, nrow=4, ncol=2)
  m[2,] = c(-1, 0)
  m[3,] = c(0, 1)
  m[4,] = c(0, -1)
  m = rbind(m, c(0,0) )
  rownames(m)=c('A','B','C','D','E')
  
  val2 = if(any(df >= spaceMax | df <= 1)) dicFrame[dicFrame$val == val, 'val2'] else val
  out = t(df) + m[val2, ]
  out = t(out)
  
  #val3 = if(any(out >= spaceMax | out <= 1)) dicFrame[dicFrame$val == val, 'val2'] else val2
  #out = t(out) + m[val2, ]
  #out = t(out)
  
  return(list(as.matrix(out),val2) )
}