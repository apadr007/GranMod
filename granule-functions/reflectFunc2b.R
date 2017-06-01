reflectFunc2b <- function(m, spaceMax, val) {
  dicFrame <- data.frame(val = c("A", "B", "C", "D", "E"), val2 = c("B", "A", "D", "C", "E"), stringsAsFactors = F)
  g = m
  #valTest <- sample(0:1, 1)
  ## x 
  #if( any(m[,1] >= spaceMax) && valTest == 1) m[,1] <- m[,1] - 2
  #if( any(m[,1] <= 1) && valTest == 1 ) m[,1] <- m[,1] + 2
  if( any(m[,1] >= spaceMax) ) m[,1] <- m[,1] - 1
  if( any(m[,1] <= 1) ) m[,1] <- m[,1] + 1
  ## y
  #if( any(m[,2] >= spaceMax) && valTest == 0) m[,2] <- m[,2] - 2
  #if( any(m[,2] <= 1) && valTest == 0) m[,2] <- m[,2] + 2
  if( any(m[,2] >= spaceMax) ) m[,2] <- m[,2] - 1
  if( any(m[,2] <= 1) ) m[,2] <- m[,2] + 1
  
  if ( all(g == m) ) val2 = val else val2 = dicFrame[dicFrame$val == val, 'val2']
  return(list(m, val2) )
}