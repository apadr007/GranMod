##need to make a df.list inside reflectFunc3 so as to not modify the actual values of layout.old!

reflectFunc3 = function(spaceMax, val, ...){
  #create a layout to modify in this function
  layout.faux <- cbind( values(hx, all.Letters[mover.m], USE.NAMES=FALSE), values(hy, all.Letters[mover.m]) )
  #order by increasing order
  layout.faux <- layout.faux[ order(as.numeric( dimnames(layout.faux)[[1]] ),decreasing = FALSE), ]  
  #remember original position
  layout.faux.orig <- layout.faux
  
  dicFrame <- data.frame(val = c("A", "B", "C", "D", "E"), val2 = c("B", "A", "D", "C", "E"), stringsAsFactors = F)
  
  ## x 
  if( any(layout.faux[,1] >= spaceMax) ) layout.faux[,1] <- layout.faux[,1] - 2
  if( any(layout.faux[,1] <= 1) ) layout.faux[,1] <- layout.faux[,1] + 2
  ## y
  if( any(layout.faux[,2] >= spaceMax) ) layout.faux[,2] <- layout.faux[,2] - 2
  if( any(layout.faux[,2] <= 1) ) layout.faux[,2] <- layout.faux[,2] + 2
  
  if ( all(layout.faux.orig == layout.faux) ) val2 = val else val2 = dicFrame[dicFrame$val == val, 'val2']
  
  return(list(layout.faux, val2) )
  }