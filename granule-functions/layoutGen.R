layoutGen = function(node_number, spaceMax){
  node.position = matrix(0, nrow = node_number, ncol = 2)
  gridz = matrix(0, nrow = spaceMax, ncol= spaceMax) 
  gridz.idx = which(gridz == 0, arr.ind = TRUE)
  
  
  x = 0
  samplesize = node_number + 1
  while (length(x) < samplesize){
    x1 = sample(1:nrow(gridz.idx), size = node_number)
    x = c(x, x1)
    x = unique(x)
  }
  x = x[-1]
  
  for (i in x){
    selector = gridz.idx[i,]
    gridz[selector[1], selector[2] ] = i
  }
  
  gridz.true = which(gridz > 0, arr.ind = TRUE)
  gridz.true = cbind(x=gridz.true[,1], y=gridz.true[,2], node=x)
  
  #gridz = t(gridz[1:spaceMax,spaceMax:1])
  
  
  return(list(gridz.true,gridz))
}