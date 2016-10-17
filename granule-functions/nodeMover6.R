nodeMover6 = function(layout.old, g, node_number, spaceMax){
  v = NA; v.val2 = 0; v.val= integer(0)
  out = list()
  
  v.val = valency_func(g, node_number)
  v.val2 = which(v.val[,2] >= 2)
  
  possibleMoves = 1:4
  
  for (move in 1:node_number) {
    out2 = list(); truefalse = matrix(); o.avail = list(); 
    x = NA; y = NA; v.old = NA; x.old = NA; y.old = NA; 
    x1 = NA; x2 = NA; y1 = NA; y2 = NA; v = NA; v3 = NA
    direction = NA; A=NA;B=NA;C=NA;D=NA
    
    v.old <- layout.old[move, ]
    x = layout.old[move,1]
    y = layout.old[move, 2] 
    
    x1 = x + 0.4
    x2 = x - 0.4
    y1 = y + 0.4
    y2 = y - 0.4
    
    A = c(x1, y)
    B = c(x2, y)
    C = c(x, y1)
    D = c(x, y2)
    E = c(x, y)
    
    out[[1]] = as.matrix(rbind(A, B, C, D, E))
    m <- out[[1]][sample(nrow(out[[1]])),] #randomize the positions
    m = as.matrix(m)
    
    v.val3 = integer(0)
    v.val3 = valency_func(g, node_number)[move, 2]
    
    if (v.val3 >= 1 ) { 
      
      #this sets up movement for granule containing by correcting for the amount I'm moving in the above. 
      x1 = x + 0.1
      x2 = x - 0.1
      y1 = y + 0.1
      y2 = y - 0.1
      
      A = c(x1, y)
      B = c(x2, y)
      C = c(x, y1)
      D = c(x, y2)
      E = c(x, y)
      
      out[[1]] = as.matrix(rbind(A, B, C, D, E))
      out[[1]] <- out[[1]][sample(nrow(out[[1]])),]
      m <- out[[1]]
      m = as.matrix(m)
    }
    
    m <- boundaryFunc(m, spaceMax)
    
    truefalse = apply(m,MARGIN = 1, layoutOverlapFinder, m = layout.old) #true means I will overlap with another node if I move there! 
    
    ## from the possible directions I can move toward, I search for the ones that are occupied & print available
    x10 = m[truefalse == FALSE,]
    #if (isEmpty(x10) == TRUE) { next }
    
    layout.old[move, ] = x10[1,]
    
    #dupCheck <- which(duplicated(layout.old)==TRUE )
    #if ( isEmpty(dupCheck)==FALSE ) {
    #  layout.old[move, ] <- v.old 
    #  }
    }
  return(layout.old)
  }