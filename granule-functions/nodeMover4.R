nodeMover4 = function(layout.old, g, node_number, spaceMax){
  v = NA; v.val2 = 0
  v.val = valency_func(g, node_number)
  v.val2 = which(v.val[,2] >= 2)
  out = list()
  possibleMoves = 1:4
  granuleNodes = layout.old
  
  for (move in 1:node_number) {
    out = list(); out2 = list(); truefalse = matrix(); o.avail = list(); 
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
    
    out[[1]] = as.matrix(rbind(A, B, C, D))
    out[[1]] <- out[[1]][sample(nrow(out[[1]])),] #randomize the positions
    
    direction = sample(1:4, 1) # randomly select index for out <- out has directions to move toward
    #this is promising#
    
    v.val = integer(0)
    v.val = valency_func(g, node_number)[move, 2]
    granuleNodes = layout.old
    if ( isEmpty(v.val2) == FALSE && v.val >= 1 ) { 
      
      #this sets up movement for granule containing by correcting for the amount I'm moving in the above. 
      #v.old <- layout.old[move, ]
      #x = layout.old[move,1]
      #y = layout.old[move, 2] 
      
      x1 = x + 0.1
      x2 = x - 0.1
      y1 = y + 0.1
      y2 = y - 0.1
      
      A = c(x1, y)
      B = c(x2, y)
      C = c(x, y1)
      D = c(x, y2) 
      
      out[[1]] = as.matrix(rbind(A, B, C, D))
      out[[1]] <- out[[1]][sample(nrow(out[[1]])),]
      }
      
      for (i in 1:nrow(out[[1]]) ) {
        oldposit = out[[1]][i, ]
        
        if (out[[1]][i,1] < 0) {
          out[[1]][i, ] = c(0, oldposit[2] )
        }
        else if (out[[1]][i,1] > spaceMax){
          out[[1]][i, ] = c(spaceMax, oldposit[2] )
        }
        else if (out[[1]][i, 2] < 0){
          out[[1]][i, ] = c(oldposit[1], 0)
        }
        else if (out[[1]][i, 2] > spaceMax){
          out[[1]][i, ] = c(oldposit[1], spaceMax )
        } else { out[[1]][i, ] <- oldposit }
        
        v = out[[1]][i,]
        
        if (is.na(v[i]) == FALSE){
          truefalse[i] = layoutOverlapFinder(layout.old, v) 
        } 
        else if (is.na(v[i] == TRUE)) {
          
      }
        
      }
      
       ##else {
      
      ##this tests for valency. if it's connected, it will not move for now…
      #for (i in 1:length(v.val2) ){ 
      #v3 <- v.val2[i]
      #layout.old[v3,] <- as.numeric(granuleNodes[i,])
      
     ## for (i in 1:nrow(out[[1]]) ) {
      ##  oldposit = out[[1]][i, ]
        
      ##  if (out[[1]][i,1] < 0) {
      ##    out[[1]][i, ] = c(0, oldposit[2] )
      ##  }
      ##  else if (out[[1]][i,1] > spaceMax){
      ##    out[[1]][i, ] = c(spaceMax, oldposit[2] )
      ##  }
      ##  else if (out[[1]][i, 2] < 0){
      ##    out[[1]][i, ] = c(oldposit[1], 0)
      ##  }
      ##  else if (out[[1]][i, 2] > spaceMax){
      ##    out[[1]][i, ] = c(oldposit[1], spaceMax )
      ##  } else { out[[1]][i, ] <- oldposit }
        
      ##  v = out[[1]][i,]
        
      ##  if (is.na(v[i]) == FALSE){
      ##    truefalse[i] = layoutOverlapFinder(layout.old, v) 
      ##  } 
      ##  else if (is.na(v[i] == TRUE)) {
          
      ##  }
        
      #}
      
      
      ## from the possible directions I can move toward, 
      ##I search for the ones are occupied & print available
      x10 = which(FALSE%in%truefalse)
      if ( !isTRUE(isEmpty(x10)) ) { 
        
        for (i in 1:length(x10)){
          o.avail[[i]] <- as.numeric(out[[1]][x10[i], ])        
        }
        dirMove <- sample(1:length(o.avail), 1)
        v = as.numeric(o.avail[[dirMove]])
        
        layout.old[move, ] = v
        dupCheck <- which(duplicated(layout.old)==TRUE )
        if ( isEmpty(dupCheck)==FALSE ) {
          
          layout.old[move, ] <- v.old }
      } else { 
        
        layout.old[move, ] = v.old }
      
      
      #  else if (layout.old[move,1] < 1 ) {
      #        layout.old[move, 1] = 1 }
      #  if (layout.old[move,2] > spaceMax ) { ##this sets up my boundaries for the y axis
      #       layout.old[move, 2] = spaceMax }
      #  else if (layout.old[move,2] < 1 ) {
      #        layout.old[move, 2] = 1 }
    } 
  return(layout.old)
}
