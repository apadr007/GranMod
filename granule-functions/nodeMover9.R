nodeMover9 = function(layout.old, g, node_number, spaceMax){
  wc = walktrap.community(g)
  wc$membership
  names(wc$membership) <- 1:node_number
  gran.length <- length(table(wc$membership)[table(wc$membership) >= 1])
  gran.members <- lapply(1:gran.length, function(x) as.numeric(names(wc$membership[wc$membership == x])))
  
  #gran.members <- get.edges(g, E(g) )
  #if ( !isEmpty(gran.members)){
  #  val <- 1:node_number
  #  gran.members <- c(lapply(seq_len(nrow(gran.members)), function(i) gran.members[i,]),  as.list( val[-gran.members] )) #two functions here. I'm cating together dimers and single nodes into a list
  #} else { gran.members <- as.list(val) }
  
  index.m <- matrix(1:5)
  rownames(index.m) <- c("A","B","C","D", "E")
  
  out = list()
  
  #v.val = valency_func(g, node_number)
  #v.val2 = which(v.val[,2] >= 2)
  
  possibleMoves = 1:4
  
  vx = 1
  vy = 1
  t = 1e-3
  sigma.x = vx*t
  sigma.y = vy*t
  
  for (move in 1:gran.length) {
    
    out2 = list(); truefalse = matrix(); o.avail = list() 
    x = NA; y = NA; v.old = NA; x.old = NA; y.old = NA 
    x1 = NA; x2 = NA; y1 = NA; y2 = NA; v = NA; v3 = NA
    direction = NA; A=NA;B=NA;C=NA;D=NA; truefalse = NA
    x10.dupName = NA; idx.val3=NA
    
    #print(move)
    #sigma.test.x <- sample(c(TRUE,FALSE), size = 1, prob = c(sigma.x, 1-sigma.x))
    #sigma.test.y <- sample(c(TRUE,FALSE), size = 1, prob = c(sigma.y, 1-sigma.y))
    
    #sigma.test = c(sigma.test.x, sigma.test.y)
    #if (all(sigma.test)==FALSE) next
    
    #val.x <- ifelse(sigma.test.x, 1, 0)
    #val.y <- ifelse(sigma.test.y, 1, 0)
    
    mover <- gran.members[[move]][1]
    v.old <- layout.old[mover, ]
    x = layout.old[mover,1]
    y = layout.old[mover,2]  
    
    x1 = x + 1
    x2 = x - 1
    y1 = y + 1
    y2 = y - 1
    
    A = c(x1, y)
    B = c(x2, y)
    C = c(x, y1)
    D = c(x, y2)
    E = c(x, y)
    
    m = as.matrix(rbind(A, B, C, D))
    
    v.val3 = valency_func(g, node_number)[mover, 2]
    
    #m <- boundaryFunc.rnp(m, spaceMax)
    #}
    m <- reflectFunc(m, spaceMax)
    #m <- boundaryFunc(m, spaceMax)
    
    #m <- unique(m)
    layout.Gran <- layout.old[-gran.members[[move]],] ## this removes members of the same community when I'm checking for overlap! 
    
    if (length(gran.members[[move]]) > 1){ truefalse = layoutOverlapFinder.m(layout.Gran, m)  } else { truefalse = layoutOverlapFinder.m(layout.old, m) }
    df.truefalse = as.vector( as.logical(truefalse) )
    df.truefalse[length(df.truefalse)+1] = FALSE
    names(df.truefalse) = c("A", "B", "C", "D", "E")
    m <- rbind(m, E)
    ## from the possible directions I can move toward, 
    ##I search for the ones that are occupied & print available
    
    x10 <- m[df.truefalse == FALSE, ]
    x10 <- matrix(x10, ncol=2)
    if ( length(gran.members[[move]]) > 1 ) {
      
      x10 <- x10[sample(nrow(x10)),]
      #print('test 1')
      if (  isTRUE(nrow(x10) > 1) && isEmpty(x10[duplicated(x10), ]) == FALSE) { x10 <- x10[duplicated(x10), ] }
      #print('test 2')
      m <- unique(m)
      x10 <- matrix(x10, ncol=2)
      x10 = data.frame(x10)
      posit.move <- apply(x10, 1, function(i) rownames(m)[!colSums(t(m) != as.vector(i))])
      posit.move = as.list(posit.move)
      posit.move <- do.call(rbind, posit.move)
      posit.move = sample(posit.move)[1]
      
      df.list = reflectFunc2(layout.old[gran.members[[move]], ], spaceMax, posit.move)
      df.list_1 = df.list[[1]]
      posit.move = df.list[[2]]
      boundResp <- indexLookUp3(df.list_1, posit.move, spaceMax)
      df.list <- boundResp[[1]]
      posit.move <- boundResp[[2]]
      
      df.list <- lapply( seq_len(nrow(df.list)), function(i) df.list[i,] ) #converts df.list to list from matrix
      
      df.truefalse = matrix()
      df.list <- do.call(rbind, df.list)
      if ( is.vector(unlist(df.list)) ) { 
        df.truefalse = layoutOverlapFinder.v(layout.Gran, unlist(df.list) ) 
      } else { 
        df.truefalse = layoutOverlapFinder.m(layout.Gran, unlist(df.list) ) 
        df.truefalse = as.vector( as.logical(df.truefalse) )
        }
      
      if(any(df.truefalse) ) { next }
      if( !is.matrix(df.list) ) { df.list = do.call(rbind, df.list) }
      
      df.list <- as.data.frame(df.list)
      df.list$V3 = posit.move
      ###testr = which(df.truefalse==FALSE)
      ###if (isEmpty(testr) == FALSE ) { 
      idx.val <- index.m[posit.move, drop=FALSE,]
      idx.val2 = sample(rownames(idx.val))[1]
      idx.val3 = df.list[df.list$V3 == idx.val2,1:2]
      layout.old[gran.members[[move]], ] <- as.matrix(idx.val3) ###}
      
    } else {
      #x10 = data.frame(x10)
      ## print('in the else if…')
      #print('2')
      
      x10 <- x10[sample(nrow(x10)),]
      x10 <- matrix(x10, ncol=2)
  
      if ( isEmpty(x10[duplicated(x10), ])==FALSE) { x10 <- x10[duplicated(x10), ] }

      posit.move <- apply(matrix(x10, ncol = 2), 1, function(i) rownames(m)[!colSums(t(m) != as.vector(i))])
      posit.move = as.list(posit.move)
      posit.move <- do.call(rbind, posit.move)
      posit.move = sample(posit.move)[1]
      df.list <- indexLookUp(gran.members[[move]], posit.move,layout.old)
      df.list <- reflectFunc(df.list, spaceMax)
      
      df.list <- lapply(seq_len(nrow(df.list)), function(i) df.list[i,]) #converts df.list to list from matrix
      df.truefalse = matrix()
      
      if ( is.vector(unlist(df.list)) ) { 
        df.truefalse = layoutOverlapFinder.v(layout.Gran, unlist(df.list) ) 
        } else { df.truefalse = layoutOverlapFinder.m(layout.Gran, unlist(df.list) ) }
      
      if( any(df.truefalse)==TRUE ) {next}
      
      df.list = do.call(rbind, df.list)
      df.list <- as.data.frame(df.list)
      df.list$V3 = posit.move
      testr = which(df.truefalse==FALSE)
      if (isEmpty(testr)==FALSE ) { 
        idx.val <- index.m[posit.move, drop=FALSE,]
        idx.val2 = sample(rownames(idx.val))[1]
        idx.val3 = df.list[df.list$V3 == idx.val2,1:2]
        layout.old[gran.members[[move]], ] <- as.matrix(idx.val3) }
    }
  }
  return(layout.old)
}