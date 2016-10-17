#h = hash(keys=letters, values = 1:26)
#clear(h)
#rm(h)

#set up all values


row.names(layout.old) = MORELETTERS(1:nrow(layout.old) )
hx = hash(row.names(layout.old), layout.old[,1]  )
hy = hash(row.names(layout.old),layout.old[,2] )

wc = walktrap.community(g)
#wc$membership
names(wc$membership) <- 1:node_number
gran.length <- length(table(wc$membership)[table(wc$membership) >= 1])
gran.members <- lapply(1:gran.length, function(x) as.numeric(names(wc$membership[wc$membership == x])))
# make a hash from gran.members
MORELETTERS <- extend(LETTERS)
all.Letters = MORELETTERS(1:length(gran.members) )
names(gran.members) = all.Letters
gran.members = hash(gran.members)

nodeMover8 = function(layout.old, g, node_number, spaceMax, ...) {
  all.val = 1:nrow(layout.old)
  
  
  
  ##make a hash for layout.old. This requires hashing each col seperately
  #row.names(layout.old) = MORELETTERS(1:nrow(layout.old) )
  #hx = hash(row.names(layout.old), layout.old[,1]  )
  #hy = hash(row.names(layout.old),layout.old[,2] )
  
  index.m <- matrix(1:5)
  rownames(index.m) <- c("A","B","C","D", "E")
  
  out = list()
  
  #v.val = valency_func(g, node_number)
  #v.val2 = which(v.val[,2] >= 2)
  
  possibleMoves = 1:4
  
  for (move in 1:gran.length) {
    out2 = list(); truefalse = matrix(); o.avail = list(); 
    x = NA; y = NA; v.old = NA; x.old = NA; y.old = NA 
    x1 = NA; x2 = NA; y1 = NA; y2 = NA; v = NA; v3 = NA
    direction = NA; A=NA;B=NA;C=NA;D=NA; truefalse = NA
    x10.dupName = NA; idx.val3=NA
    
    #index the gran.members hash
    mover.m <- gran.members[[ all.Letters[move] ]]
    #granule movement is based on the first (i.e. smallest) node in the granule
    mover = mover.m[1]
    

    x = hx[[ all.Letters[mover] ]]
    y = hy[[ all.Letters[mover] ]]
    
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
    m <- reflectFunc(m, spaceMax)
    
    # prep the subset of layout.old / make a hash for layout.Gran
    gran.subset <- all.val[-mover.m]
    #val.hy = values(hy, MORELETTERS(gran.subset) )
    val.hy = hy[ all.Letters[gran.subset] ]
    val.hx = hx[ all.Letters[gran.subset] ]
    #val.hx = values(hx, MORELETTERS(gran.subset) )
    
    ##make a hash for layout.Gran
    #gran.hx = hash(keys=names(val.hx), values=as.numeric(val.hx) )
    #gran.hy = hash(keys=names(val.hy), values=as.numeric(val.hy) )
    
    # this makes layout.Gran
    layout.Gran = cbind( values(val.hx, USE.NAMES=FALSE), 
                          values(val.hy, USE.NAMES=FALSE) )
    row.names(layout.Gran) <- names(val.hx)
    #layout.Gran <- layout.Gran[order(layout.Gran[,3], decreasing = FALSE), ]
    #layout.Gran <- n(layout.Gran[,1:2])
    
    if (length(mover.m) > 1){ truefalse = apply(m,MARGIN = 1, layoutOverlapFinder.v, m = layout.Gran)  } else { truefalse = apply(m,MARGIN = 1, layoutOverlapFinder.v, m = layout.old)}
    
    truefalse[length(truefalse)+1] = FALSE
    names(truefalse) = c("A", "B", "C", "D", "E")
    m <- rbind(m, E)
    ## from the possible directions I can move toward, 
    ## I search for the ones that are occupied & print available
    
    x10 <- m[truefalse == FALSE, ]
    x10 <- matrix(x10, ncol=2)
    if ( length(mover.m) > 1 ) {
      
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
      #posit.move <- reflectFunc2(gran.members[[move]], posit.move, layout.old, spaceMax)
      #df.list <- indexLookUp3(df.list, posit.move, layout.old)
      # df.list <- reflectFunc(df.list, spaceMax)
      #df.list = indexLookUp(gran.members[[move]], posit.move, layout.old)
      df.list = reflectFunc3(spaceMax, posit.move, hx, hy)
      df.list_1 = df.list[[1]]
      posit.move = df.list[[2]]
      boundResp <- indexLookUp3(df.list_1, posit.move, spaceMax)
      df.list <- boundResp[[1]]
      posit.move <- boundResp[[2]]
      
      #if ( !isTRUE( apply(df.list,MARGIN = 1, layoutOverlapFinder, m = layout.Gran)) ) { print('test') }
      #if(any(apply(df.list,MARGIN = 1, layoutOverlapFinder, m = layout.Gran))) { 
      #  next
      #df.list <- layout.old[ gran.members[[move]], ] 
      #  }
      
      df.list <- lapply( seq_len(nrow(df.list)), function(i) df.list[i,] ) #converts df.list to list from matrix
      
      df.truefalse = matrix()
  
      if ( is.vector(df.list) ){
        df.truefalse[i] <- layoutOverlapFinder.v(layout.old, df.list)
      } else {  df.truefalse[i] <- sapply(layoutOverlapFinder.m(layout.old, x), v.tester)  }
      
      
      #for ( i in 1:length(df.list) ){ df.truefalse[i] = layoutOverlapFinder(layout.Gran, df.list[[i]] ) }
      
      #for ( i in 1:length(df.list) ){ df.truefalse2[i] = layoutOverlapFinder(layout.old, df.list[[i]] ) }
      ###if ( isTRUE(df.truefalse)==FALSE ) { next }
      
      if(any(df.truefalse) ) { next }
      #if(any(df.truefalse2) ) { next }
      
      df.list = do.call(rbind, df.list)
      df.list <- as.data.frame(df.list)
      df.list$V3 = posit.move
      ###testr = which(df.truefalse==FALSE)
      ###if (isEmpty(testr) == FALSE ) { 
      idx.val <- index.m[posit.move, drop=FALSE,]
      idx.val2 = sample(rownames(idx.val))[1]
      idx.val3 = df.list[df.list$V3 == idx.val2, 1:2]
      
      
     #layout.old[gran.members[[move]], ] <- as.matrix(idx.val3) ###}
      
      hx[ all.Letters[mover.m] ] <- idx.val3[, 1]
      hy[ all.Letters[mover.m] ] <- idx.val3[, 2]
      layout.old = cbind( values(hx, USE.NAMES=FALSE), 
                          values(hy, USE.NAMES=FALSE) )
      row.names(layout.old) <- names(hx)
      #layout.old <- layout.old[order(layout.old[,3], decreasing = FALSE), ]
      #layout.old <- layout.old[,1:2]
    
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
      df.list <- indexLookUp4(mover.m, posit.move, hx, hy)
      df.list <- reflectFunc(df.list, spaceMax)
      
      df.truefalse = matrix()
      
      if ( is.vector(df.list) ){
        df.truefalse[i] <- layoutOverlapFinder.v(layout.old, df.list)
      } else {  df.truefalse[i] <- sapply(layoutOverlapFinder.m(layout.old, x), v.tester)  }
      
      #for ( i in 1:length(df.list) ){ df.truefalse[i] = layoutOverlapFinder.v(layout.Gran, df.list[[i]] ) }

      
      if( any(df.truefalse)==TRUE ) {next}
      
      df.list = as.data.frame( matrix(df.list, ncol=2) )
      df.list$V3 = posit.move
      testr = which(df.truefalse==FALSE)
      if (isEmpty(testr)==FALSE ) { 
        idx.val <- index.m[posit.move, drop=FALSE,]
        idx.val2 = sample(rownames(idx.val))[1]
        idx.val3 = df.list[df.list$V3 == idx.val2,1:2]
        hx[ all.Letters[mover.m] ] <- idx.val3[, 1]
        hy[ all.Letters[mover.m] ] <- idx.val3[, 2]
        #hx[[ paste(mover.m, '', sep = '') ]] <- idx.val3[, 1]
        #hy[[ paste(mover.m, '', sep = '') ]] <- idx.val3[, 2]
        
        layout.old = cbind( values(hx, USE.NAMES=FALSE), 
                            values(hy, USE.NAMES=FALSE) )
        row.names(layout.old) <- names(hx)
        
        #layout.old = cbind( values(hx, USE.NAMES=FALSE), 
        #                     values(hy, USE.NAMES=FALSE), as.numeric(names(hx)) )    
        #layout.old <- layout.old[order(layout.old[,3], decreasing = FALSE), ]
        #layout.old <- layout.old[,1:2]
        
        #layout.old[gran.members[[move]], ] <- as.matrix(idx.val3) }
      }
    }
  }
  
  return(layout.old) 
}
clear(gran.members)
rm(gran.members)
clear(hx); clear(hy)
rm(hy); rm(hx)
clear(gran.hx); clear(gran.hy)
rm(gran.hx); rm(gran.hy)
