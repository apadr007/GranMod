is.integer0 <- function(x){
  is.integer(x) && length(x) == 0L
}

drawErdos = function(n, p){
    g = erdos.renyi.game(n, p)
  missing <- which(degree(g) == 0)
  if (!isEmpty(missing)){
    n = 1:n
    n = n[!n == missing]
    for (i in 1:length(missing)){
      newEdge <- sample(n, size=1, replace=TRUE)
      g = add.edges(g, c(missing[i], newEdge))
    }
    suppressWarnings(return(g))
  } else {
    return(g)}
}
# toss the coins
toss.create = function(freq) {
  tossing = NULL
  for (i in 1:freq ) tossing[i] = sample(coins, 1, rep=TRUE, prob=probabilities)
  tossing = sum(tossing)
  return (tossing)
}
toss.decay = function(freq) {
  tossing = NULL
  for (i in 1:freq ) tossing[i] = sample(coins, 1, rep=TRUE, prob=probabilities.decay)
  tossing = sum(tossing)
  return (tossing)
}

# function to update the diffusers as suggested by NIT  <- the approach I tend to favor... 
update_diffusers = function(diffusers){
  nearest_neighbors = data.frame(table(unlist(neighborhood(g, 1, diffusers))))
  nearest_neighbors = subset(nearest_neighbors, !(nearest_neighbors[,1]%in%diffusers))
  x = data.frame(table(unlist(diffusers)))
  
  keep = unlist(lapply(nearest_neighbors[,2], toss.create))
  loss = unlist(lapply(x[,2], toss.decay))
  
  new_infected = as.numeric(as.character(nearest_neighbors[,1][keep >= 1]))
  loss_infected = as.numeric(as.character(x[,1][loss >= 1])) #getting decayed nodes
  diffusers = diffusers[!diffusers %in% loss_infected] #decaying the nodes 
  diffusers = unique(c(diffusers, new_infected))
  return(diffusers)
}

### option 2 of 2 for the sponteanous addons as suggested by lab ###
update_diffusers = function(diffusers){
  myNeighbors = c(1,2,3)                     # the neighborhood repeats if its set to order=3  <- redo this. 
  output = list()
  set1 = neighborhood(g, 1, diffusers)
  set2 = neighborhood(g, 2, diffusers)
  set3 = neighborhood(g, 3, diffusers)
  
  delete2 = which(set1[[1]]%in%set2[[1]])
  set2 = set2[[1]][-delete2]
  delete3 = which(set2[[1]]%in%set3[[1]])
  set3 = set3[[1]][-delete3]
  
  sets = list()
  sets[[1]] <- set1
  sets[[2]] <- set2
  sets[[3]] <- set3
  
  for (i in 1:3){
    nearest_neighbors = data.frame(table(unlist(sets[[i]])))
    nearest_neighbors = subset(nearest_neighbors, !(nearest_neighbors[,1]%in%diffusers))  
    keep = unlist(lapply(nearest_neighbors[,2], toss, prob_neighbor=c(probabilities[i,1], probabilities[i,2])))
    new_infected = as.numeric(as.character(nearest_neighbors[,1][keep >= 1]))
    diffusers = unique(c(diffusers, new_infected))
    diffusers 
    return(diffusers)
  }
}

# plot the results:
plot_time_series = function(infected, m){
  num_cum = unlist(lapply(1:m, 
                          function(x) length(infected[[x]]) ))
  p_cum = num_cum/node_number
  p = diff(c(0, p_cum))
  time = 1:m
  plot(p_cum~time, type = "b", 
       ylab = "Fraction of interacting mRNPs", xlab = "n",
       xlim = c(0,total_time), ylim =c(0,1))
}

#this returns the cumulative probabilities without plotting them
plot_time_series2 = function(infected, m){
  num_cum = unlist(lapply(1:m,
                          function(x) length(infected[[x]]) ))
  p_cum = num_cum/node_number
  p = diff(c(0, p_cum))
  return(p_cum)
}

#this returns a hash table of current valencies for a graph
valency_func = function(g, node_number){  ## integrate this function into edgeAdder with conditionals
  hash_list = NULL; black_list = NULL
  list = NULL; blacks = NULL;list2 = NULL
  hash_list = cbind(1:node_number, 0)
  
  list <- as.vector(get.edges(g, E(g)[which(E(g)$color=='black')]))
  for (j in 1:length(list)){
    list2 = list[j]
    hash_list[list2,2] <- hash_list[list2,2] + 1
  }
  return(hash_list)
}

#this returns a probabilistic true or false statement based on the valency of each node -- its input is the hash table from valency_func
probFunc = function(x, highestProb){
  possibilities = c('TRUE','FALSE')
  val = NULL
  
  if(x == 1){
    val <- sample(possibilities, size = 1, prob = c(highestProb/6, (1-highestProb/6)))
  } 
  else if (x == 2){
    val <- sample(possibilities, size = 1, prob = c(highestProb/5, (1-highestProb/5)))
  }
  else if (x == 3){
    val <- sample(possibilities, size = 1, prob = c(highestProb/4, (1-highestProb/4)))
  }
  else if (x == 4){
    val <- sample(possibilities, size = 1, prob = c(highestProb/3, (1-highestProb/3)))
  }
  else if (x == 5){
    val <- sample(possibilities, size = 1, prob = c(highestProb/2, (1-highestProb/2)))
  }
  else if (x == 6){
    val <- sample(possibilities, size = 1, prob = c(highestProb, (1-highestProb)))
  }
  else { 
    return("FALSE") 
  }
  return(val)
}
probFunc = function(x, highestProb){
  possibilities = c('TRUE','FALSE')
  val = NULL
  
  if(x == 1){
    val <- sample(possibilities, size = 1, prob = c(highestProb, (1-highestProb)))
  } 
  else if (x == 2){
    val <- sample(possibilities, size = 1, prob = c(highestProb, (1-highestProb)))
  }
  else if (x == 3){
    val <- sample(possibilities, size = 1, prob = c(highestProb, (1-highestProb)))
  }
  else if (x == 4){
    val <- sample(possibilities, size = 1, prob = c(highestProb, (1-highestProb)))
  }
  else if (x == 5){
    val <- sample(possibilities, size = 1, prob = c(highestProb, (1-highestProb)))
  }
  else if (x == 6){
    val <- sample(possibilities, size = 1, prob = c(highestProb, (1-highestProb)))
  }
  else { 
    return("FALSE") 
  }
  return(val)
}
probFunc = function(x, highestProb){
  possibilities = c('TRUE','FALSE')
  val = NULL
  
  if(x == 1){
    val <- sample(possibilities, size = 1, prob = c(highestProb, (1-highestProb)))
  } 
  else if (x == 2){
    val <- sample(possibilities, size = 1, prob = c(highestProb/2, (1-highestProb/2)))
  }
  else if (x == 3){
    val <- sample(possibilities, size = 1, prob = c(highestProb/3, (1-highestProb/3)))
  }
  else if (x == 4){
    val <- sample(possibilities, size = 1, prob = c(highestProb/4, (1-highestProb/4)))
  }
  else if (x == 5){
    val <- sample(possibilities, size = 1, prob = c(highestProb/5, (1-highestProb/5)))
  }
  else if (x == 6){
    val <- sample(possibilities, size = 1, prob = c(highestProb/6, (1-highestProb/6)))
  }
  else { 
    return("FALSE") 
  }
  return(val)
}

# add edges 
edgeAdder3 = function(g, edgeA, edgeB.vector, infectors, node_number){
  edge_init = NULL; edge_inst = NULL
  probFilter = NULL; valency.out = NULL
  edge_init = which(E(g)$color == 'black')
  edgeB.vector = unlist(edgeB.vector)
  
  valency.out <- valency_func(g, node_number)
  valency.out <- valency.out[edgeB.vector, ]        # subset and match order of edgeB.vector 
  probFilter <- sapply(valency.out[,2], probFunc, highestProb = transmission_rate)   # apply probFunc() based on each nodes valency
  
  edgeB.vector_pass <- edgeB.vector[probFilter == TRUE]  ## remove the neighbors that are false in valency_func
  
  if (isTRUE(!is.numeric(edgeA)) | isTRUE(!is.numeric(edgeB.vector))) {
    return(edge_init)
  } else{
    edgeB.vector = unlist(edgeB.vector)
    edgeB.vector_filtered = edgeB.vector%in%infectors
    edgeB.vector_pass = NULL
    for (i in 1:length(edgeB.vector_filtered)){
      if(isTRUE(edgeB.vector_filtered[i])){
        edgeB.vector_pass[i] <- edgeB.vector[i] }  
      edgeB.vector_pass <- as.vector(edgeB.vector_pass)}
    edgeB.vector_pass[1] <- NA
    edgeB.vector_pass <- edgeB.vector_pass[!is.na(edgeB.vector_pass)]
    if(!isEmpty(edgeB.vector_pass)){
      for (i in 1:length(edgeB.vector_pass)){    
        E(g, c(edgeA, c(edgeB.vector_pass[i])))$color = 'black'
      }  
      edge_inst <- which(E(g)$color == 'black')
      edge_all <- c(edge_init, edge_inst)
      edge_all <- unique(edge_all)
      return(edge_all) 
    }
  } 
}

## the original add edges
edgeAdder = function(g, edgeA, edgeB.vector, infectors){
  edge_init = NULL
  edge_init = which(E(g)$color == 'black')
  edge_inst = NULL
  edgeB.vector = unlist(edgeB.vector)
  
  if (isTRUE(!is.numeric(edgeA)) | isTRUE(!is.numeric(edgeB.vector))) {
    return(edge_init)
  } else{
    edgeB.vector = unlist(edgeB.vector)
    edgeB.vector_filtered = edgeB.vector%in%infectors
    edgeB.vector_pass = NULL
    for (i in 1:length(edgeB.vector_filtered)){
      if(isTRUE(edgeB.vector_filtered[i])){
        edgeB.vector_pass[i] <- edgeB.vector[i] }  
      edgeB.vector_pass <- as.vector(edgeB.vector_pass)}
    edgeB.vector_pass[1] <- NA
    edgeB.vector_pass <- edgeB.vector_pass[!is.na(edgeB.vector_pass)]
    if(!isEmpty(edgeB.vector_pass)){
      for (i in 1:length(edgeB.vector_pass)){    
        E(g, c(edgeA, c(edgeB.vector_pass[i])))$color = 'black'
      }  
      edge_inst <- which(E(g)$color == 'black')
      edge_all <- c(edge_init, edge_inst)
      edge_all <- unique(edge_all)
      return(edge_all) 
    }
  } 
}

edgeFinder = function(g, v1,v2){
  if (g[v1,v2] == 1){
    return(TRUE)  } 
  else { return(FALSE) }
}

#valencyMax = function(g, maxVal, node_number){
  neighbor_nodes = NULL; delete_edge = NULL; delete_these = NULL
  total_nodes = 1:node_number
  
  for (i in 1:node_number){
    neighbor_nodes <- unlist(neighborhood(g, order = 1, nodes = total_nodes[i] ))
    neighbor_nodes[1] <- NA
    neighbor_nodes <- neighbor_nodes[!is.na(neighbor_nodes)]
    neighbor_nodes <- sample(neighbor_nodes, maxVal)
    
    delete_edge <- which(total_nodes%in%neighbor_nodes == FALSE)
    delete_edge <- which(delete_edge%in%total_nodes[i] == FALSE)
    
    for (j in 1:length(delete_edge)){
      #print(i); print(delete_edge)
      if (isTRUE(edgeFinder(g, total_nodes[i], delete_edge[j]))){
        g = delete.edges(g, E(g, P=c(i, delete_edge[j])))
      }
    }
  }
  return(g)
}

isEmpty = function(x) {
  return(length(x)==0)
}

#valencyMax_attempt2 = function(g, maxVal, node_number){
  neighbor_nodes = NULL; delete_edge = NULL; nodeDegree = NULL
  delete_these = NULL ; edgesToDelete = NULL; total_nodes = NULL
  total_nodes = 1:node_number
  
  x = sample(total_nodes, node_number) # randomly starting at some vertice
  for ( i in 1:length(x) ) {           # main for loop going through every vertice, randomly selecting starting position, though... 
    nodeDegree[i] <- degree(g, v = x[i])  # getting degree of each vertice in graph
   if (nodeDegree[i] > maxVal) {          # conditional statement if vertice has larger than maxVal edges
     
     neighbor_nodes <- unlist(neighborhood(g, order = 1, nodes = x[i] ))
     neighbor_nodes[1] <- NA
     neighbor_nodes <- neighbor_nodes[!is.na(neighbor_nodes)]
     
     edgesToDelete <- nodeDegree - maxVal # how many edges to delete from each vertice in order to get too the upperbound of maxVal
     neighbor_nodes <- sample(neighbor_nodes, edgesToDelete, replace = TRUE) # randomly pick the difference of edges to remove
     
     for (j in 1:length(neighbor_nodes)){   # starting loop for the adjacent edges I'll be removing. 
       if (isTRUE(edgeFinder(g, x[i], neighbor_nodes[j]))){ # conditional 
         g = delete.edges(g, (E(g, P=c(x[i], neighbor_nodes[j]))))
       }
     } 
   } 
  }
   return(g)
}


removeOne = function(x){
  x[1] <- NA
  x <- x[!is.na(x)]
  return(x)
}

## this function will produce a list of active edge indices
#given a list of active vertices (mRNPs) and a 
activeEdges = function(g, infected){
  m=1
  edge_inst2 = NULL; edge_inst = NULL; x = NULL
  #layout.old = layout.fruchterman.reingold(g)
  E(g)$color = 'white'; V(g)$name = ''
  valency_out = NULL
  while(m <= length(infected)){
    V(g)$color = "white"
    V(g)$color[V(g)%in%infected[[m]]] = "red"
    #adding edges
    for (j in 1:length(infected[[m]])){
      neighbors = NULL; infectors = NULL
      neighbors <- neighborhood(graph = g, order = 1, nodes = infected[[m]][j])
      neighbors <- as.vector(unlist(neighbors))
      infectors <- as.vector(infected[[m]])
      edgeA <- as.vector(infected[[m]][j])
      x[[m]] = unique(c(x[[m-1]], edgeAdder3(g, edgeA, neighbors, infectors, node_number))) ### wtf did I write here....
      edge.output <-  edgeAdder3(g, edgeA, neighbors, infectors, node_number) ### this repeats in the line above...
      E(g)$color[edge.output] = 'black'
      valency_out[[m]] <- valency_func(g, node_number)
    }
    #plot(g, layout =layout.old)
    m = m + 1
  }
  #return(list(x,valency_out))
  return(x)
}      


activeGraph = function(g, infected, drawActiveGraph){
  m=1
  edge_inst2 = NULL; edge_inst = NULL; x = NULL
  #layout.old = layout.fruchterman.reingold(g)
  E(g)$color = 'white'; V(g)$name = ''
  while(m <= length(infected)){
    V(g)$color = "white"
    V(g)$color[V(g)%in%infected[[m]]] = "red"
    #adding edges
    for (j in 1:length(infected[[m]])){
      neighbors = NULL; infectors = NULL
      neighbors <- neighborhood(graph = g, order = 1, nodes = infected[[m]][j])
      neighbors <- as.vector(unlist(neighbors))
      infectors <- as.vector(infected[[m]])
      edgeA <- as.vector(infected[[m]][j])
      edge.output <-  edgeAdder3(g, edgeA, neighbors, infectors, node_number)
      E(g)$color[edge.output] = 'black'
    }
    if ( drawActiveGraph == TRUE){
      plot(g, layout =layout.old)
      m = m + 1
    } else { 
      m = m + 1
    }
  }
  return(g)
}

activeGraph.step = function(g, infected, activeEdge, n, node_number, edge.OR.node){
  edge_inst2 = NULL; edge_inst = NULL; x = NULL
  E(g)$color = 'black'; V(g)$name = ''
  total_node = 1:node_number
  totalEdge = nrow(get.edges(g, E(g)))
  totalEdge = 1:totalEdge
  activeEdge <- is.integer0(activeEdge)
  
  active_edge = activeEdge[[n]]
  active_node = infected[[n]]
  
  delete_edge <- which(totalEdge%in%active_edge == FALSE)
  delete_node <- which(total_node%in%active_node == FALSE)
  
  if (edge.OR.node == "edge"){
    g = delete.edges(g, delete_edge)
    return(g)
  }
  else if (edge.OR.node == 'node'){
    g = delete.vertices(g, delete_node)
    return(g)
  }
}
 
is.integer0 <- function(x) {
  for(i in 1:length(x)){  
    if (isTRUE(is.integer(x[[i]])) & length(x[[i]]) == 0L) {
      x[[i]] <- 0 
    }
  }
  return(x)
}  ### this is a modified in.interger0() for edge functions... 

density.graph = function(g, node_number){
  x = (node_number*(node_number-1))/2
  current.edges = nrow(as.matrix(E(g)))
  density=current.edges/x
  return(density)
}


edgeAdder4 = function(g, edgeA, edgeB.vector, infectors, node_number){
  edge_init = NULL; edge_inst = NULL
  probFilter = NULL; valency.out = NULL
  edge_init = which(E(g)$color == 'black')
  edgeB.vector = unlist(edgeB.vector)
  
  valency.out <- valency_func(g, node_number)
  valency.out <- valency.out[edgeB.vector, ]        # subset and match order of edgeB.vector 
  probFilter <- sapply(valency.out[,2], probFunc, highestProb = 0.128)   # apply probFunc() based on each nodes valency
  
  edgeB.vector_pass <- edgeB.vector[probFilter == TRUE]  ## remove the neighbors that are false in valency_func
  
  if (isTRUE(!is.numeric(edgeA)) | isTRUE(!is.numeric(edgeB.vector))) {
    return(edge_init)
  } else{
    edgeB.vector = unlist(edgeB.vector)
    edgeB.vector_filtered = edgeB.vector%in%infectors
    edgeB.vector_pass = NULL
    for (i in 1:length(edgeB.vector_filtered)){
      if(isTRUE(edgeB.vector_filtered[i])){
        edgeB.vector_pass[i] <- edgeB.vector[i] }  
      edgeB.vector_pass <- as.vector(edgeB.vector_pass)}
    edgeB.vector_pass[1] <- NA
    edgeB.vector_pass <- edgeB.vector_pass[!is.na(edgeB.vector_pass)]
    if(!isEmpty(edgeB.vector_pass)){
      for (i in 1:length(edgeB.vector_pass)){    
        E(g, c(edgeA, c(edgeB.vector_pass[i])))$color = 'black'
      }  
      edge_inst <- which(E(g)$color == 'black')
      edge_all <- c(edge_init, edge_inst)
      edge_all <- unique(edge_all)
      return(edge_all) 
    }
  } 
}
  
del_edge = function(g, valToDel){
  delme = NULL; delme.edge = NULL
  delme = neighborhood(g, 1, nodes = valToDel)
  delme[[1]][1] = NA
  delme = delme[[1]][!is.na(delme[[1]])]
  delme.edge = sample(delme, size = 1)
  g = delete.edges(g, E(g, P=c(valToDel, delme.edge)) )
  return(g)
}

rxnrate=function(t,c,parms){
  # rate constant passed through a list called parms
  k1=parms$k1
  k2=parms$k2
  k1_off=parms$k1_off
  k2_off=parms$k2_off
  # c is the concentration of species
  
  # derivatives dc/dt are computed below
  r=rep(0,length(c))
  r[1]=-k1*c["A"] + k1_off*c["B"] #dcA/dt
  r[2]=k1*c["A"] + k2_off*c["C"] - k2*c["B"]^2 #dcB/dt
  r[3]=k2*c["B"] - k2_off*c["C"] #dcC/dt
   
  # the computed derivatives are returned as a list
  # order of derivatives needs to be the same as the order of species in c
  return(list(r))
}


nodeMover = function(layout.old, g, node_number, spaceMax){
  v = valency_func(g, node_number)
  v2 = which(v[,2] >= 1)
  
  for (move in 1:node_number) {
    #    if (isEmpty(intersect(move, v2))==FALSE){ ##this is if I want to move slower as a granule forms…implement later
    x = sample(c(-1,1), size = 1)
    y = sample(c(-1,1), size = 1)
    layout.old[move,1] = layout.old[move,1] + x
    layout.old[move,2] = layout.old[move,2] + y    
    #   } else {
    #    x = sample(c(-0.01,0.01), size = 1)
    #    y = sample(c(-0.01,0.01), size = 1)
    #   layout.old[move,1] = layout.old[move,1] + x
    #    layout.old[move,2] = layout.old[move,2] + y
    if (layout.old[move,1] > spaceMax ) { ##this sets up my boundaries for the x axis
      layout.old[move, 1] = spaceMax }
    else if (layout.old[move,1] < 1 ) {
      layout.old[move, 1] = 1 }
    
    if (layout.old[move,2] > spaceMax ) { ##this sets up my boundaries for the y axis
      layout.old[move, 2] = spaceMax }
    else if (layout.old[move,2] < 1 ) {
      layout.old[move, 2] = 1 }
  }
  return(layout.old)
}

findNearestMovers2 = function(t){ ## I can think about adding a size criteria to adjust as a parameter the size of my nodes and the "neighboring range" I'm consider
  x=matrix(); x2=list()
  y = matrix(); y2 = list()
  for (i in 1:nrow(t)){
    for (j in 1:nrow(t)){
      if (abs(t[i,1] - t[j,1]) <= 0.005) {
        x[j] = j
      } else { x[j] = NA }
      if (abs(t[i,2] - t[j,2]) <= 0.005) {
        y[j] = j
      } else { y[j] = NA }
    }
    x2[[i]] = x
    y2[[i]] = y
  }
  for (i in 1:length(x2)){
    x2[[i]] = x2[[i]][!x2[[i]] == i]
    y2[[i]] = y2[[i]][!y2[[i]] == i]
  }
  x2 = lapply(x2, function(x) x[!is.na(x)])
  y2 = lapply(y2, function(x) x[!is.na(x)])
  
  #this intersects the x and y axis to find factors that are close in both axes
  z = list()
  for (i in 1:length(x2)){
    z[[i]] = intersect(unlist(x2[[i]]), unlist(y2[[i]]))
  }
  return(z)
}

cleaner.valency2 = function(nearest.mRNA_list, P.int.on){
  output = list()
  
  for (i in 1:length(nearest.mRNA_list)){
    
    filter <- sample(0:1, 1, prob = c(1-P.int.on, P.int.on) )
    if (filter == 1) {
      valency.node = sample(valency, size = 1, replace = TRUE, 
                                     prob = c(P.int.on/1, 
                                     P.int.on/2, P.int.on/3, 
                                     P.int.on/4, P.int.on/5, P.int.on/6))
      
      new.edge = nearest.mRNA_list[[i]]
      
      #valTest = valency_func(g, node_number)
      #valTestdel = which(valTest[,2] > 6)   
      #del = NA
      
      #if (isEmpty(valTestdel)==FALSE) {
      #  del = which(new.edge%in%valTestdel)
      #  if(isEmpty(del)== FALSE) { 
      #    new.edge = new.edge[-del]
      #  }
      #}
      if (isEmpty(new.edge) == TRUE) {
        output[[i]] = new.edge
      } else if (valency.node > length(new.edge)){
        output[[i]] = new.edge
      } else if (length(new.edge) == 1) {
        output[[i]] = new.edge
      } else if (valency.node <= length(new.edge)){
        new.edgeToForm <- new.edge[1:valency.node]
        output[[i]] = new.edgeToForm }
    } else { output[[i]] <- 0 }
  }
  return(output)
}

#rates = function(rate, k1, k1.rev, k2, k2.rev, t, A.0, B.0, C.0){
  if(rate == 1){
    A = A.0*exp(-k1*t) 
    B = B.0*exp(-k1.rev*t)
    out = list(A=A, B=B)
  }
  else if (rate ==2){
    B = B.0/(1 + k2*t*B.0)
    C = C.0*exp(-k2.rev*t) 
    out = list(B=B, C=C)
  }
  return(out)
}

layoutGen = function(node_number, spaceMax){
  node.position = matrix(0, nrow = node_number, ncol = 2)
  gridz = matrix(0, nrow = spaceMax, ncol= spaceMax) 
  gridz.idx = which(gridz == 0, arr.ind = T)
  
  
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
  
  
  return(gridz.true)
}


layoutOverlapFinder = function(m, x){
  #m is layout.old matrix
  #x is a 2D position of the node you are querying about moving
  val = NA
  #if ( isTRUE(is.na(x)) ){ 
  #    FALSE }
  val <- any(apply(m,1,function(n,x) all(n==x),x=x))
  if (is.na(val)==TRUE){ val <- FALSE}
  return(val)
}

edgeRemoverSpace= function(g, layout.old){
  for (v in 1:nrow(layout.old)){
    n = neighbors(g, v)
    f = findNearestMovers2(layout.old)[[v]]
    b = which( !n%in%f)
    if (isEmpty(b)==FALSE){
      for ( i in 1:length(b) ) {
        b.idx = b[i]
        n.rdy = n[ b.idx ]
        g = delete.edges(g, E(g, P=c(v, n.rdy)) )
      }
    }
  }
  return(g)
}


###this is the updated nodeMover to use <- #3 ! :D 
nodeMover3 = function(layout.old, g, node_number, spaceMax){
  v = NA; v.val2 = 0
  v.val = valency_func(g, node_number)
  v.val2 = which(v.val[,2] >= 1)
  out = list()
  possibleMoves = 1:4
  granuleNodes = layout.old
  
  for (move in 1:node_number) {
    out = list(); truefalse = matrix(); o.avail = list(); 
    x = NA; y = NA; v.old = NA; x.old = NA; y.old = NA; 
    x1 = NA; x2 = NA; y1 = NA; y2 = NA; v = NA; v3 = NA
    direction = NA; A=NA;B=NA;C=NA;D=NA
    
    
    v.old <- layout.old[move, ]
    x = layout.old[move,1]
    y = layout.old[move, 2] 
    
    x1 = x + 1
    x2 = x - 1
    y1 = y + 1
    y2 = y - 1
    
    A = c(x1, y)
    B = c(x2, y)
    C = c(x, y1)
    D = c(x, y2) 
    
    
    
    
    out[[1]] = as.matrix(rbind(A, B, C, D))
    out[[1]] <- out[[1]][sample(nrow(out[[1]])),] #randomize the positions
    ##
    direction = sample(1:4, 1) # randomly select index for out <- out has directions to move toward
    #this is promising#
    
    v.val = integer(0)
    v.val = valency_func(g, node_number)[move, 2]
    granuleNodes = layout.old
    if ( isEmpty(v.val2) == FALSE && v.val >= 1 ) { next } else {
      
      ##this tests for valency. if it's connected, it will not move for now…
      #for (i in 1:length(v.val2) ){ 
      #v3 <- v.val2[i]
      #layout.old[v3,] <- as.numeric(granuleNodes[i,])
      
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
  }
  return(layout.old)
}


communityfast.algo = function(g){
  
  mRNA.return = length(which(V(g)$color == 'lightblue'))
  fg= fastgreedy.community(g, modularity = TRUE, membership = TRUE)
  V(g)$size = 1
  E(g)$count = 1
  comm.graph <- contract.vertices(g, fg$membership, vertex.attr.comb=list(size="sum", "ignore"))
  which(vertex.attributes(comm.graph)$size>=1)
  
  for(i in 1:length(which(vertex.attributes(comm.graph)$size >= 1)) ){
    V(comm.graph)$color[i] = 'orange' 
  }
  
  E(comm.graph)$color = 'white'
  comm.graph = add.vertices(comm.graph, mRNA.return)
  V(comm.graph)$color[which(V(comm.graph)%in%which(V(comm.graph)$color == 'orange') == FALSE)] = 'lightblue'
  V(comm.graph)$size[which(V(comm.graph)%in%which(V(comm.graph)$size >= 1) == FALSE)] = 1
  
  #x = communityfast.algo(g)[ which(communityfast.algo(g) >= 3) ]
  x = vertex.attributes(comm.graph)$size
  return(x)
}