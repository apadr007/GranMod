drawErdos = function(n, p){
    g = erdos.renyi.game(n, p)
  missing <- which(degree(g) == 0)
  if (!isEmpty(missing)){
    n = 1:n
    n = n[!n == missing]
    for (i in 1:length(missing)){
      newEdge <- sample(n, size=1, replace=T)
      g = add.edges(g, c(missing[i], newEdge))
    }
    suppressWarnings(return(g))
  } else {
    return(g)}
}
# toss the coins
toss = function(freq, prob_neighbor) {
  tossing = NULL
  for (i in 1:freq ) tossing[i] = sample(coins, 1, rep=TRUE, prob=prob_neighbor)
  tossing = sum(tossing)
  return (tossing)
}
# function to update the diffusers as suggested by NIT  <- the approach I tend to favor... 
update_diffusers = function(diffusers){
  nearest_neighbors = data.frame(table(unlist(neighborhood(g, 1, diffusers))))
  nearest_neighbors = subset(nearest_neighbors, !(nearest_neighbors[,1]%in%diffusers))
  # if (isTRUE(nrow(nearest_neighbors == 0))){
  #  randEdge = sample(V(g), size = 1, replace = T)
  # add.edges(g, edges = c(infected[[1]], randEdge))
  #}
  keep = unlist(lapply(nearest_neighbors[,2], toss, prob_neighbor=probabilities))
  new_infected = as.numeric(as.character(nearest_neighbors[,1][keep >= 1]))
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
       ylab = "CDF", xlab = "n",
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
    val <- sample(possibilities, size = 1, prob = c(highestProb, (1-highestProb)))
  } 
  else if (x == 2){
    val <- sample(possibilities, size = 1, prob = c(highestProb/2, (1-highestProb)))
  }
  else if (x == 3){
    val <- sample(possibilities, size = 1, prob = c(highestProb/3, (1-highestProb)))
  }
  else if (x == 4){
    val <- sample(possibilities, size = 1, prob = c(highestProb/4, (1-highestProb)))
  }
  else if (x == 5){
    val <- sample(possibilities, size = 1, prob = c(highestProb/5, (1-highestProb)))
  }
  else if (x == 6){
    val <- sample(possibilities, size = 1, prob = c(highestProb/6, (1-highestProb)))
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
     neighbor_nodes <- sample(neighbor_nodes, edgesToDelete, replace = T) # randomly pick the difference of edges to remove
     
     for (j in 1:length(neighbor_nodes)){   # starting loop for the adjacent edges I'll be removing. 
       if (isTRUE(edgeFinder(g, x[i], neighbor_nodes[j]))){ # conditional 
         g = delete.edges(g, (E(g, P=c(x[i], neighbor_nodes[j]))))
       }
     } 
   } 
  }
   return(g)
}
