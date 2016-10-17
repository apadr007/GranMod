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
