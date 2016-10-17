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
