node_number = 50
spaceMax = 25

g=graph.empty(node_number,directed = FALSE)
valency = c(1,2,3,4,5,6)
total_time = 1
V(g)$color = 'lightblue'
P.int.on = 0.6
P.int.off = 0.4

layout.old = layoutGen(node_number, spaceMax)
names = 1:node_number
layout.old = layout.old[[1]][,1:2]

plot(g, layout = layout.old,vertex.size=0.1,
     rescale = FALSE, axes=TRUE, xlim=c(1,spaceMax+1), ylim=c(1,spaceMax+1), asp = 0 )

mRNP_pop = list()
mRNA_pop = list()
gran_pop = list()

strt<-Sys.time()

total_time = 1
while(total_time <= 10){ 
  mRNA = 1:node_number 
  mRNA = which(!mRNA%in%which(degree(g) > 6)) 
  
  layout.old = nodeMover7e(layout.old, g, node_number, spaceMax)
  #if ( any(duplicated(layout.old)) ) { break }
  nearest.mRNA2 = findNearestMovers4(layout.old) # this finds the nearest nodes to every other node
  
  nearest.mRNA3 = cleaner.valency2(nearest.mRNA2, P.int.on)
  g <- edgeRemoverSpace(g, layout.old ) #this deletes an edge if its a neighbor BUT not found in findNearestMovers2
  
  for (m in 1:length(mRNA)) {    
    if ( isTRUE(nearest.mRNA3[[m]] != 0) || !isEmpty(nearest.mRNA3[[m]])) {
      for (j in 1:length(nearest.mRNA3[[m]])) {                                          ### this adds edges
        
        if (isEmpty(nearest.mRNA3[[m]]) == FALSE && isTRUE(degree(g)[nearest.mRNA3[[m]][j]] < 6) ) {  
          g = add.edges(g, c(mRNA[m], nearest.mRNA3[[m]][j])) } 
      }
    } else { next }
  }
  g = simplify(g, remove.loops = TRUE)
  E(g)$color = 'black'
  
  del.index = which(degree(g) > 8) 
  if (isEmpty(del.index) == FALSE){
    for(k in 1:length(del.index)){
      g = del_edge(g, valToDel = del.index[k])
    }
  }
  #  print(which(E(g)$color=='black'))
  
  ##this is important coloring code!!!! 150915 ####
  ##                                           ####
  ##                                           ####
  colorme = which(degree(g) >= 1)
  staysamecolor = which(degree(g) < 1)
  if (isEmpty(colorme) == FALSE){
    for (i in 1:length(colorme)){
      x = colorme[i]
      V(g)$color[x] = 'orange'
    }
  }
  
  
  if (isEmpty(staysamecolor) == FALSE){
    for (i in 1:length(staysamecolor)){
      x = staysamecolor[i]
      V(g)$color[x] = 'lightblue'
    }
  }
  #plot(g, layout = layout.old, vertex.label='',
  #  rescale = F, axes=TRUE, xlim=c(0,spaceMax), ylim=c(0,spaceMax), asp = 0 )
  
  #plot(g, layout = layout.old,
   #    rescale = FALSE, axes=TRUE, xlim=c(1,spaceMax), ylim=c(1,spaceMax), asp = 0 )
  plot(g, layout = layout.old,vertex.size=0.1,
       rescale = FALSE, axes=TRUE, xlim=c(1,spaceMax), ylim=c(1,spaceMax), asp = 0 )
  
  gran_pop[[total_time]] = communityfast.algo(g)
  mRNP_pop[[total_time]] = which(degree(g) >= 1)
  mRNA_pop[[total_time]] = which(degree(g) < 1)
  total_time = total_time + 1  
  print(total_time)
}


#########
print(Sys.time()-strt)


colorme = which(degree(g) >= 1)
for (i in 1:length(colorme)){
  x = colorme[i]
  V(g)$color[x] = 'orange'
}

set.seed(2014)
plot(g, layout = layout.old,vertex.label='',
     rescale = F, axes=TRUE, xlim=c(0,spaceMax), ylim=c(0,spaceMax), asp = 0 )
mRNA.return = length(which(V(g)$color == 'lightblue'))
g = delete.vertices(g, which(V(g)$color == 'lightblue'))
fg= fastgreedy.community(g, modularity = TRUE, membership = TRUE)
V(g)$size = 1
E(g)$count = 1
comm.graph <- contract.vertices(g, fg$membership, vertex.attr.comb=list(size="sum", "ignore"))
comm.graph <- simplify(comm.graph, remove.loops=TRUE, edge.attr.comb=list(count="sum", "ignore"))
which(vertex.attributes(comm.graph)$size>=1)

for(i in 1:length(which(vertex.attributes(comm.graph)$size >= 1)) ){
  V(comm.graph)$color[i] = 'orange' 
}
E(comm.graph)$color = 'white'
comm.graph = add.vertices(comm.graph, mRNA.return)
V(comm.graph)$color[which(V(comm.graph)%in%which(V(comm.graph)$color == 'orange') == FALSE)] = 'lightblue'
V(comm.graph)$size[which(V(comm.graph)%in%which(V(comm.graph)$size >= 1) == FALSE)] = 1
plot(comm.graph, vertex.label='', layout=layout.old)