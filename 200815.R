require('RColorBrewer')
node_number = 50
spaceMax = 15
g=graph.empty(node_number,directed = FALSE)
valency = c(1,2,3,4,5,6)
total_time = 1
V(g)$color = 'lightblue'
P.int.on = 0.8
P.int.off = 0.2
layout.old = layoutGen(node_number, spaceMax)
layout.old = layout.old[[1]]
names = 1:node_number
layout.old = layout.old[,1:2]

plot(g, layout = layout.old, vertex.label='',
     rescale = F, axes=TRUE, xlim=c(0,spaceMax), ylim=c(0,spaceMax), asp = 0 )

mRNP_pop = list()
mRNA_pop = list()

layout.array = array(0, dim = c(node_number,2,100))
#prefilterNearest = list(); clearnerStep = list()

#FunTest = function(g, node_number, valency, total_time, P.on, P.off, layout.old){
#  P.int.on = P.on
#  P.int.off = P.off
#  transmission_rate = P.int.on*(1 - P.int.off)
#  decay_rate = P.int.off*(1 - P.int.on)
#  coins = c(1,0)
#  coins2 = c(P.int.on, (1-P.int.on)) 
#  probabilities = c(transmission_rate, 1 - transmission_rate)                   # for NTI suggestion 
#  probabilities.decay = c(decay_rate, 1 - decay_rate)

strt<-Sys.time()
total_time = 1
while(total_time <= 20){ 
  mRNA = 1:node_number 
  mRNA = which(!mRNA%in%which(degree(g) > 6)) 
  
  coinFlips = sample(1:0, size = 1, replace = TRUE, prob = probabilities)
  # if (isTRUE((coinFlips == 1 ))){                                             
  #  g = add.vertices(g, 1) 
  #  node_number = node_number + 1
  #  mRNA = 1:node_number
  #  mRNA = which(!mRNA%in%which(degree(g) > 6))                             ### reverse select mRNAs with low val
  
  # }                                                                        ### increase population of mRNAs
  #      if (isTRUE(coinFlips == 1) ){                                              ### this increases the population of mRNPs
  #        current_mRNAs = which(V(g)$color=='lightblue')
  #        new.mRNP = sample(current_mRNAs, size = 1)
  #        V(g)$color[new.mRNP] = 'orange'
  #      }
  
  #      coinFlips3 = sample(1:0, size = 1, replace = T,                             
  #                          prob = probabilities.decay)
  
  #if (isTRUE((coinFlips3 == 1))){                                          ### this deletes mRNAs
  #  g = delete.vertices(g, sample(mRNA, size = 1))
  #  node_number = node_number - 1
  #  mRNA = 1:node_number
  #  mRNA = which(!mRNA%in%which(degree(g) > 6))
  #}
  
  #      if (isTRUE((coinFlips3 == 1 & isEmpty(which(V(g)$color=='orange')) == FALSE))){                                            ### this deletes mRNPs
  ##        mRNP_del = which(V(g)$color=='orange')
  #        mRNP_del2 = sample(mRNP_del, size=1)
  ##        if (degree(g)[mRNP_del2] > 0){
  ##          mRNP_connections = unlist(neighborhood(g, order = 1, nodes = mRNP_del2))
  #          mRNP_connections[1] <- NA; mRNP_connections = mRNP_connections[!is.na(mRNP_connections)]
  #          for (y in 1:length(mRNP_connections)){
  #            g = delete.edges(g, E(g, P=c(mRNP_del2, mRNP_connections[y])) )
  #            V(g)$color[mRNP_del2] = 'lightblue'
  #          }
  #        } else { V(g)$color[mRNP_del2] = 'lightblue' }
  #        V(g)$color[which(degree(g) == 0)] = 'lightblue'
  #      }
  layout.old = nodeMover3(layout.old, g, node_number, spaceMax)
  #layout.old = nodeMover(layout.old, g, node_number) # this performs a 2D random walk on each node
  layout.array[,,total_time] = layout.old      
  
  nearest.mRNA2 = findNearestMovers2(layout.old) # this finds the nearest nodes to every other node
  #prefilterNearest[[total_time]] = nearest.mRNA2
  nearest.mRNA3 = cleaner.valency2(nearest.mRNA2, P.int.on)
  g <- edgeRemoverSpace(g, layout.old) #this deletes an edge if its a neighbor BUT not found in findNearestMovers2
  
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
  
  del.index = which(degree(g) > 6) 
  if (isEmpty(del.index) == FALSE){
    for(k in 1:length(del.index)){
      g = del_edge(g, valToDel = del.index[k])
    }
  }
  #  print(which(E(g)$color=='black'))

##this is important coloring code!!!! 150915 ####
##                                           ####
##                                           ####
#                colorme = which(degree(g) >= 1)
#               staysamecolor = which(degree(g) < 1)
#                if (isEmpty(colorme) == FALSE){
#                 for (i in 1:length(colorme)){
#                    x = colorme[i]
#                    V(g)$color[x] = 'orange'
#                  }
#                }
#                if (isEmpty(staysamecolor) == FALSE){
#                  for (i in 1:length(staysamecolor)){
#                    x = staysamecolor[i]
#                    V(g)$color[x] = 'lightblue'
#                  }
#                }
 # plot(g, layout = layout.old,vertex.label='',
#     rescale = F, axes=TRUE, xlim=c(0,spaceMax), ylim=c(0,spaceMax), asp = 0 )
  
  mRNP_pop[[total_time]] = which(degree(g) >= 1)
  mRNA_pop[[total_time]] = which(degree(g) < 1)
  total_time = total_time + 1  
  
}
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
fg= fastgreedy.community(g, modularity = T, membership = T)
V(g)$size = 1
E(g)$count = 1
comm.graph <- contract.vertices(g, fg$membership, vertex.attr.comb=list(size="sum", "ignore"))
comm.graph <- simplify(comm.graph, remove.loops=T, edge.attr.comb=list(count="sum", "ignore"))
which(vertex.attributes(comm.graph)$size>=1)

for(i in 1:length(which(vertex.attributes(comm.graph)$size >= 1)) ){
  V(comm.graph)$color[i] = 'orange' 
}
E(comm.graph)$color = 'white'
comm.graph = add.vertices(comm.graph, mRNA.return)
V(comm.graph)$color[which(V(comm.graph)%in%which(V(comm.graph)$color == 'orange') == FALSE)] = 'lightblue'
V(comm.graph)$size[which(V(comm.graph)%in%which(V(comm.graph)$size >= 1) == FALSE)] = 1
plot(comm.graph, vertex.label='', layout=layout.old)




