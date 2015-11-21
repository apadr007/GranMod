require('RColorBrewer')
node_number = 45
spaceMax = 10
g=graph.empty(node_number,directed = FALSE)
valency = c(1,2,3,4,5,6)
total_time = 1
V(g)$color = 'lightblue'
P.int.on = 0.3
P.int.off = 0.7
layout.old = layoutGen(node_number, spaceMax)
layout.old = layout.old[[1]]
names = 1:node_number
layout.old = layout.old[,1:2]
deg = list()
granules = list()

plot(g, layout = layout.old, vertex.label='',
     rescale = F, axes=TRUE, xlim=c(0,spaceMax), ylim=c(0,spaceMax), asp = 0 )

mRNPc_pop = list()
mRNP_pop = list()
mRNA_pop = list()

#layout.array = array(0, dim = c(node_number,2,100))
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

total_time2 = 0.001
total_time = 1
while(total_time2 <= 0.01){ 
  mRNA = 1:node_number 
  mRNA = which(!mRNA%in%which(degree(g) > 6)) 
  
  coinFlips = sample(1:0, size = 1, replace = TRUE, prob = probabilities)
 
  layout.old = nodeMover3(layout.old, g, node_number, spaceMax)
  nearest.mRNA2 = findNearestMovers2(layout.old) # this finds the nearest nodes to every other node
  
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
  # plot(g, layout = layout.old,vertex.label='',
  #     rescale = F, axes=TRUE, xlim=c(0,spaceMax), ylim=c(0,spaceMax), asp = 0 )
  
  mRNPc_pop[[total_time]] = which(degree(g) >= 3)
  mRNP_pop[[total_time]] = length(which(V(g)$color == 'orange' & degree(g) < 3))
  mRNA_pop[[total_time]] = which(degree(g) < 1)
  
  deg[[total_time]] = which( degree(g) >= 2 )
  
  comm.rnp = communityfast.algo(g)
  granules[[total_time]] = comm.rnp[which(comm.rnp >= 3)]
  total_time = total_time + 1  
  total_time2 = total_time2 + 0.001

}

print(Sys.time()-strt)

mRNPc = c(0,unlist(lapply(mRNPc_pop, length)))
mRNP = c(0,unlist(mRNP_pop))
mRNA = c(node_number,unlist(lapply(mRNA_pop, length)))
grans = unlist(lapply(granules, length))

#colorme = which(degree(g) >= 1)
#for (i in 1:length(colorme)){
#  x = colorme[i]
#  V(g)$color[x] = 'orange'
#}

#plot(g, layout = layout.old,vertex.label='',
#     rescale = F, axes=TRUE, xlim=c(0,spaceMax), ylim=c(0,spaceMax), asp = 0 )

#lec <- leading.eigenvector.community(g)
#lec

#lec2=leading.eigenvector.community(g, start=membership(lec))

