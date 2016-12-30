source('~/granule-functions/layoutOverlapFinder_v.R')
node_number = 30
spaceMax = 10

g=graph.empty(node_number,directed = FALSE)
valency = c(1,2,3,4,5,6)
total_time = 1
V(g)$color = 'lightblue'
P.int.on = .1
P.int.off = 0.1

#movement
vx = 1
t = 1e-3
sigma.x = vx*t

#reaction
k_on = 1e-04
t = 1e-3
sigma.kon = k_on*t


layout.old = layoutGen(node_number, spaceMax)
names = 1:node_number
layout.old = layout.old[[1]][,1:2]

plot(g, layout = layout.old,
rescale = F, axes=TRUE, xlim=c(0,spaceMax), ylim=c(0,spaceMax), asp = 0 )

mRNP_pop = list()
mRNA_pop = list()
#layout_list = list()[1:100]

total_time = 1

#Rprof('myFunction.out')
while(total_time < 50){
    print(total_time)
    mRNA = 1:node_number
    
    
    mRNA = which(!mRNA%in%which(degree(g) > 4))
    valency = c(1,2,3,4,5,6)
    
    #sigma.test <- sample(c(TRUE,FALSE), size = 1, prob = c(sigma, 1-sigma))
    #if (sigma.test == TRUE){
    layout.old = nodeMover7d(layout.old, g, node_number, spaceMax)
    # }
    # this finds the nearest nodes to every other node
    nearest.mRNA2 = myDistOne(layout.old)
    if ( !is.list(nearest.mRNA2) ) nearest.mRNA2 = list(nearest.mRNA2)
    nearest.mRNA3 = cleaner.valency2(nearest.mRNA2, P.int.on)
    #this deletes an edge if its a neighbor BUT not found in findNearestMovers
    g <- edgeRemoverSpace(g, layout.old)
    
    
    #####
    ##### placing these for loops inside a function called node_adder(mRNA, nearest.mRNA3)
    for (m in 1:length(nearest.mRNA3)) {
        if ( isTRUE(nearest.mRNA3[[m]] != 0) || !isEmpty(nearest.mRNA3[[m]])) {
            ### this adds edges
            for (j in 1:length(nearest.mRNA3[[m]])) {
                
                if (isEmpty(nearest.mRNA3[[m]]) == FALSE && isTRUE(degree(g)[nearest.mRNA3[[m]][j]] < 6) ) {
                    g = add.edges(g, c(mRNA[m], nearest.mRNA3[[m]][j])) }
            }
        } else { next }
    }
    #####
    #####
    
    #reaction
    #sigma.kon.test <- sample(c(TRUE,FALSE), size = 1, prob = c(sigma.kon, 1-sigma.kon))
    #if (sigma.kon.test == TRUE){
    
    #
    # g <- node_adder(mRNA, nearest.mRNA3)
    #  }
    
    g = simplify(g, remove.loops = TRUE)
    E(g)$color = 'black'
    
    #del.index = which(degree(g) > 8)
    #if (isEmpty(del.index) == FALSE){
    #  for(k in 1:length(del.index)){
    #    g = del_edge(g, valToDel = del.index[k])
    #  }
    #}
    
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
    
    plot(g, layout = layout.old,
    rescale = F, axes=TRUE, xlim=c(0,spaceMax), ylim=c(0,spaceMax), asp = 0 )
    
    
    #gran_pop[[total_time]] = communityfast.algo(g)
    mRNP_pop[[total_time]] = which(degree(g) >= 1)
    names(mRNP_pop) <- 'mRNP'
    mRNA_pop[[total_time]] = which(degree(g) < 1)
    names(mRNA_pop) <- 'mRNA'
    
    #layout_list[[total_time]] <- layout.old
    
    total_time = total_time + 1
}
#Rprof(NULL)
#summaryRprof('myFunction.out')


val.RNP <- (c(0,unlist(lapply(mRNP_pop, length))))
val.mRNA <- c(node_number,unlist(lapply(mRNA_pop, length)))

#plot mRNP over time
plot(val.RNP/max(val.RNP), cex=0.5, pch=19, col='blue',
ylab=c('Fraction of Dimers'), xlab=c('Timestep'))
#plot mRNAs over time
plot(val.mRNA/max(val.mRNA),cex=0.5, pch=19, col='red',
ylab=c('Fraction of of mRNA'), xlab=c('Timestep') )

# find densely connected regions
colorme = which(degree(g) >= 1)
for (i in 1:length(colorme)){
    x = colorme[i]
    V(g)$color[x] = 'orange'
}

fg = walktrap.community(g, modularity = TRUE, membership = TRUE)
fg$membership
mRNA.return = length(which(V(g)$color == 'lightblue'))
g = delete.vertices(g, which(V(g)$color == 'lightblue'))
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
plot(comm.graph, vertex.label='')

# get size of granule (how many nodes inside each densely connected region of the graph)
gran.size <- V(comm.graph)$size[V(comm.graph)$size > 1]
# get number of granules
gran.num <- length(gran.size)
