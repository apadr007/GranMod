source('~/granule-functions/layoutOverlapFinder_v.R')
source('~/granule-functions/matrixIndex.R')
library(igraph)
library(reshape2)

load.functions <- list.files('~/granule-functions', full.names = TRUE)
for (i in 1:length(load.functions)){
  source(load.functions[i])
}


set.seed(10)
start.time <- Sys.time()
node_number = 80
spaceMax = 10


g=graph.empty(node_number,directed = FALSE)


total_time = 1
V(g)$color = 'lightblue'
#P.int.on = 0.1

#time 
t = 1e-3

#movement
vel = 5e-1
sigma = vel*t

#reaction
k_on = 1e-2
sigma.kon = k_on*t
P.int.on <- sigma.kon
#reaction
P.int.off = 1e-2
P.int.off <- P.int.off * t

layout.old = layoutGen(node_number, spaceMax)
names = 1:node_number
layout.old = layout.old[[1]][,1:2]

plot(g, layout = layout.old,
       rescale = F, axes=TRUE, xlim=c(0,spaceMax), ylim=c(0,spaceMax), asp = 0 )

# This makes 500 time steps
#Run <- rep(1:20, 25)
#get total length 
#RunLength = length(Run)
#nodez <- Run
RunLength = 200

gran_pop = vector('list', RunLength)
mRNP_pop = vector("list", RunLength ) 
mRNA_pop = vector("list", RunLength ) 
layout_list = vector("list", RunLength ) 
g_list = vector("list", RunLength ) 


total_time = 1

while(total_time <= RunLength){ 
  
  whatToDo <- sample(1:3, 1, prob = c(sigma, P.int.on, P.int.off) )
  
  if (whatToDo == 1){
    layout.old = nodeMover10(layout.old, g, node_number, spaceMax)
   }
  
  if (whatToDo == 2){
    
    #this randomly deletes edge 
    DegreeType <- sample(0:4, size = 1, 
                         prob = c(P.int.on/5,P.int.on/4, P.int.on/3, 
                                  P.int.on/2, P.int.on/1) )
    
    toBind.index <-  which(degree(g) == DegreeType)
    if ( !isEmpty(toBind.index) && DegreeType >= 0 ) {
      toBind <- sample( toBind.index, 1)
      mRNA <- sample(1:node_number, 1)
    
      x <- indexLookUp(mRNA, c('A', 'B', 'C', 'D'), layout.old)
      optionsToBind <- x[layoutOverlapFinder.m(layout.old, x) == 1,]
      optionsToBind <- `if`(isEmpty(optionsToBind), layout.old[mRNA, ], optionsToBind)
    
    if ( is.null( nrow(optionsToBind) ) ) {
      indexFound <- indexFinder(layout.old, optionsToBind)
      g = add.edges(g, c(mRNA, indexFound) )
      g = simplify(g, remove.loops = TRUE) } else { 
        # this matches position to index using Rcpp. Function needs to be compiled every time
        NodeOptions <- matrixIndex(optionsToBind, layout.old)
        NodeOptions <- NodeOptions[NodeOptions != 0]
        
        #this finds the degrees of the available mRNPs to bind to
        degreeTypesAvailable <- degree(g)[NodeOptions]
        
        # this preferentially selects higher degrees than lower degree nearby nodes to form
        degreeTypesAvailable.selected <- sample(degreeTypesAvailable, 1, 
               prob = c( mRNP_selector_prob(degreeTypesAvailable, P.int.on) ) )
        degreeTypesAvailable.indx <- sample(which(degreeTypesAvailable == degreeTypesAvailable.selected, arr.ind = TRUE), 1)
        
        NewEdgeNode <- NodeOptions[degreeTypesAvailable.indx]
        g = add.edges(g, c(mRNA, NewEdgeNode) )
        g = simplify(g, remove.loops = TRUE)
      }
    }
  }
    
     #this deletes an edge if its a neighbor BUT not found in findNearestMovers
    #g <- edgeRemoverSpace(g, layout.old )
    
  
  if (whatToDo == 3){
    #this randomly deletes edge 
    DegreeType <- sample(0:4, size = 1, 
                         prob = c(P.int.off/1,P.int.off/2, P.int.off/3, 
                                  P.int.off/4, P.int.off/5) )
    toDelete.index <- which(degree(g) == DegreeType)
    if ( !isEmpty(toDelete.index) && DegreeType > 0 ) {
      toDelete <- sample( toDelete.index, 1)
      g <- delete_edges(g, E(g)[ from(toDelete) ] )
    }
  }
  
  E(g)$color = 'black'
  
  V(g)$color[degree(g) >= 1] <- 'orange'
  V(g)$color[degree(g) < 1] <- 'lightblue'

  g <- edgeRemoverSpace(g, layout.old ) 
  
  #plot(g, layout = layout.old,
  #     rescale = F, axes=TRUE, xlim=c(0,spaceMax), ylim=c(0,spaceMax), asp = 0 )
  
  gran_pop[[total_time]] = walktrap.community(g)
  mRNP_pop[[total_time]] = which(degree(g) >= 1)
  names(mRNP_pop) <- 'mRNP'
  mRNA_pop[[total_time]] = which(degree(g) < 1)
  names(mRNA_pop) <- 'mRNA'
  layout_list[[total_time]] = layout.old
  g_list[[total_time]] <- g
  
  print(total_time)
  
  total_time = total_time + 1  

  } 

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

Num.Community=list()
size.Community.mean=list()
size.Community.sd=list()
for (i in 1:RunLength) {
  names(gran_pop[[i]]$membership) <- 1:node_number
  gran.length <- length(table(gran_pop[[i]]$membership)[table(gran_pop[[i]]$membership) ])
  gran.members <- lapply(1:gran.length, function(x) as.numeric(names(gran_pop[[i]]$membership[gran_pop[[i]]$membership == x])))
  Num.Community[[i]] <- length(unlist(lapply(gran.members, length))[unlist(lapply(gran.members, length)) ])
  size.Community.mean[[i]] <- mean(unlist(lapply(gran.members,length)))
  size.Community.sd[[i]] <- sd(unlist(lapply(gran.members,length)))
  }

val.gran <- c(node_number, unlist(Num.Community))
val.RNP <- (c(0,unlist(lapply(mRNP_pop, length))))
val.mRNA <- c(node_number,unlist(lapply(mRNA_pop, length)))
val.comm.size <- c(unlist(size.Community.mean))

par(pin=c(2,2), tcl=0.25, ps=9, family="Helvetica")
#plot number of communities over time
pdf('~/granule-model/plot1_test.pdf')
plot(val.gran,lwd=2, pch=19, col='magenta', 
     ylab=c('Number of Communities'), xlab=c('Timestep'), type='l' )
dev.off()
#plot community size over time
par(pin=c(2,2), tcl=0.25, ps=9, family="Helvetica")
pdf('~/granule-model/plot2_test.pdf')
plot(val.comm.size, lwd=2, cex=0.3, pch=19, 
     ylab=c('Mean Community Size'), xlab=c('Timestep'), type='l')
dev.off()
#plot mRNP over time
par(pin=c(2,2), tcl=0.25, ps=9, family="Helvetica")
pdf('~/granule-model/plot3_test.pdf')
plot(val.RNP, lwd=2, pch=19, col='blue', 
     ylab=c('Number of Dimers'), xlab=c('Timestep'), type='l', ylim=c(0, node_number) )
dev.off()
#plot mRNAs over time
par(pin=c(2,2), tcl=0.25, ps=9, family="Helvetica")
pdf('~/granule-model/plot4_test.pdf')
plot(val.mRNA,lwd=2, pch=19, col='red', 
     ylab=c('Number of mRNA'), xlab=c('Timestep'), type='l', ylim=c(0, node_number) )
dev.off()