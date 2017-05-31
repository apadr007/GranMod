source('~/granule-functions/layoutOverlapFinder_v.R')
library(igraph)
library(reshape2)

set.seed(10)
node_number = 20
spaceMax = 10

g=graph.empty(node_number,directed = FALSE)


total_time = 1
V(g)$color = 'lightblue'
#P.int.on = 0.1

#time 
t = 1e-3

#movement
vel = 1
sigma = vel*t

#reaction
k_on = 1e-2
sigma.kon = k_on*t
P.int.on <- sigma.kon
#reaction
P.int.off = 1e-5
P.int.off <- P.int.off * t

layout.old = layoutGen(node_number, spaceMax)
names = 1:node_number
layout.old = layout.old[[1]][,1:2]

plot(g, layout = layout.old,
       rescale = F, axes=TRUE, xlim=c(0,spaceMax), ylim=c(0,spaceMax), asp = 0 )

RunLength = 500

mRNP_pop = vector("list", RunLength ) 
mRNA_pop = vector("list", RunLength ) 
layout_list = vector("list", RunLength ) 
g_list = vector("list", RunLength ) 


total_time = 1

while(total_time <= RunLength){ 
  mRNA = 1:node_number
  
  mRNA = which(!mRNA%in%which(degree(g) > 4))
  
  whatToDo <- sample(1:3, 1, prob = c(sigma, P.int.on, P.int.off) )
  
  #sigma.test <- sample(c(TRUE,FALSE), size = 1, prob = c(sigma, 1-sigma))
  #if (sigma.test == TRUE){
  
  if (whatToDo == 1){
    layout.old = nodeMover9(layout.old, g, node_number, spaceMax)
   }
  
  if (whatToDo == 2){
    # this finds the nearest nodes to every other node
    nearest.mRNA2 = myDistOneFast(layout.old) 
  
    nearest.mRNA2[ sapply(nearest.mRNA2, is.null) ] <- NA
    
    nearest.mRNA2 <- as.list(nearest.mRNA2)
  
    nearest.mRNA3 = cleaner.valency4(nearest.mRNA2, P.int.on)
    #this deletes an edge if its a neighbor BUT not found in findNearestMovers
    g <- edgeRemoverSpace2(g, layout.old ) 
    
    
    nearest.mRNA3 <- melt(nearest.mRNA3, na.rm = TRUE)
    nearest.mRNA3 <- as.matrix(nearest.mRNA3)
    g = add.edges(g, t(nearest.mRNA3) )
  
    g = simplify(g, remove.loops = TRUE)
  }
  
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

  g <- edgeRemoverSpace2(g, layout.old ) 
  
  plot(g, layout = layout.old,
       rescale = F, axes=TRUE, xlim=c(0,spaceMax), ylim=c(0,spaceMax), asp = 0 )
  
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

Num.Community=list()
size.Community.mean=list()
size.Community.sd=list()
for (i in 1:500) {
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
plot(val.gran,lwd=2, pch=19, col='magenta', 
     ylab=c('Number of Communities'), xlab=c('Timestep'), type='l' )
#plot community size over time
plot(val.comm.size, lwd=2, cex=0.3, pch=19, 
     ylab=c('Fraction of Mean Community Size'), xlab=c('Timestep'), type='l')
#plot mRNP over time
plot(val.RNP, lwd=2, pch=19, col='blue', 
     ylab=c('Number of Dimers'), xlab=c('Timestep'), type='l' )
#plot mRNAs over time
plot(val.mRNA,lwd=2, pch=19, col='red', 
     ylab=c('Number of mRNA'), xlab=c('Timestep'), type='l' )
