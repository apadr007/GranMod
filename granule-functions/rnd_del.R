rnd.del = function(g, valency, mRNP.rate){
  #weighted selection of the valency
  valency.node = sample(valency, size = 1, replace = TRUE, 
                        prob = c(mRNP.rate/1, 
                                 mRNP.rate/2, mRNP.rate/3, 
                                 mRNP.rate/4) )
  #get list of nodes with a randomly selected valency
  sub = which(valency.node==degree(g) )
  
  if ( isEmpty(sub) ) {
    g } else {
      #this gets the index of the edge for the node selected by searching for the node# rather than index. Necessary! 
      sub2 <- sample(sub, 1)
      rnd_del <- sample( which( get.edgelist(g)==sub2 , 1) )
      
      #rnd_del <- sample( sub, size = 1)
      delete_edges(g, rnd_del) 
      }
}
