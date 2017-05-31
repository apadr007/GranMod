cleaner.valency4 <- function(nearest.mRNA_list, mRNP.rate){
  output =  vector("list", length(nearest.mRNA_list) ) 
  valency=1:4
  if( isEmpty(nearest.mRNA_list) ){ return(0) }
  
  for (i in 1:length(nearest.mRNA_list)){
    
    filter <- sample(1:0, 1, prob = c(1-mRNP.rate, mRNP.rate) )
    if (filter == 1) {
      valency.node = sample(valency, size = 1, 
                            prob = c(mRNP.rate/4, 
                                     mRNP.rate/3, mRNP.rate/2, 
                                     mRNP.rate/1) )
      
      new.edge = nearest.mRNA_list[[i]]
      
      #valTest = valency_func(g, node_number)
      #valTestdel = which(valTest[,2] > 6)   
      #del = NA
      
      #if (isEmpty(valTestdel)==FALSE) {
      #  del = which(new.edge%in%valTestdel)
      #  if(isEmpty(del)== FALSE) { 
      #    new.edge = new.edge[-del]
      #  }
      #}
      if (isEmpty(new.edge) == TRUE) {
        output[[i]] = new.edge
      } else if (valency.node > length(new.edge)){
        output[[i]] = new.edge
      } else if (length(new.edge) == 1) {
        output[[i]] = new.edge
      } else if (valency.node <= length(new.edge)){
        new.edgeToForm <- new.edge[1:valency.node]
        output[[i]] = new.edgeToForm }
    } else { output[[i]] <- integer(0) }
  }
  return(output)
}