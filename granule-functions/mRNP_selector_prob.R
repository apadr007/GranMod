mRNP_selector_prob <- function(x, p){
  x[x == 0 ] <- p/4
  x[x == 1 ] <- p/3
  x[x == 2 ] <- p/2
  x[x == 3 ] <- p
  x[x == 4 ] <- 0
  return(x)
}