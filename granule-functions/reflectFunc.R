reflectFunc=function(m, val) {
  m[m > val] <- val
  m[m < 1] <- 1
  return(m)
}