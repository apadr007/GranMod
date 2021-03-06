myDistOne = function(m) {
  v1 <- m[,1L]; v2 <- m[,2L]
  rs <- rowSums(m)
  lapply(seq_along(rs), function(x) {
    t1 <- which(abs(rs[x] - rs) == 1)
    t2 <- t1[which(abs(v1[x] - v1[t1]) <= 1)]
    t2[which(abs(v2[x] - v2[t2]) <= 1)]
  })
}

myfun = function(x){
  doo <- function(j) {
    j.mat <- matrix(rep(j, length = length(x)), ncol = ncol(x), byrow = TRUE)
    j.abs <- abs(j.mat - x)
    return(which(rowSums(j.abs) == 1))
  }
  return(apply(x, 1, doo))
}