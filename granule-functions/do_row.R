do_row <- function(r) {
  r_mat <- matrix(rep(r, length = length(x)), ncol = ncol(x), byrow = TRUE)
  abs_dist <- abs(r_mat - x)
  return(which(rowSums(abs_dist) == 1))
}