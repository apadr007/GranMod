indexLookUp4 = function(mover.m, posit.move, ...){
  out = list()
  x = hx[[ all.Letters[mover.m] ]]
  y = hy[[ all.Letters[mover.m] ]]
  
  A = c(x + 1, y)
  B = c(x - 1, y)
  C = c(x, y + 1)
  D = c(x, y - 1)
  E = c(x, y)
  m = as.matrix(rbind(A, B, C, D, E))
  
  out <- m[ posit.move, ]
  return(out)
}