myDistOneFast <- function(m) {
  v1 <- m[,1L]; v2 <- m[,2L]
  origrs <- rowSums(m)
  mySort <- order(origrs)
  rs <- origrs[mySort]
  myDiff <- c(0L, diff(rs))
  brks <- which(myDiff > 0L)
  lenB <- length(brks)
  n <- nrow(m)
  myL <- vector("list", length = n)
  
  findRows <- function(v, s, r, u1, u2) {
    lapply(v, function(x) {
      sx <- s[x]
      tv1 <- s[r]
      tv2 <- tv1[which(abs(u1[sx] - u1[tv1]) <= 1)]
      tv2[which(abs(u2[sx] - u2[tv2]) <= 1)]
    })
  }
  
  t1 <- brks[1L]; t2 <- brks[2L]
  ## setting first index in myL
  myL[mySort[1L:(t1-1L)]] <- findRows(1L:(t1-1L), mySort, t1:(t2-1L), v1, v2)
  k <- t0 <- 1L
  
  while (k < (lenB-1L)) {
    t1 <- brks[k]; t2 <- brks[k+1L]; t3 <- brks[k+2L]
    vec <- t1:(t2-1L)
    if (myDiff[t1] == 1L) {
      if (myDiff[t2] == 1L) {
        myL[mySort[vec]] <- findRows(vec, mySort, c(t0:(t1-1L), t2:(t3-1L)), v1, v2)
      } else {
        myL[mySort[vec]] <- findRows(vec, mySort, t0:(t1-1L), v1, v2)
      }
    } else if (myDiff[t2] == 1L) {
      myL[mySort[vec]] <- findRows(vec, mySort, t2:(t3-1L), v1, v2)
    }
    if (myDiff[t2] > 1L) {
      if (myDiff[t3] > 1L) {
        k <- k+2L; t0 <- t2
      } else {
        k <- k+1L; t0 <- t1
      }
    } else {k <- k+1L; t0 <- t1}
  }
  
  ## setting second to last index in myL
  if (k == lenB-1L) {
    t1 <- brks[k]; t2 <- brks[k+1L]; t3 <- n+1L; vec <- t1:(t2-1L)
    if (myDiff[t1] == 1L) {
      if (myDiff[t2] == 1L) {
        myL[mySort[vec]] <- findRows(vec, mySort, c(t0:(t1-1L), t2:(t3-1L)), v1, v2)
      } else {
        myL[mySort[vec]] <- findRows(vec, mySort, t0:(t1-1L), v1, v2)
      }
    } else if (myDiff[t2] == 1L) {
      myL[mySort[vec]] <- findRows(vec, mySort, t2:(t3-1L), v1, v2)
    }
    k <- k+1L; t0 <- t1
  }
  
  t1 <- brks[k]; vec <- t1:n
  if (myDiff[t1] == 1L) {
    myL[mySort[vec]] <- findRows(vec, mySort, t0:(t1-1L), v1, v2)
  }
  
  myL
}