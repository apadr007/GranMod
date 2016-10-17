findNearestMovers3 = function(t){
x=matrix(); x2=list()
y = matrix(); y2 = list()
for (i in 1:nrow(t)){
    for (j in 1:nrow(t)){
        if (abs(t[i,1] - t[j,1]) <= 0.4) {
            x[j] = j
        } else { x[j] = NA }
        if (abs(t[i,2] - t[j,2]) <= 0.4) {
            y[j] = j
        } else { y[j] = NA }
    }
    x2[[i]] = x
    y2[[i]] = y
}
for (i in 1:length(x2)){
    x2[[i]] = x2[[i]][!x2[[i]] == i]
    y2[[i]] = y2[[i]][!y2[[i]] == i]
}
x2 = lapply(x2, function(x) x[!is.na(x)])
y2 = lapply(y2, function(x) x[!is.na(x)])

#this intersects the x and y axis to find factors that are close in both axes
z = list()
for (i in 1:length(x2)){
    z[[i]] = intersect(unlist(x2[[i]]), unlist(y2[[i]]))
}
return(z)}
