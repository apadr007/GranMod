layoutGetter = function(layout.old, nodes.you.want){
    del = nrow(layout.old) - nodes.you.want
i = 1
while (i <= del){
    x = sample(1:nrow(layout.old), size = 1)
    layout.old = layout.old[-x,]
    i = i + 1
}
return(layout.old)}
