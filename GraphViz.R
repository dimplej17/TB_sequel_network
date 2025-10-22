library(igraph)
library(tidyverse)


################## functions
# continuous colors
range01 <- function(x) {
    isPos <- x > 0
    pos <- c(0, x[isPos])
    neg <- c(0, x[!isPos])
    pos <- (pos - min(pos)) / (max(pos) - min(pos))
    neg <- (neg - min(neg)) / (max(neg) - min(neg))
    xx <- rep(-1, length(x))
    xx[isPos] <- pos[2:length(pos)] / 2 + 0.5
    xx[!isPos] <- neg[2:length(neg)] / 2
    return(xx)
}

getColor <- function(x, poscol = "red", midcol = "white", negcol = "blue") {
    pal <- colorRampPalette(c(poscol, midcol, negcol))(101)
    xx <- round(range01(x), 2) * 100 + 1 # added 1 as indexing from 0 does not work in R
    return(pal[xx])
}


########## main

adjlist = data.frame(node1 = sample(LETTERS, 100, replace = TRUE), 
        node2 = sample(LETTERS, 100, replace = TRUE), 
        weight = runif(100,-1,1))

lfcTab = data.frame(node = LETTERS, 
        lfc = runif(length(LETTERS),-5,5)) %>% 
        column_to_rownames('node')

moduleTab = data.frame(node = LETTERS, 
        mod = sample(c('turquoise','grey','red','darkgreen'), 
                length(LETTERS), replace = T) %>% as.factor()
        ) %>% 
        column_to_rownames('node')

g <- graph_from_data_frame(adjlist, directed = TRUE)
g = set_vertex_attr(g, "color", value = lfcTab[V(g),])

markgroups = moduleTab %>% rownames_to_column('node') %>%
        group_by(mod) %>% 
        group_map(function(x,i){return(x$node)})
names(markgroups) = moduleTab %>% rownames_to_column('node') %>%
        group_by(mod) %>% 
        group_map(function(x,i){ return(i %>% pull)}) %>% unlist()

V(g)$color <- getColor(V(g)$color)


plot(g, mark.groups = markgroups, mark.col = c(), mark.border = names(markgroups))

# using additional library
library(qgraph)

coords = moduleTab[V(g)$name,'mod'] %>% sapply( function(x){
        switch(as.character(x),
        turquoise = c(0,0),
        red = c(1,0),
        grey =c(1,1),
        darkgreen = c(0,1))
}) %>% t() %>% as.data.frame()

for(row in rownames(coords)){
        coords[row, 1] = coords[row, 1] + runif(1,-0.2,0.2)
        coords[row, 2] = coords[row, 2] + runif(1, -0.2, 0.2)
}

coords = coords %>% as.matrix

e <- as_edgelist(g, names = FALSE)
#  niter = 40, max.delta = 0.01 are the parameters that prevent the nodes from moving too far out
lo = qgraph.layout.fruchtermanreingold(e, vcount = vcount(g), init = coords, niter = 40, max.delta = 0.01)
pdf("/home/obaranov/GraphViz.pdf")
p = plot(g, layout = lo, mark.groups = markgroups, mark.col = c(), mark.border = names(markgroups))
dev.off()