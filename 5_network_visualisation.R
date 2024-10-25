######################################
# Network Visualisation using iGraph #
######################################

source('1_libraries.R')
source('2_differential_gene_expression_analysis.R')
source('4_coexpression_network_analysis.R')

thr = 0.65
adjacencyMat <- adjacency(datExpr, power = 1, type="signed hybrid");
m = lower.tri(adjacencyMat,diag = FALSE)
m = melt(m)
net = melt(adjacencyMat)
net = net[which(m$value==TRUE),]
names(net) = c("Protein1","Protein2","Correlation")
net = net[which(net$Correlation>=thr),]

##################
### Functions ###
#################

#####################
# Continuous colors #
#####################

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
  pal <- colorRampPalette(c(negcol, midcol, poscol))(101)
  xx <- round(range01(x), 2) * 100 + 1 # added 1 as indexing from 0 does not work in R
  return(pal[xx])
}

#########
# main #
########
adjlist = data.frame(net, row.names = NULL)

lfcTab = data.frame(node = MWtest_results$UniProt, lfc = MWtest_results$estimate)


moduleTab = data.frame(node = geneInfo$Proteins, 
                       mod = geneInfo$moduleColor) %>%
  column_to_rownames('node')

g <- graph_from_data_frame(adjlist, directed = TRUE)

moduleTab = moduleTab[V(g)$name,"mod", drop=FALSE]
moduleTab$mod <- as.factor(moduleTab$mod)
levels(moduleTab$mod)

markgroups = moduleTab %>% rownames_to_column("node") %>%
  group_by(mod) %>% 
  group_map(function(x,i){return(x$node)})
names(markgroups) = moduleTab %>% rownames_to_column("node") %>%
  group_by(mod) %>% 
  group_map(function(x,i){ return(as.character(i %>% pull))})

sublfcTab = lfcTab %>% column_to_rownames("node")
V(g)$color <- getColor(sublfcTab[ V(g)$name,"lfc"])

initcoords = moduleTab[V(g)$name,"mod"] %>% 
  sapply(function(x){
    switch(as.character(x), 
           red = c(1.,0.), 
           blue = c(0.,0.), 
           turquoise = c(0.,1.), 
           brown = c(1.,1.))}) %>% as.data.frame()


for(x in colnames(initcoords)){
  initcoords[1,x] = initcoords[1,x] + rnorm(1) * 0.01
  initcoords[2,x] = initcoords[2,x] + rnorm(1) * 0.01
}


coords = layout_with_fr(g, coords = t(initcoords))

e = as_edgelist(g, names = F)
coords = qgraph.layout.fruchtermanreingold(e, vcount = vcount(g), groups = markgroups,
                                           init = t(initcoords), repulse.rad = 10000)

plot(g,
     mark.groups = markgroups,
     mark.col = c(),
     mark.border = names(markgroups),
     layout = coords,
     vertex.size = 10,
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     edge.width = 1,
     edge.arrow.size = 0.2)