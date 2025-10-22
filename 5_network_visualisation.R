######################################
# Network Visualisation using iGraph #
######################################

library(tidyverse)
library(vroom)
library(reshape2)
library(WGCNA)
library(igraph)
library(ggplot2)
library(qgraph)

basedir = '/home/obaranov/projects/TBSequel/DimplesProject/'

# adjacencyMat = readRDS(paste0(basedir, '/OutputData/adjacency.RDS'))
MWtest_results = vroom(paste0(basedir, '/OutputData/pairedM0vsM6_wilcox.csv'))
geneInfo = vroom(paste0(basedir, '/OutputData/geneModuleAssociation.csv') )
datExpr = vroom(paste0(basedir, '/FilteredData/scaledExpressionData_wo3222.csv'))

# thr = 0.2
thr = 0.45
# this step actually calculates a unsigned correlation matrix rather than a network structure.... 
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

# replace grey with something more visible

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

# module colors change :(
initcoords = moduleTab[V(g)$name,"mod"] %>% 
  sapply(function(x){
    switch(as.character(x), 
           grey = c(0.5,0.5), 
           blue = c(0.25,0.), 
           black = c(0.75,0), 
           brown = c(0.,0.5),
           green = c(1.,0.5), 
           pink = c(0.,1.), 
           red = c(0.75,1.)
           )
        }) %>% as.data.frame()



for(x in colnames(initcoords)){
  initcoords[1,x] = initcoords[1,x] + rnorm(1) * 0.01
  initcoords[2,x] = initcoords[2,x] + rnorm(1) * 0.01
}


coords = layout_with_fr(g, coords = t(initcoords))
dim(initcoords)

e = as_edgelist(g, names = F)
coords = qgraph.layout.fruchtermanreingold(e, vcount = vcount(g), groups = markgroups,
                                           init = t(initcoords),
                                           niter = 100, max.delta = 0.01, repulse.rad = 1)

pdf(paste0(basedir, '/OutputData/graph.pdf'))
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
dev.off()

#########################################
### visualisation of the initial data ###
#########################################

sfnetwork = readRDS(paste0(basedir,'/FilteredData/adjacency.RDS'))
mods = readRDS(paste0(basedir,'/FilteredData/modules.RDS'))
# sfg = graph_from_adjacency_matrix(sfnetwork)
sfg = graph_from_adjacency_matrix(sfnetwork > 0)
sfg = simplify(sfg)

for(edge in E(sfg)){
  tail = igraph::tail_of(sfg, edge)$name
  head = igraph::head_of(sfg, edge)$name
  set_edge_attr(sfg, edge,'weight',sfnetwork[head, tail])
}

tmpprot = c()
tmpcol = c()
for(mod in names(mods)){
  tmpprot = c(tmpprot,mods[[mod]])
  tmpcol = c(tmpcol, rep(mod, length(mods[[mod]])))
}
names(tmpcol) = tmpprot

moduleTab = data.frame(node = geneInfo$Proteins, 
                       mod = geneInfo$moduleColor) %>%
  column_to_rownames('node')

moduleTab = moduleTab[V(sfg)$name,"mod", drop=FALSE]

markgroups = moduleTab %>% rownames_to_column("node") %>%
  group_by(mod) %>% 
  group_map(function(x,i){return(x$node)})
names(markgroups) = moduleTab %>% rownames_to_column("node") %>%
  group_by(mod) %>% 
  group_map(function(x,i){ return(as.character(i %>% pull))})


initcoords = moduleTab[V(sfg)$name,"mod"] %>% 
  sapply(function(x){
    switch(as.character(x), 
           grey = c(0.5,0.5) + runif(1, min = -0.1, max = 0.1), 
           blue = c(0.25,0.) + runif(1, min = -0.1, max = 0.1), 
           black = c(0.75,0) + runif(1, min = -0.1, max = 0.1), 
           brown = c(0.,0.5) + runif(1, min = -0.1, max = 0.1),
           green = c(1.,0.5) + runif(1, min = -0.1, max = 0.1), 
           pink = c(0.,1.) + runif(1, min = -0.1, max = 0.1), 
           red = c(0.75,1.) + runif(1, min = -0.1, max = 0.1)
           )
        }) %>% as.data.frame()

e = as_edgelist(sfg, names = F)
coords = qgraph.layout.fruchtermanreingold(e, edge_attr(sfg, 'weight'), vcount = vcount(sfg), groups = markgroups)
                                          #  init = t(initcoords),
                                          #  niter = 100, max.delta = 0.01, repulse.rad = 1)

pdf(paste0(basedir, '/OutputData/fullGraph.pdf'))
# jpeg(paste0(basedir, '/OutputData/fullGraph.jpg'))

plot(sfg,
     mark.groups = markgroups,
     mark.col = markgroups,
     mark.border = names(markgroups),
     layout = coords,
     vertex.size = 10,
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     edge.width = 1,
     edge.arrow.size = 0.2)
dev.off()


################################
### export network for gephi ###
################################

mods=readRDS(paste0(basedir,'/FilteredData/modules.RDS'))
adjlist=readRDS(paste0(basedir,'/FilteredData/adjacency.RDS'))
lfc=readRDS(paste0(basedir,'/FilteredData/MWtest_results_paired.RData')) %>% 
      column_to_rownames('UniProt')

adjmat = apply(adjlist, 2, function(x){x>0})
nw = graph_from_adjacency_matrix(adjmat)
nw = igraph::set_vertex_attr(nw, 'name',V(nw), colnames(adjlist) )

weightlist = c()
for(edge in E(nw)){
  target = tail_of(nw,edge)
  source = head_of(nw,edge)
  weightlist = c(weightlist, adjlist[target, source])
}

nw = set_edge_attr(nw, 'weight',E(nw), weightlist)
tmp = unlist(mods)
modcol = names(tmp) %>% sapply(function(x){gsub('[0-9]','',x)}) %>% unname
nm = unname(tmp)
modtib = tibble(name = nm,module = modcol) %>% column_to_rownames('name')

nw = set_vertex_attr(nw, 'module', V(nw), modtib[V(nw)$name,])
nw = set_vertex_attr(nw, 'symbol', V(nw), lfc[V(nw)$name,'Assay'])
nw = set_vertex_attr(nw, 'lfc', V(nw), lfc[V(nw)$name,'estimate'])
nw = set_vertex_attr(nw, 'adjP', V(nw), lfc[V(nw)$name,'Adjusted_pval'])

write_graph(simplify(nw),paste0(basedir,'/OutputData/graph.graphml'), format = 'graphml')