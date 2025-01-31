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

net_with_symbol <- read.csv("network_correlations_with_symbols_only", stringsAsFactors = TRUE)
geneInfo_with_symbol <- read.csv("geneInfo1.csv", stringsAsFactors = TRUE)

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
adjlist2 = data.frame(net_with_symbol, row.names = NULL)

lfcTab2 = data.frame(node = MWtest_results$Assay, lfc = MWtest_results$estimate)


moduleTab2 = data.frame(node = geneInfo_with_symbol$Symbol, 
                       mod = geneInfo_with_symbol$moduleColor) %>%
  column_to_rownames('node')

g2 <- graph_from_data_frame(adjlist2, directed = TRUE)

moduleTab2 = moduleTab2[V(g2)$name,"mod", drop=FALSE]
moduleTab2$mod <- as.factor(moduleTab2$mod)
levels(moduleTab$mod)


markgroups2 = moduleTab2 %>% rownames_to_column("node") %>%
  group_by(mod) %>% 
  group_map(function(x,i){return(x$node)})
names(markgroups2) = moduleTab2 %>% rownames_to_column("node") %>%
  group_by(mod) %>% 
  group_map(function(x,i){ return(as.character(i %>% pull))})

sublfcTab2 = lfcTab2 %>% column_to_rownames("node")
V(g2)$color <- getColor(sublfcTab2[V(g2)$name,"lfc"])

initcoords2 = moduleTab2[V(g2)$name,"mod"] %>% 
  sapply(function(x){
    switch(as.character(x), 
           red = c(1.,0.), 
           blue = c(0.,0.), 
           turquoise = c(0.,1.), 
           brown = c(1.,1.))}) %>% as.data.frame()




for(x in colnames(initcoords2)){
  initcoords2[1,x] = initcoords2[1,x] + rnorm(1) * 0.01
  initcoords2[2,x] = initcoords2[2,x] + rnorm(1) * 0.01
}


coords2 = layout_with_fr(g2, coords = t(initcoords2))
dim(initcoords2)  # Check dimensions

e2 = as_edgelist(g2, names = F)
coords2 = qgraph.layout.fruchtermanreingold(e2, vcount = vcount(g2), groups = markgroups2,
                                           init = t(initcoords2), repulse.rad = 10000)

plot(g2,
     mark.groups = markgroups2,
     mark.col = c(),
     mark.border = names(markgroups2),
     layout = coords2,
     vertex.size = 10,
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     edge.width = 1,
     edge.arrow.size = 0.2)
