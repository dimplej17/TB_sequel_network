##################################
# Co-expression Network Analysis #
##################################

source('1_libraries.R')
source('2_differential_gene_expression_analysis.R')

####################################################################################
# 1. Detecting and removing samples that are outliers in the dataset - perform PCA #
####################################################################################

temp <- numeric_data %>% as.data.frame() %>% column_to_rownames("UniProt")
# temp <- as.numeric(unlist(temp)) 

# Converting the data frame into a matrix
matrix_pca <- as.matrix(temp) 

# Removing coloumns containing missing values (coloumns -> samples)
matrix_pca <-  matrix_pca[ , colSums(is.na(matrix_pca))==0]

# Transpose the matrix so that the outliers are patients and not proteins
matrix_pca <- t(matrix_pca)

# Check for missing values
sum(is.na(matrix_pca)) 
dim(matrix_pca) 

# Standardize the Data. PCA requires data standardization. Standardize the variables in the dataset.
standardized_data <- scale(matrix_pca)

# Perform PCA
pca_result <- prcomp(standardized_data) 

# Check PCA results
summary(pca_result) 

# Visualise PCA results
autoplot(pca_result, label = TRUE) 

merged_OD_IF_data %>% 
  filter(!str_detect(SampleID, 'CONTROL_SAMPLE')) %>% 
  olink_pca_plot(df = .,
                 color_g = "QC_Warning", byPanel = TRUE)  

# Dendogram
# Compute distances and hierarchical clustering
dd <- dist(scale(matrix_pca), method = "euclidean")
hc <- hclust(dd, method = "ward.D2")

# Putting the labels at the same height (hang = -1)
plot(hc, hang = -1, cex = 0.6) 

###########################################
# 2. Choosing the soft-thresholding power #
###########################################

# Run pickSoftThreshold to identify the optimal power

sft <- pickSoftThreshold(
  matrix_pca, 
  dataIsExpr = TRUE,
  powerVector = c(seq(1, 10, by = 1), seq(2, 20, by = 2)), 
  networkType = "unsigned")

# Extract the results
fitIndices <- sft$fitIndices

# Plotting the results
par(mfrow = c(1, 2)) 

# Plot the scale-free topology fit index as a function of the soft-thresholding power
plot(fitIndices[,1], -sign(fitIndices[,3])*fitIndices[,2],
     xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = "Scale independence",
)
text(fitIndices[,1], -sign(fitIndices[,3])*fitIndices[,2], 
     labels = powers, 
     col = "red")
abline(h = 0.85, col = "blue", lty = 2)

# Plot the mean connectivity as a function of the soft-thresholding power
plot(fitIndices[,1], fitIndices[,5],
     xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean connectivity")
text(fitIndices[,1], fitIndices[,5], 
     labels = powers, 
     col = "red")

#################################################################################
# 3. Co-expression network construction and clustering the network into modules #
#################################################################################

# "In the signed hybrid network, only positive correlations among proteins expressions are included in the network and all negative correlations are removed."

# Calculate network adjacency

datExpr <- matrix_pca
softPower <- 4
adjacencyMat <- adjacency(datExpr, 
                          type = "signed hybrid", 
                          power = softPower)

# Convert adjacency matrix to Topological Overlap matrix (TOM)
TOM <- TOMsimilarity(adjacencyMat, TOMType="signed")

# Call the hierarchical clustering function
geneTree = flashClust(as.dist(1 - TOM), method = "average")

# Plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# We want large modules, so the minimum module size is set relatively high:
minModuleSize = 5
# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree, distM = 1 - TOM,
                            deepSplit = 4, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Plot the dendrogram and colors underneath
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

############################
# 4. Merging close modules #
############################

# "All modules are clustered according to the correlations between their module eigengenes and all modules which are close to each other, can be merged into only one module."

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)

# Cluster module eigengenes
METree = flashClust(as.dist(MEDiss), method = "average")

# Plot the result
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# Plot the cut line into the dendrogram
MEDissThres = 0.25
abline(h=MEDissThres, col = "red")

# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 0)

# The merged module colors
mergedColors = merge$colors

# Eigengenes of the new merged modules
mergedMEs = merge$newMEs

MEList = moduleEigengenes(datExpr, colors = mergedColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)

# Cluster module eigengenes
METree = flashClust(as.dist(MEDiss), method = "average")

# Plot the result
sizeGrWindow(12, 9)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Rename to moduleColors
moduleColors = mergedColors

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder) - 1
MEs = mergedMEs

############################################
# 5. Quantifying module-trait associations #
############################################

# "In this step, the amount of correlation between module eigengenes and the trait of interest is computed for all modules. For this dataset, the trait vector is defined as a vector including 0 and 1 values. For m0 and m6 samples, we consider 0 and 1 in the trait vector respectively."

# Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# Create the trait vector
sample_grps_df = merged_OD_IF_data %>% select(c('SampleID','Groups'))
duplicated(sample_grps_df)
# sample_grps_df[!duplicated(sample_grps_df),]
sample_grps_df = sample_grps_df[!duplicated(sample_grps_df),] %>% column_to_rownames('SampleID')
datTraits = sample_grps_df[rownames(MEs),]
moduleTraitCor = cor(MEs, as.numeric(datTraits), use = "all.obs")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
# par(mar = c(6, 8.5, 3, 3)); --> figure margins too large?
# Display the correlation values within a heatmap plot
dev.new(width = 10, height = 8)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = "TB Progression",
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               cex.main = 0.7,
               cex.axis=0.7,
               cex.lab=0.9,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# Define variable group containing the group column of datTrait
Group = as.data.frame(datTraits);
names(Group) = "Group"

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, as.numeric(Group$Group), use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Group), sep="");
names(GSPvalue) = paste("p.GS.", names(Group), sep="");
probes = colnames(datExpr)

geneInfo0 = data.frame(Proteins = probes,
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
modOrder = order(-abs(cor(MEs, as.numeric(Group$Group), use = "p")));
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Group));
geneInfo = geneInfo0[geneOrder, ]
# resultsDirectory <- "./Results/"
# if (!dir.exists(resultsDirectory)) {
# dir.create(resultsDirectory, recursive = TRUE)
# }

# write.csv(geneInfo, file = paste0(resultsDirectory,"geneInfo.csv"))
# save(geneInfo,file=paste0(resultsDirectory,"geneInfo.RData"))


filterGenes = function(x)
{
  module = x
  columnMM = paste0("MM.",module)
  columnPMM = paste0("p.MM.",module)
  geneInfo_Module = geneInfo[which(geneInfo$moduleColor==module),]
  selectedGenes = intersect(which(geneInfo_Module[,columnMM]>0),
                            which(geneInfo_Module[,columnPMM]<0.05))
  row.names(geneInfo_Module[selectedGenes,])
}
modules = lapply(modNames,function(x) filterGenes(x))
save(modules,modNames,file=paste0(resultsDirectory,"modules.RData"))

###########################################
# 6. Visualizing modules with PC1 vs. PC2 #
###########################################

# "In this section, we perform principal component analysis (PCA) on the expression profile of proteins of each module and use the first and the second principal components (PC1 and PC2) to show differences in protein levels."

# for the turquoise module
load(file = paste0(resultsDirectory,"modules.RData"))
module = "turquoise"
column = match(module, modNames);
moduleGenes = modules [[column]]

pca <- prcomp(standardized_data[,moduleGenes], scale=TRUE)
pca_x = pca$x[,1:2]
pca_x = as.data.frame(pca_x)
pca_x$SampleID = row.names(pca_x)
pca_x$Groups <- "m0"
pca_x$Groups[grepl("m6$", pca_x$SampleID)] <- "m6"
pca_x$Groups <- factor(pca_x$Groups, levels = c("m0", "m6"))

color_by_group <- function(pca_x, xaxis, yaxis) {
  p <- ggplot(pca_x, aes(y = .data[[yaxis]], x = .data[[xaxis]]))
  p + geom_point(aes(color=Groups,size=4)) + theme_bw() +
    theme(legend.title = element_blank(), legend.position="none", axis.title = element_blank(),
          axis.text = element_blank(),panel.border = element_rect(colour = module, fill=NA, size=3)) +
    labs(x=xaxis, y=yaxis)+
    geom_text(label=pca_x$SampleID)
}
color_by_group(pca_x, "PC1", "PC2") 
# pdf(paste0("../R Studio_LMU_Student Assistant_TB/Results/",module,".pdf"))
# color_by_group(pca_x, "PC1", "PC2") 
# dev.off()

# the same process is repeated for each module