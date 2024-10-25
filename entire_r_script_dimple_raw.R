# Load the necessary libraries
library(OlinkAnalyze)
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(stringr)
library(readxl)
library(limma)
library(ggrepel)
library(RColorBrewer)
library(ggvenn)
library(reshape)
library(reshape2)
library(cowplot)
library(gridExtra)
library(grid)
library(clusterProfiler)
library(qgraph)
library(igraph)
library(GEOquery)
library(EnhancedVolcano) #i don't think i used this




# 1. Preprocessing

# 2. Differential Expression Analysis

# Reading the data using the read_NPX()
OD_data <- read_NPX("/Users/dimplejanardhan/Downloads/LMU Clinic_Student Assistant/R Studio_LMU_Student Assistant_TB/TB Sequel OD_NPX.xlsx")
IF_data <- read_NPX("/Users/dimplejanardhan/Downloads/LMU Clinic_Student Assistant/R Studio_LMU_Student Assistant_TB/TB Sequel Inflammation_NPX.xlsx")

# Merging/stacking the 2 datasets
merged_OD_IF_data <- rbind(OD_data_pass, IF_data_pass)

# Removing all controls and observations that have not passed the QC, that are below the LOD, and observations with no UniProt IDs. 
controls_to_remove <- c("CGctrl", "IPC", "IPC-3", "KHctrl", "SC", "SC_02", "SRctrl", "IPC_03", "CBctrl")

merged_OD_IF_data <- merged_OD_IF_data %>% 
  filter(QC_Warning == "Pass") %>%
  filter(UniProt != "-") %>%
  filter(MaxLOD >= 0.5) %>%
  filter(!SampleID %in% controls_to_remove)

# Creating the coloumn "Groups" - m0 & m6
merged_OD_IF_data$Groups <- "m0"

# Assign 'm6' to sample IDs ending with 'm6'
merged_OD_IF_data$Groups[grepl("m6$", merged_OD_IF_data$SampleID)] <- "m6"

# Convert 'Groups' to a factor variable
merged_OD_IF_data$Groups <- factor(merged_OD_IF_data$Groups, levels = c("m0", "m6"))

# Unpaired MW test 
MWtest_results <- olink_wilcox(df = merged_OD_IF_data,
            variable = 'Groups')
MWtest_results %>% filter(Adjusted_pval < 0.05) %>% select(UniProt) %>% print(n = 92) # 92

# Volcano Plot
# Select names of proteins to show
top_10_name_unpairedMWtest <- MWtest_results %>%
  slice_head(n = 10) %>%
  pull(OlinkID)

MWtest_volcano <- olink_volcano_plot(MWtest_results,
                                    x_lab = "log2FC",
                                    olinkid_list = top_10_name_unpairedMWtest)

# Create SubjectID coloumn from the SampleID
merged_OD_IF_data$SubjectID <- substr(merged_OD_IF_data$SampleID, 1, nchar(merged_OD_IF_data$SampleID) - 2)

# Paired MW test
paired_MWtest_results <- olink_wilcox(df = merged_OD_IF_data,
                             variable = 'Groups',
                             pair_id = 'SubjectID')  # 126
paired_MWtest_results %>% filter(Adjusted_pval < 0.05) %>% select(UniProt) %>% print(n = 97) # 97

# Paired MW test Volcano Plot
top_10_name_pairedMWtest <- paired_MWtest_results %>%
  slice_head(n = 10) %>%
  pull(OlinkID)

pairedMWtest_volcano <- olink_volcano_plot(paired_MWtest_results,
                                          x_lab = "log2FC",
                                          olinkid_list = top_10_name_pairedMWtest)

# Both MW tests Volcano Plots
grid.arrange(MWtest_volcano, pairedMWtest_volcano, ncol = 2)

# p-value distribution
ggplot(data = MWtest_results) +
  geom_histogram(mapping = aes(x = p.value))

ggplot(data = paired_MWtest_results) +
  geom_histogram(mapping = aes(x = p.value))

# Add site to the dataframe
merged_OD_IF_data$Site <- NA

merged_OD_IF_data$Site[substr(merged_OD_IF_data$SubjectID, 1, 1) == "2"] <- "INS Mzb"

merged_OD_IF_data$Site[substr(merged_OD_IF_data$SubjectID, 1, 1) == "3"] <- "MMRC Tz"
# because SubjectIDs from Mzb begin from 2 and from Tz they begin from 3

merged_OD_IF_data$Site <- factor(merged_OD_IF_data$Site, levels = c("MMRC Tz", "INS Mzb"))
# factor: The function factor is used to encode a vector as a factor 

#########################################################################################
# Differential gene expression analysis between m0 & m6 of Tanzanian patients
Tz_m0_m6_MWtest_results <- olink_wilcox(df = filter(merged_OD_IF_data, Site == "MMRC Tz"),
                                             variable = 'Groups',
                                             pair_id = 'SubjectID')  # 126




Tz_m0_m6_MWtest_results %>% filter(Adjusted_pval < 0.05) %>% select(UniProt) %>% print(n= 97) # 97

# Volcano Plot for MW test of m0 & m6 Tanzanian patients
volcano_plott <- olink_volcano_plot(Tz_m0_m6_MWtest_results,
                   x_lab = "log2FC")




# Define the proteins of interest
proteins_of_interest <- c("CCL19", "BTC", "CCL20", "PRKRA", "CCL23", "PSMA1", "CXCL9",
                          "RARRES1", "CXCL10", "CXCL11", "EN-RAGE", "FGF21", "HGF", 
                          "MCP3", "OSM", "PD-L1", "SCF", "SIRT2", "TNFSF14", "VEGFA")

# Filter the dataset to keep only these proteins
Tz_m0_m6_MWtest_results$color <- ifelse(
  Tz_m0_m6_MWtest_results$Adjusted_pval >= 0.05, "grey", Tz_m0_m6_MWtest_results$Panel)

volcano_plott <- ggplot(Tz_m0_m6_MWtest_results, aes(x = estimate, y = -log10(Adjusted_pval))) +
  geom_point(data = subset(Tz_m0_m6_MWtest_results, Adjusted_pval >= 0.05 | abs(estimate) <= 1), 
             aes(color = NULL), color = "grey", size = 2, show.legend = FALSE) +
  geom_point(data = subset(Tz_m0_m6_MWtest_results, Adjusted_pval < 0.05 & abs(estimate) > 1), 
             aes(color = Panel), size = 2) +  
  scale_color_manual(values = c("blue","orange")) +
  geom_text_repel(data = subset(Tz_m0_m6_MWtest_results, Assay %in% proteins_of_interest),
                  aes(label = Assay), size = 3, box.padding = 0.3, point.padding = 0.3) +  # Add labels
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # P-value cut-off
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +  # Fold change cut-offs
  labs(x = "log2FC", y = "-log10(p-value)", title = "Expression Over Time", color = "Panel") +
  theme_minimal()



print(volcano_plott)


#########################################################################################


# p-value distribution
ggplot(data = Tz_m0_m6_MWtest_results) +
  geom_histogram(mapping = aes(x = p.value))

# Differential gene expression analysis between Tanzania & Mozambique of m0 patients
m0_Tz_Mzb_MWtest_results <- olink_wilcox(df = filter(merged_OD_IF_data, Groups == "m0"),
                                        variable = 'Site')

m0_Tz_Mzb_MWtest_results %>% filter(Adjusted_pval < 0.05) %>% select(UniProt) %>% print(n = 28) # 28

# Volcano Plot for unpaired Tz & Mzb MW U test
olink_volcano_plot(m0_Tz_Mzb_MWtest_results,
                   x_lab = "log2FC")

# p-value distribution
ggplot(data = m0_Tz_Mzb_MWtest_results) +
  geom_histogram(mapping = aes(x = p.value))

# 3. Pathway Enrichment 

numeric_data <- pivot_wider(merged_OD_IF_data, names_from = 'SampleID', id_cols = 'UniProt', values_from = c('NPX'))

# Removing variables not required for pathway enrichment
npx_df <- merged_OD_IF_data[ , -c(6, 8, 9, 14, 15)]

utest_results <- olink_wilcox(
  df = npx_df,
  variable = "Groups")

try({
  gsea_results <- olink_pathway_enrichment(data = merged_OD_IF_data, test_results = utest_results)
  ora_results <- olink_pathway_enrichment(
    data = merged_OD_IF_data,
    test_results = utest_results, method = "ORA")
}, silent = TRUE)

olink_pathway_heatmap(gsea_results, utest_results)
olink_pathway_heatmap(ora_results, utest_results, method = "ORA")

olink_pathway_visualization(gsea_results) 
olink_pathway_visualization(ora_results, method = "ORA") 

# Saving GSEA & ORA results into csv files
write.table(gsea_results,file="gsea_results.csv",sep='\t')
write.table(ora_results,file="ora_results.csv",sep='\t')

# 4. Co-expression Network Analysis
# 4.1 Detecting and removing samples that are outliers in the dataset - perform PCA

library(ggfortify)
library(stats) 

temp <- numeric_data %>% as.data.frame() %>% column_to_rownames("UniProt")
# temp <- as.numeric(unlist(temp)) 

# Converting the data frame into a matrix
matrix_pca <- as.matrix(temp) 

# Removing coloumns containing missing values (coloumns -> samples)
matrix_pca <-  matrix_pca[ , colSums(is.na(matrix_pca))==0]

# Transpose the matrix so that the outliers are patients and not proteins
matrix_pca <- t(matrix_pca)

# Check for missing values
sum(is.na(matrix_pca)) # should be zero?
dim(matrix_pca) # should be 118 126?

# Standardize the Data. PCA requires data standardization. Standardize the variables in the dataset.
standardized_data <- scale(matrix_pca)

# Perform PCA
pca_result <- prcomp(standardized_data) # prcomp or princomp? Both functions implement PCA, however the princomp() function uses the spectral decomposition approach, whereas the prcomp() function uses singular value decomposition (SVD)?????

# Check PCA results
summary(pca_result) # This output provides the standard deviations of each principal component and the proportion of total variance explained?

# Visualise PCA results
autoplot(pca_result, label = TRUE) # no clear outliers

# OR - use Olink PCA plot - BUT there are no controls in our dataset???????
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
# no clear outliers, don't remove the 2 you thought could be outliers (2516D1m0, 3036R1m0) because they're not sticking out of the main branch or anything like that so don't remove anything.

# 4.2 Choosing the soft-thresholding power
# "Most of the biological networks have the scale-free topology. So, we try different powers factor in this step to choose the most suitable power factor which makes a network with the scale-free property."

# Run pickSoftThreshold to identify the optimal power
library(WGCNA)

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

# 4.2 Co-expression network construction and clustering the network into modules

# "We can construct a signed hybrid network using different correlation measures, like pearson and bicor, to see which network makes better sense in biology. In the signed hybrid network, only positive correlations among proteins expressions are included in the network and all negative correlations are removed."

# Calculate network adjacency
library(xfun)
library(WGCNA)

datExpr <- matrix_pca
softPower <- 4
adjacencyMat <- adjacency(datExpr, 
          type = "signed hybrid", 
          power = softPower)

# Convert adjacency matrix to Topological Overlap matrix (TOM)
TOM <- TOMsimilarity(adjacencyMat, TOMType="signed")

# Call the hierarchical clustering function
library(flashClust)
geneTree = flashClust(as.dist(1 - TOM), method = "average")

# Plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high:
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

# 4.3 Merging close modules
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
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

# 4.4 Quantifying module-trait associations
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

write.csv(geneInfo, file = paste0(resultsDirectory,"geneInfo.csv"))
save(geneInfo,file=paste0(resultsDirectory,"geneInfo.RData"))


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


# 4.5 Visualizing modules with PC1 vs. PC2
# "In this section, we perform principal component analysis (PCA) on the expression profile of proteins of each module and used the first and the second principal components (PC1 and PC2) to show differences in protein levels."

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

# # for the red module
# load(file = paste0(resultsDirectory,"modules.RData"))
# module = "red"
# column = match(module, modNames);
# moduleGenes = modules [[column]]
# 
# pca <- prcomp(standardized_data[,moduleGenes], scale=TRUE)
# pca_x = pca$x[,1:2]
# pca_x = as.data.frame(pca_x)
# pca_x$SampleID = row.names(pca_x)
# pca_x$Groups <- "m0"
# pca_x$Groups[grepl("m6$", pca_x$SampleID)] <- "m6"
# pca_x$Groups <- factor(pca_x$Groups, levels = c("m0", "m6"))
# 
# color_by_group <- function(pca_x, xaxis, yaxis) {
#   p <- ggplot(pca_x, aes(y = .data[[yaxis]], x = .data[[xaxis]]))
#   p + geom_point(aes(color=Groups,size=4)) + theme_bw() +
#     theme(legend.title = element_blank(), legend.position="none", axis.title = element_blank(),
#           axis.text = element_blank(),panel.border = element_rect(colour = module, fill=NA, size=3)) +
#     labs(x=xaxis, y=yaxis)+
#     geom_text(label=pca_x$SampleID)
# }
# color_by_group(pca_x, "PC1", "PC2") 
# pdf(paste0("../R Studio_LMU_Student Assistant_TB/Results/",module,".pdf"))
# color_by_group(pca_x, "PC1", "PC2") 
# dev.off()
# 
# 
# # for the green module
# load(file = paste0(resultsDirectory,"modules.RData"))
# module = "green"
# column = match(module, modNames);
# moduleGenes = modules [[column]]
# 
# pca <- prcomp(standardized_data[,moduleGenes], scale=TRUE)
# pca_x = pca$x[,1:2]
# pca_x = as.data.frame(pca_x)
# pca_x$SampleID = row.names(pca_x)
# pca_x$Groups <- "m0"
# pca_x$Groups[grepl("m6$", pca_x$SampleID)] <- "m6"
# pca_x$Groups <- factor(pca_x$Groups, levels = c("m0", "m6"))
# 
# color_by_group <- function(pca_x, xaxis, yaxis) {
#   p <- ggplot(pca_x, aes(y = .data[[yaxis]], x = .data[[xaxis]]))
#   p + geom_point(aes(color=Groups,size=4)) + theme_bw() +
#     theme(legend.title = element_blank(), legend.position="none", axis.title = element_blank(),
#           axis.text = element_blank(),panel.border = element_rect(colour = module, fill=NA, size=3)) +
#     labs(x=xaxis, y=yaxis)+
#     geom_text(label=pca_x$SampleID)
# }
# color_by_group(pca_x, "PC1", "PC2") 
# pdf(paste0("../R Studio_LMU_Student Assistant_TB/Results/",module,".pdf"))
# color_by_group(pca_x, "PC1", "PC2") 
# dev.off()
# 
# # for the black module
# load(file = paste0(resultsDirectory,"modules.RData"))
# module = "black"
# column = match(module, modNames);
# moduleGenes = modules [[column]]
# 
# pca <- prcomp(standardized_data[,moduleGenes], scale=TRUE)
# pca_x = pca$x[,1:2]
# pca_x = as.data.frame(pca_x)
# pca_x$SampleID = row.names(pca_x)
# pca_x$Groups <- "m0"
# pca_x$Groups[grepl("m6$", pca_x$SampleID)] <- "m6"
# pca_x$Groups <- factor(pca_x$Groups, levels = c("m0", "m6"))
# 
# color_by_group <- function(pca_x, xaxis, yaxis) {
#   p <- ggplot(pca_x, aes(y = .data[[yaxis]], x = .data[[xaxis]]))
#   p + geom_point(aes(color=Groups,size=4)) + theme_bw() +
#     theme(legend.title = element_blank(), legend.position="none", axis.title = element_blank(),
#           axis.text = element_blank(),panel.border = element_rect(colour = module, fill=NA, size=3)) +
#     labs(x=xaxis, y=yaxis)+
#     geom_text(label=pca_x$SampleID)
# }
# color_by_group(pca_x, "PC1", "PC2") 
# pdf(paste0("../R Studio_LMU_Student Assistant_TB/Results/",module,".pdf"))
# color_by_group(pca_x, "PC1", "PC2") 
# dev.off()
# 
# # for the brown module
# load(file = paste0(resultsDirectory,"modules.RData"))
# module = "brown"
# column = match(module, modNames);
# moduleGenes = modules [[column]]
# 
# pca <- prcomp(standardized_data[,moduleGenes], scale=TRUE)
# pca_x = pca$x[,1:2]
# pca_x = as.data.frame(pca_x)
# pca_x$SampleID = row.names(pca_x)
# pca_x$Groups <- "m0"
# pca_x$Groups[grepl("m6$", pca_x$SampleID)] <- "m6"
# pca_x$Groups <- factor(pca_x$Groups, levels = c("m0", "m6"))
# 
# color_by_group <- function(pca_x, xaxis, yaxis) {
#   p <- ggplot(pca_x, aes(y = .data[[yaxis]], x = .data[[xaxis]]))
#   p + geom_point(aes(color=Groups,size=4)) + theme_bw() +
#     theme(legend.title = element_blank(), legend.position="none", axis.title = element_blank(),
#           axis.text = element_blank(),panel.border = element_rect(colour = module, fill=NA, size=3)) +
#     labs(x=xaxis, y=yaxis)+
#     geom_text(label=pca_x$SampleID)
# }
# color_by_group(pca_x, "PC1", "PC2") 
# pdf(paste0("../R Studio_LMU_Student Assistant_TB/Results/",module,".pdf"))
# color_by_group(pca_x, "PC1", "PC2") 
# dev.off()
# 
# # for the blue module
# load(file = paste0(resultsDirectory,"modules.RData"))
# module = "blue"
# column = match(module, modNames);
# moduleGenes = modules [[column]]
# 
# pca <- prcomp(standardized_data[,moduleGenes], scale=TRUE)
# pca_x = pca$x[,1:2]
# pca_x = as.data.frame(pca_x)
# pca_x$SampleID = row.names(pca_x)
# pca_x$Groups <- "m0"
# pca_x$Groups[grepl("m6$", pca_x$SampleID)] <- "m6"
# pca_x$Groups <- factor(pca_x$Groups, levels = c("m0", "m6"))
# 
# color_by_group <- function(pca_x, xaxis, yaxis) {
#   p <- ggplot(pca_x, aes(y = .data[[yaxis]], x = .data[[xaxis]]))
#   p + geom_point(aes(color=Groups,size=4)) + theme_bw() +
#     theme(legend.title = element_blank(), legend.position="none", axis.title = element_blank(),
#           axis.text = element_blank(),panel.border = element_rect(colour = module, fill=NA, size=3)) +
#     labs(x=xaxis, y=yaxis)+
#     geom_text(label=pca_x$SampleID)
# }
# color_by_group(pca_x, "PC1", "PC2") 
# pdf(paste0("../R Studio_LMU_Student Assistant_TB/Results/",module,".pdf"))
# color_by_group(pca_x, "PC1", "PC2") 
# dev.off()
# 
# # for the grey module
# load(file = paste0(resultsDirectory,"modules.RData"))
# module = "grey"
# column = match(module, modNames);
# moduleGenes = modules [[column]]
# 
# pca <- prcomp(standardized_data[,moduleGenes], scale=TRUE)
# pca_x = pca$x[,1:2]
# pca_x = as.data.frame(pca_x)
# pca_x$SampleID = row.names(pca_x)
# pca_x$Groups <- "m0"
# pca_x$Groups[grepl("m6$", pca_x$SampleID)] <- "m6"
# pca_x$Groups <- factor(pca_x$Groups, levels = c("m0", "m6"))
# 
# color_by_group <- function(pca_x, xaxis, yaxis) {
#   p <- ggplot(pca_x, aes(y = .data[[yaxis]], x = .data[[xaxis]]))
#   p + geom_point(aes(color=Groups,size=4)) + theme_bw() +
#     theme(legend.title = element_blank(), legend.position="none", axis.title = element_blank(),
#           axis.text = element_blank(),panel.border = element_rect(colour = module, fill=NA, size=3)) +
#     labs(x=xaxis, y=yaxis)+
#     geom_text(label=pca_x$SampleID)
# }
# color_by_group(pca_x, "PC1", "PC2") 
# pdf(paste0("../R Studio_LMU_Student Assistant_TB/Results/",module,".pdf"))
# color_by_group(pca_x, "PC1", "PC2") 
# dev.off()

# 5. Network Visualisation using igraph

thr = 0.65
adjacencyMat <- adjacency(datExpr, power = 1, type="signed hybrid");
m = lower.tri(adjacencyMat,diag = FALSE)
m = melt(m)
net = melt(adjacencyMat)
net = net[which(m$value==TRUE),]
names(net) = c("Protein1","Protein2","Correlation")
net = net[which(net$Correlation>=thr),]

unique(net$Protein2)

library(igraph)


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
  pal <- colorRampPalette(c(negcol, midcol, poscol))(101)
  xx <- round(range01(x), 2) * 100 + 1 # added 1 as indexing from 0 does not work in R
  return(pal[xx])
}

########## main

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
?layout_with_fr

e = as_edgelist(g, names = F)
coords = qgraph.layout.fruchtermanreingold(e, vcount = vcount(g), groups = markgroups,
                                           init = t(initcoords), repulse.rad = 10000)
# init = t(initcoords)

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

#############
# graphviz2 #
#############

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

adjlist = data.frame(net, row.names = NULL)

lfcTab = data.frame(node = MWtest_results$UniProt, lfc = MWtest_results$estimate)


moduleTab = data.frame(node = geneInfo$Proteins, 
                       mod = geneInfo$moduleColor) %>%
  column_to_rownames('node')

g <- graph_from_data_frame(adjlist, directed = TRUE)

moduleTab = moduleTab[V(g)$name,"mod", drop=FALSE]
moduleTab$mod <- as.factor(moduleTab$mod)

g = set_vertex_attr(g, "color", value = lfcTab[V(g),] %>% column_to_rownames('node'))


markgroups = moduleTab %>% rownames_to_column("node") %>%
  group_by(mod) %>% 
  group_map(function(x,i){return(x$node)})
names(markgroups) = moduleTab %>% rownames_to_column("node") %>%
  group_by(mod) %>% 
  group_map(function(x,i){ return(as.character(i %>% pull))})

sublfcTab = lfcTab %>% column_to_rownames("node")
V(g)$color <- getColor(sublfcTab[ V(g)$name,"lfc"])

plot(g, mark.groups = markgroups, mark.border = names(markgroups))

# using additional library
library(qgraph)


coords = moduleTab[V(g)$name,'mod'] %>% sapply( function(x){
  switch(as.character(x),
         turquoise = c(0,0),
         red = c(1,0),
         brown =c(1,1),
         blue = c(0,1))
})  %>% t() %>% as.data.frame()


for(row in rownames(coords)){
  coords[row, 1] = coords[row, 1] + runif(1,-0.2,0.2)
  coords[row, 2] = coords[row, 2] + runif(1, -0.2, 0.2)
}

coords = coords %>% as.matrix

e <- as_edgelist(g, names = FALSE)
#  niter = 40, max.delta = 0.01 are the parameters that prevent the nodes from moving too far out
lo = qgraph.layout.fruchtermanreingold(e, vcount = vcount(g), init = coords, niter = 40, max.delta = 0.01)
?qgraph.layout.fruchtermanreingold
#pdf("./GraphViz.pdf")
p = plot(g, layout = lo, mark.groups = markgroups, mark.col = c(), mark.border = names(markgroups))
#dev.off()


# Network visaulisation with Protein names instead of UniProtIDs


############################
# 5. Signature Refinement #
###########################

# 5.1 Lower respiratory infections datasets

# "We performed differential expression analysis using multiple transcriptomic datasets including either viral or bacterial lower respiratory tract infections (LRTI) or systemic inflammatory response syndrome (SIRS). To this end, we used three transcriptomic datasets, GSE42026, GSE40012 and GSE60244, each containing various types of lower respiratory infections."


library(GEOquery)#directlycall
eSet <- getGEO("GSE42026",
               destdir = '.',
               getGPL = F)

GSE42026 <- read.delim("GSE42026_series_matrix.txt", header = TRUE, sep = '\t', dec = ".")
View(GSE42026)

varnames = grepl("^!", GSE42026[,1])
GSE42026[varnames,1]








