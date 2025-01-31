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
library(EnhancedVolcano)




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
  filter(MaxLOD >= 0.2) %>%
  filter(!SampleID %in% controls_to_remove)

# Creating the coloumn "Groups" - m0 & m6
merged_OD_IF_data$Groups <- "m0"

# Assign 'm6' to sample IDs ending with 'm6'
merged_OD_IF_data$Groups[grepl("m6$", merged_OD_IF_data$SampleID)] <- "m6"

# Convert 'Groups' to a factor variable
merged_OD_IF_data$Groups <- factor(merged_OD_IF_data$Groups, levels = c("m0", "m6"))

# Create SubjectID coloumn from the SampleID
merged_OD_IF_data$SubjectID <- substr(merged_OD_IF_data$SampleID, 1, nchar(merged_OD_IF_data$SampleID) - 2)

# Add site to the dataframe
merged_OD_IF_data$Site <- NA

merged_OD_IF_data$Site[substr(merged_OD_IF_data$SubjectID, 1, 1) == "2"] <- "INS Mzb"

merged_OD_IF_data$Site[substr(merged_OD_IF_data$SubjectID, 1, 1) == "3"] <- "MMRC Tz"
# because SubjectIDs from Mzb begin from 2 and from Tz they begin from 3

merged_OD_IF_data$Site <- factor(merged_OD_IF_data$Site, levels = c("MMRC Tz", "INS Mzb"))

###############################################################################
# Differential gene expression analysis between m0 & m6 of Tanzanian patients #
###############################################################################
Tz_m0_m6_MWtest_results <- olink_wilcox(df = filter(merged_OD_IF_data, Site == "MMRC Tz"),
                                        variable = 'Groups',
                                        pair_id = 'SubjectID')  #135 (126 for MaxLOD >= 0.5)

olinkvolcano_plot <- olink_volcano_plot(Tz_m0_m6_MWtest_results,
                                    x_lab = "log2FC")


# Define the proteins of interest
proteins_of_interest <- c("CCL19", "BTC", "CCL20", "PRKRA", "CCL23", "PSMA1", "CXCL9",
                          "RARRES1", "CXCL10", "CXCL11", "EN-RAGE", "FGF21", "HGF", 
                          "MCP3", "OSM", "PD-L1", "SCF", "SIRT2", "TNFSF14", "VEGFA")

# Filter the dataset to keep only these proteins
Tz_m0_m6_MWtest_results$color <- ifelse(
  Tz_m0_m6_MWtest_results$Adjusted_pval >= 0.05, "grey", Tz_m0_m6_MWtest_results$Panel)

volcano_plot <- ggplot(Tz_m0_m6_MWtest_results, aes(x = estimate, y = -log10(Adjusted_pval))) +
  geom_point(data = subset(Tz_m0_m6_MWtest_results, Adjusted_pval >= 0.05 | abs(estimate) <= 1), 
             aes(color = NULL), color = "grey", size = 2, show.legend = FALSE) +
  geom_point(data = subset(Tz_m0_m6_MWtest_results, Adjusted_pval < 0.05 & abs(estimate) > 1), 
             aes(color = "blue"), size = 2, show.legend = FALSE) +  
  #geom_text_repel(data = Tz_m0_m6_MWtest_results,
                  #aes(label = Assay), size = 3, box.padding = 0.3, point.padding = 0.3) +  
  geom_text_repel(data = subset(Tz_m0_m6_MWtest_results, Adjusted_pval < 0.05 & abs(estimate) > 1),
                  aes(label = Assay), size = 3, box.padding = 0.3, point.padding = 0.3) +  # Add labels
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # P-value cut-off
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +  # Fold change cut-offs
  labs(x = "log2FC", y = "-log10(p-value)", title = "Expression Over Time") +
  theme_minimal()



print(volcano_plot)



