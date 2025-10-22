#####################################
# Differential Expression Analysis #
####################################

library(OlinkAnalyze)
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(gridExtra)

basedir = '/home/obaranov/projects/TBSequel/DimplesProject/'

# Reading the data using the read_NPX()
OD_data <- read_NPX(paste0(basedir,"RawData/TB_Sequel_OD_NPX.xlsx"))
IF_data <- read_NPX(paste0(basedir,"RawData/TB_Sequel_Inflammation_NPX.xlsx"))
OD_data['kit'] = 'OD'
IF_data['kit'] = 'IF'


#checking LOD fractions
# OD_data['NPX'] - OD_data['PlateLOD'] < 0
# around 8000 rows have NPX < LOD

# Merging/stacking the 2 datasets
merged_OD_IF_data <- rbind(OD_data, IF_data)

# Removing all controls and observations that have not passed the QC, that are below the LOD, and observations with no UniProt IDs. 
controls_to_remove <- c("CBctrl","CGctrl","IPC","IPC_02","IPC_03","IPC-2","IPC-3","KHctrl","NC","NC_02","NC_03","NC-2","NC-3","SC","SC_02","SC-2","SRctrl","XXctrl")


merged_OD_IF_data <- merged_OD_IF_data %>% 
  # filter(QC_Warning == "Pass") %>%   # commented out to make it consistent with K+N analysis
  filter(UniProt != "-") %>%
  # filter(MaxLOD >= 0.5) %>%   # commented out as Olink is now suggesting not to remove the below LOD values
  filter(!SampleID %in% controls_to_remove)

# Creating the coloumn "Groups" - m0 & m6
merged_OD_IF_data$Groups <- "m0"

# Assign 'm6' to sample IDs ending with 'm6'
merged_OD_IF_data$Groups[grepl("m6$", merged_OD_IF_data$SampleID)] <- "m6"

# Convert 'Groups' to a factor variable
merged_OD_IF_data$Groups <- factor(merged_OD_IF_data$Groups, levels = c("m0", "m6"))

write_delim(merged_OD_IF_data, paste0(basedir, '/FilteredData/mergedKits.csv'), delim = '\t')



# Unpaired MW test 
MWtest_results <- olink_wilcox(df = merged_OD_IF_data,
                               variable = 'Groups')
saveRDS(MWtest_results, file = paste0(basedir,"/FilteredData/MWtest_results_unpaired.RData"))

# how is the significance distributed across the kits?
# sigs = MWtest_results %>% filter(Adjusted_pval < 0.05) %>% pull('Assay')
# merged_OD_IF_data %>% filter(Assay %in% sigs) %>% select(c('Assay','kit')) %>% unique() %>% pull('kit') %>% table
# ~ 65 per kit, similar number for IF & OD


# Volcano Plot

MWtest_volcano <- olink_volcano_plot(MWtest_results,
                                     x_lab = "log2FC")

# Create SubjectID coloumn from the SampleID
merged_OD_IF_data$SubjectID <- substr(merged_OD_IF_data$SampleID, 1, nchar(merged_OD_IF_data$SampleID) - 2)

# Paired MW test
paired_MWtest_results <- olink_wilcox(df = merged_OD_IF_data,
                                      variable = 'Groups',
                                      pair_id = 'SubjectID')
# Subjects that don't have a pair are removed automatically
saveRDS(paired_MWtest_results, file = paste0(basedir,"/FilteredData/MWtest_results_paired.RData"))



# Paired MW test Volcano Plot
pairedMWtest_volcano <- olink_volcano_plot(paired_MWtest_results,
                                           x_lab = "log2FC")

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

merged_OD_IF_data$Site <- factor(merged_OD_IF_data$Site, levels = c("MMRC Tz", "INS Mzb"))



# Differential gene expression analysis between Tanzania & Mozambique of m0 patients
m0_Tz_Mzb_MWtest_results <- olink_wilcox(df = filter(merged_OD_IF_data, Groups == "m0"),
                                         variable = 'Site')

m0_Tz_Mzb_MWtest_results %>% filter(Adjusted_pval < 0.05)

# Volcano Plot for unpaired Tz & Mzb MW U test
olink_volcano_plot(m0_Tz_Mzb_MWtest_results,
                   x_lab = "log2FC")

# p-value distribution
ggplot(data = m0_Tz_Mzb_MWtest_results) +
  geom_histogram(mapping = aes(x = p.value))

# overlap of genes that are DEGs in paired m0 vs m6 analysis and the Tz vs Mz analysis


sigTvM = m0_Tz_Mzb_MWtest_results %>% filter(Adjusted_pval < 0.05) %>% pull(Assay)
sig0v6 = paired_MWtest_results %>% filter(Adjusted_pval < 0.05) %>% pull(Assay)

drop = intersect(sig0v6, sigTvM)

paired_MWtest_results['difBtwSites'] = FALSE
paired_MWtest_results[paired_MWtest_results$Assay %in% drop,'difBtwSites'] = TRUE
write_delim(paired_MWtest_results, paste0(basedir, "/OutputData/pairedM0vsM6_wilcox.csv"), delim = '\t')