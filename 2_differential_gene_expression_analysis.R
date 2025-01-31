#####################################
# Differential Expression Analysis #
####################################

source('1_libraries.R')

# Reading the data using the read_NPX()
OD_data <- read_NPX("/Users/dimplejanardhan/Downloads/LMU_Clinic_Student Assistant/IMPORTANT_files_for_Rproject/TB_Sequel_OD_NPX.xlsx")
IF_data <- read_NPX("/Users/dimplejanardhan/Downloads/LMU_Clinic_Student Assistant/IMPORTANT_files_for_Rproject/TB_Sequel_Inflammation_NPX.xlsx")

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
# save(MWtest_results, file = "MWtest_results_unpaired.RData")
MWtest_results %>% filter(Adjusted_pval < 0.05) %>% select(UniProt) %>% print(n = 92)

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
                                      pair_id = 'SubjectID')
paired_MWtest_results %>% filter(Adjusted_pval < 0.05) %>% select(UniProt) %>% print(n = 97)

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

merged_OD_IF_data$Site <- factor(merged_OD_IF_data$Site, levels = c("MMRC Tz", "INS Mzb"))

# Differential gene expression analysis between m0 & m6 of Tanzanian patients
Tz_m0_m6_MWtest_results <- olink_wilcox(df = filter(merged_OD_IF_data, Site == "MMRC Tz"),
                                        variable = 'Groups',
                                        pair_id = 'SubjectID')

Tz_m0_m6_MWtest_results %>% filter(Adjusted_pval < 0.05) %>% select(UniProt) %>% print(n = 97)

# Volcano Plot for MW test of m0 & m6 Tanzanian patients
olink_volcano_plot(Tz_m0_m6_MWtest_results,
                   x_lab = "log2FC")
# p-value distribution
ggplot(data = Tz_m0_m6_MWtest_results) +
  geom_histogram(mapping = aes(x = p.value))

# Differential gene expression analysis between Tanzania & Mozambique of m0 patients
m0_Tz_Mzb_MWtest_results <- olink_wilcox(df = filter(merged_OD_IF_data, Groups == "m0"),
                                         variable = 'Site')

m0_Tz_Mzb_MWtest_results %>% filter(Adjusted_pval < 0.05) %>% select(UniProt) %>% print(n = 28) 

# Volcano Plot for unpaired Tz & Mzb MW U test
olink_volcano_plot(m0_Tz_Mzb_MWtest_results,
                   x_lab = "log2FC")

# p-value distribution
ggplot(data = m0_Tz_Mzb_MWtest_results) +
  geom_histogram(mapping = aes(x = p.value))