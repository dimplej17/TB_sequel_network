#########################################
# Volcano plot - combined severity - m0 #
########################################

library(OlinkAnalyze)
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(stringr)
library(readxl)
library(ggrepel)

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

# Create PersID coloumn from the SampleID
merged_OD_IF_data$PersID <- substr(merged_OD_IF_data$SampleID, 1, nchar(merged_OD_IF_data$SampleID) - 2)

# Add site to the dataframe
merged_OD_IF_data$Site <- NA

merged_OD_IF_data$Site[substr(merged_OD_IF_data$PersID, 1, 1) == "2"] <- "INS Mzb"

merged_OD_IF_data$Site[substr(merged_OD_IF_data$PersID, 1, 1) == "3"] <- "MMRC Tz"
# because PersIDs from Mzb begin from 2 and from Tz they begin from 3

merged_OD_IF_data$Site <- factor(merged_OD_IF_data$Site, levels = c("MMRC Tz", "INS Mzb"))

# Keeping samples only from Tanzania and m0
merged_OD_IF_data <- merged_OD_IF_data %>% 
  filter(Site == "MMRC Tz") %>%
  filter(Groups == "m0")

# Read the clinical data
clinical_data <- read_excel("TBS_ClinData.xlsx")

# Keeping samples only from Tanzania
clinical_data <- clinical_data %>% 
  filter(SiteID == "MMRC Tz")

# Limiting to the timepoint "m0" and removing variables not required for the volcano plot
clinical_data <- clinical_data[ , c(1, 59)] # 1 -PersID; 59 - comb_sev_Mth_06; 56 - comb_sev_Baseline
clinical_data <- clinical_data %>%
  rename(comb_sev_m6 = comb_sev_Mth_06)

# Grouping Normal & Mild under "Good" and Moderate & Severe under "Bad"
clinical_data$comb_sev_m6[grepl("Moderate|Severe", clinical_data$comb_sev_m6)] <- "Bad"
clinical_data$comb_sev_m6[grepl("Normal|Mild", clinical_data$comb_sev_m6)] <- "Good"

# Convert 'comb_sev_m6' to a factor variable
clinical_data$comb_sev_m6 <- factor(clinical_data$comb_sev_m6, levels = c("Bad", "Good"))

# Merge the two data frames by 'SampleID'
merged_OD_IF_data <- merge(merged_OD_IF_data, clinical_data[, c("PersID", "comb_sev_m6")], by = "PersID", all.x = TRUE)

##########################################
# Differential gene expression analysis #
#########################################
m0_MWtest_results <- olink_wilcox(df = merged_OD_IF_data,
                                        variable = 'comb_sev_m6')  

olinkvolcano_plot <- olink_volcano_plot(m0_MWtest_results,
                                        x_lab = "log2FC")
olinkvolcano_plot +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "red") +
  geom_text_repel(data = subset(m0_MWtest_results, abs(estimate) > 0.5), aes(label = Assay), size = 3, box.padding = 0.3  , point.padding = 0.3) +
  labs(title = "Differential Expression Analysis at m0 by Lung Function Severity")





olinkvolcano_plot <- ggplot(m0_MWtest_results, aes(x = estimate, y = -log10(Adjusted_pval))) +
  geom_point(data = subset(m0_MWtest_results, abs(estimate) <= 0.5), color = "grey", size = 2, show.legend = FALSE) +
  # Blue points for significant points and large estimates
  geom_point(data = subset(m0_MWtest_results, abs(estimate) > 0.5), color = "blue", size = 2, show.legend = FALSE) +    
  #geom_text_repel(data = m0_MWtest_results, aes(label = Assay), size = 3, box.padding = 0.3, point.padding = 0.3) +  
  geom_text_repel(data = subset(m0_MWtest_results, abs(estimate) > 0.5), aes(label = Assay), size = 3, box.padding = 0.3, point.padding = 0.3) +  # Add labels
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # P-value cut-off
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "red") +  # Fold change cut-offs
  labs(x = "log2FC", y = "-log10(p-value)", title = "Differential Expression Analysis at m0 by Lung Function Severity") +
  theme_minimal()


print(olinkvolcano_plot)




library(writexl)
write_xlsx(m0_MWtest_results, "/Users/dimplejanardhan/Downloads/LMU Clinic_Student Assistant/R Studio_LMU_Student Assistant_TB/m0_MWtest_results.xlsx")

