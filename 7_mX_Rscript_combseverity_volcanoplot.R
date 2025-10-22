#########################################
# Volcano plot - combined severity - m0 #
########################################

library(vroom)
library(OlinkAnalyze)
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(stringr)
library(readxl)
library(ggrepel)

basedir = '/home/obaranov/projects/TBSequel/DimplesProject/'

merged_OD_IF_data = vroom(paste0(basedir, '/FilteredData/mergedKits_wo3222.csv'))

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
  filter(Groups == "m6")

# Read the clinical data
clinical_data <- read_excel(paste0(basedir,"/RawData/CurTabForOlink.xlsx"))
clinical_data = clinical_data %>% filter(visit_number_unified == "Mth 06")
clinical_data['spiroSeverity'] = ''
mixsev = clinical_data$mixed_severity == 'severe'
mixsev[is.na(mixsev)] = FALSE
noimp = clinical_data$combined_severity == "none"
noimp[is.na(noimp)] = FALSE
clinical_data[mixsev,'spiroSeverity'] = 'mixedsevere'
clinical_data[noimp,'spiroSeverity'] = 'none'
# clinical_data = clinical_data[clinical_data$spiroSeverity != "",]

# Convert 'comb_sev_m6' to a factor variable
clinical_data$spiroSeverity <- factor(clinical_data$spiroSeverity, levels = c("none", "mixedsevere"))
clinical_data['broadSevGroup'] = 'good'
clinical_data[grepl('severe|moderate',clinical_data$combined_severity),'broadSevGroup'] = 'bad'

# Merge the two data frames by 'SampleID'
mergedkits_meta <- merge(merged_OD_IF_data, clinical_data[, c("person_id_complete", "broadSevGroup", "spiroSeverity")], by.x = "PersID", by.y = "person_id_complete", all.x = TRUE)



mergedkits_meta = mergedkits_meta[!is.na(mergedkits_meta$spiroSeverity),]
##########################################
# Differential gene expression analysis #
#########################################
m0_MWtest_results <- olink_wilcox(df = mergedkits_meta,
                                        variable = 'spiroSeverity')  

m0_MWtest_results = merge( m0_MWtest_results, 
      merged_OD_IF_data[c('Assay','kit')] %>% unique(), by = 'Assay')
m0_MWtest_results['logp'] = m0_MWtest_results$p.value %>%
          sapply(function(x){-1*log10(x)})


ggplot(data = m0_MWtest_results, aes(x = estimate, y = logp, label = Assay, color = kit)) + 
    geom_point() + geom_text_repel() 


olinkvolcano_plot <- olink_volcano_plot(m0_MWtest_results,
                                        x_lab = "log2FC")

olinkvolcano_plot +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "blue") +
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

