#######################
# Pathway Enrichment #
######################

library('vroom')
library(tidyverse)

basedir = '/home/obaranov/projects/TBSequel/DimplesProject/'

data = OD_data <- read_NPX(paste0(basedir,"RawData/TB_Sequel_OD_NPX.xlsx"))
DEGs = vroom(paste0(basedir, '/OutputData/pairedM0vsM6_wilcox.csv'))

gsea_results <- olink_pathway_enrichment(data = data, test_results = DEGs)
ora_results <- olink_pathway_enrichment(data = data,
    test_results = DEGs, method = "ORA")
write_delim(gsea_results, paste0(basedir, "/OutputData/pairedM0vsM6_gsea.csv"), delim = '\t')
write_delim(ora_results, paste0(basedir, "/OutputData/pairedM0vsM6_ora.csv"), delim = '\t')


olink_pathway_heatmap(gsea_results, DEGs)
olink_pathway_heatmap(ora_results, DEGs, method = "ORA")


olink_pathway_visualization(gsea_results) 
olink_pathway_visualization(ora_results, method = "ORA") 

# Saving GSEA & ORA results into csv files
# write.table(gsea_results,file="gsea_results.csv",sep='\t')
# write.table(ora_results,file="ora_results.csv",sep='\t')