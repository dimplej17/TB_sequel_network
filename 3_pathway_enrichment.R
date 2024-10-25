#######################
# Pathway Enrichment #
######################

source('1_libraries.R')
source('2_differential_gene_expression_analysis.R')

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
# write.table(gsea_results,file="gsea_results.csv",sep='\t')
# write.table(ora_results,file="ora_results.csv",sep='\t')