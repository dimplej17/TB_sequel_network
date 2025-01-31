########################
# Signature Refinement #
########################

source('1_libraries.R')
source('2_differential_gene_expression_analysis.R')
source('4_coexpression_network_analysis.R')

################################################
# 1. Lower respiratory infections datasets #
################################################

# gse1 <- getGEO("GSE42026") --> after running this I got this thing pasted below so the GDS file exists already I'm loading the existing (compressed?) file instead of downloading the GDS file each time
# gds42026 <- getGEO(filename='/var/folders/2t/yrqd68g17sl423r_hbm32sz00000gn/T//Rtmphhl1rN/GPL6947.soft.gz')

# write_xlsx(gds42026@dataTable@table, path = "table_gse42026.xlsx")

# using bgx file i got the file "table_gse42026_controls.xlsx"


# raw_GSE42026 <- readBGX("/Users/dimplejanardhan/Downloads/LMU_Clinic_Student Assistant/GSE/GPL6947_HumanHT-12_V3_0_R1_11283641_A.bgx")



# reading the series matrix file for all 3 studies - Metadata

meta_GSE42026 = getGEO(filename="/Users/dimplejanardhan/Downloads/LMU_Clinic_Student Assistant/GSE/Lower_respiratory_infections_datasets/GSE42026/GSE42026_series_matrix.txt.gz")
# output: Using locally cached version of GPL6947 found here: "/var/folders/2t/yrqd68g17sl423r_hbm32sz00000gn/T//Rtmphhl1rN/GPL6947.soft.gz"

meta_GSE60244 = getGEO(filename = "/Users/dimplejanardhan/Downloads/LMU_Clinic_Student Assistant/GSE/Lower_respiratory_infections_datasets/GSE60244/GSE60244_series_matrix.txt.gz")

meta_GSE40012 = getGEO(filename = "/Users/dimplejanardhan/Downloads/LMU_Clinic_Student Assistant/GSE/Lower_respiratory_infections_datasets/GSE40012/GSE40012_series_matrix.txt.gz")
# output: Using locally cached version of GPL6947 found here: "/var/folders/2t/yrqd68g17sl423r_hbm32sz00000gn/T//Rtmphhl1rN/GPL6947.soft.gz" 


# meta_GSE1234 is the ExpressionSet object so:
exprs_matrix_gse42026 <- meta_GSE42026@assayData[["exprs"]]
exprs_matrix_gse60244 <- meta_GSE60244@assayData[["exprs"]]
exprs_matrix_gse40012 <- meta_GSE40012@assayData[["exprs"]]


# Convert the matrix to a DataFrame
exprs_df_gse42026 <- as.data.frame(exprs_matrix_gse42026, row.names = rownames(exprs_matrix_gse42026))
exprs_df_gse60244 <- as.data.frame(exprs_matrix_gse60244, row.names = rownames(exprs_matrix_gse60244))
exprs_df_gse40012 <- as.data.frame(exprs_matrix_gse40012, row.names = rownames(exprs_matrix_gse40012))
#write.csv(exprs_df_gse40012, "assaydata_gse40012.csv", row.names = TRUE)
#write.csv(exprs_df_gse42026, "assaydata_gse42026.csv", row.names = TRUE)
#write.csv(exprs_df_gse60244, "assaydata_gse60244.csv", row.names = TRUE)


# Ensuring the IluminaIDs are the same in both
ill_ids_exprs <- rownames(exprs_matrix_gse42026)
ill_ids_meta <- rownames(meta_GSE42026@featureData@data)

# Check if they match
identical(ill_ids_exprs, ill_ids_meta)

# including the symbol name in the expression matrix for the corresponding IluminIDs to collapse/average the same ones
#exprs_df_gse42026$Symbol <- meta_GSE42026@featureData@data$Symbol 

# Check for missing values in Symbol: If Symbol contains NA values, these could cause issues with the grouping. You can check this using:
sum(is.na(meta_GSE42026@featureData@data$Symbol))

# Check of Blank strings in Symbol
sum(meta_GSE42026@featureData@data$Symbol == "")


# Confirm Alignment of IlluminaIDs and Symbols
all(rownames(exprs_df_gse42026) == ill_ids_meta) # TRUE

# If they don’t match, you can find mismatched rows using:
# which(rownames(exprs_df_gse42026) != ill_ids_meta)


# Check for Duplicates in Symbol
# sum(duplicated(meta_GSE42026@featureData@data$Symbol)) # 23643

# To see which Symbols are duplicated:
# meta_GSE42026@featureData@data$Symbol[duplicated(meta_GSE42026@featureData@data$Symbol)]

# using collapseRows()
# exprs_df_gse42026_averaged <- collapseRows(datET = exprs_df_gse42026, rowID = rownames(exprs_df_gse42026), rowGroup = meta_GSE42026@featureData@data$Symbol, method = "Average")

# averaging using avereps()
exprs_df_gse42026_avereps <- avereps(exprs_df_gse42026, ID = meta_GSE42026@featureData@data$Symbol)



########################
########################
# GSE42026 #
########################
########################

# adding a new column "Type" - levels - Control, Bacterial, H1N109, RSV
df_meta_GSE42026 <- meta_GSE42026@phenoData@data
df_meta_GSE42026 <- as.data.frame(df_meta_GSE42026)
df_meta_GSE42026$source_name_ch1 <- as.factor(df_meta_GSE42026$source_name_ch1)
levels(df_meta_GSE42026$source_name_ch1)

df_meta_GSE42026 <- df_meta_GSE42026 %>%
  mutate(Type = case_when(
    source_name_ch1 == "Whole Blood from healthy control" ~ "Control",
    source_name_ch1 == "Whole Blood from patient with gram positive bacterial infection" ~ "Bacterial",
    source_name_ch1 == "Whole Blood from patient with H1N1/09 influenza infection" ~ "H1N109",
    source_name_ch1 == "Whole Blood from patient with RSV infection" ~ "RSV"
  ))

control_Samples = row.names(df_meta_GSE42026[df_meta_GSE42026$Type=="Control",,drop=FALSE])
disease_Samples = row.names(df_meta_GSE42026[df_meta_GSE42026$Type!="Control",,drop=FALSE])

groups = df_meta_GSE42026$Type
f <- factor(groups)
design <- model.matrix(~0+f)
colnames(design) <- levels(f)
cont.matrix <- makeContrasts(Bacterial-Control,H1N109-Control,RSV-Control,levels=design)

fit <- lmFit(exprs_df_gse42026_avereps, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit.b <- eBayes(fit2)
save(fit.b,file="GSE42026-DE-Results.RData")


########################
########################
# GSE60244 #
########################
########################

ill_ids_exprs <- rownames(exprs_matrix_gse60244)
ill_ids_meta <- rownames(meta_GSE60244@featureData@data)
identical(ill_ids_exprs, ill_ids_meta)

exprs_df_gse60244_avereps <- avereps(exprs_df_gse60244, ID = meta_GSE60244@featureData@data$Symbol)

df_meta_GSE60244 <- meta_GSE60244@phenoData@data
df_meta_GSE60244 <- as.data.frame(df_meta_GSE60244)
df_meta_GSE60244$source_name_ch1 <- as.factor(df_meta_GSE60244$source_name_ch1)
levels(df_meta_GSE60244$source_name_ch1)

df_meta_GSE60244 <- df_meta_GSE60244 %>%
  mutate(Type = case_when(
    source_name_ch1 == "whole blood, BACTERIA" ~ "Bacteria",
    source_name_ch1 == "whole blood, COINFECTION" ~ "Coinfection",
    source_name_ch1 == "whole blood, Healthy Control" ~ "Control",
    source_name_ch1 == "whole blood, VIRUS" ~ "Virus"
  ))

control_Samples = row.names(df_meta_GSE60244[df_meta_GSE60244$Type=="Control",,drop=FALSE])
disease_Samples = row.names(df_meta_GSE60244[df_meta_GSE60244$Type!="Control",,drop=FALSE])

groups = df_meta_GSE60244$Type
f <- factor(groups)
design <- model.matrix(~0+f)
colnames(design) <- levels(f)
cont.matrix <- makeContrasts(Bacteria-Control,Virus-Control,Coinfection-Control,levels=design)
fit <- lmFit(exprs_df_gse60244_avereps, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit.b <- eBayes(fit2)
save(fit.b,file="GSE60244-DE-Results.RData")


########################
########################
# GSE40012 #
########################
########################

ill_ids_exprs <- rownames(exprs_matrix_gse40012)
ill_ids_meta <- rownames(meta_GSE40012@featureData@data)
identical(ill_ids_exprs, ill_ids_meta)

exprs_df_gse40012_avereps <- avereps(exprs_df_gse40012, ID = meta_GSE40012@featureData@data$Symbol)

df_meta_GSE40012 <- meta_GSE40012@phenoData@data
df_meta_GSE40012 <- as.data.frame(df_meta_GSE40012)
df_meta_GSE40012$characteristics_ch1.1 <- as.factor(df_meta_GSE40012$characteristics_ch1.1)
levels(df_meta_GSE40012$characteristics_ch1.1)

df_meta_GSE40012 <- df_meta_GSE40012 %>%
  mutate(Type = case_when(
    characteristics_ch1.1 == "sample type: bacterial pneumonia" ~ "Bacterial",
    characteristics_ch1.1 == "sample type: healthy control" ~ "Control",
    characteristics_ch1.1 == "sample type: influenza A pneumonia" ~ "InfluenzaA",
    characteristics_ch1.1 == "sample type: mixed bacterial and influenza A pneumonia" ~ "Mixed",
    characteristics_ch1.1 == "sample type: systemic inflammatory response" ~ "SIR"
  ))

groups = df_meta_GSE40012$Type
f <- factor(groups)
design <- model.matrix(~0+f)
colnames(design) <- levels(f)
cont.matrix <- makeContrasts(Bacterial-Control,SIR-Control,InfluenzaA-Control,Mixed-Control,levels=design)
fit <- lmFit(exprs_df_gse40012_avereps, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit.b <- eBayes(fit2)
save(fit.b,file="GSE40012-DE-Results.RData")



############################
# 2. Sarcoidosis datasets #
###########################

#############
# GSE83456 #
############

meta_GSE83456 = getGEO(filename="/Users/dimplejanardhan/Downloads/LMU_Clinic_Student Assistant/GSE/Sarcoidosis_datasets/GSE83456/GSE83456_series_matrix.txt.gz")
# Op: Using locally cached version of GPL10558 found here: 
"/var/folders/2t/yrqd68g17sl423r_hbm32sz00000gn/T//Rtmphhl1rN/GPL10558.soft.gz"  

exprs_matrix_gse83456 <- meta_GSE83456@assayData[["exprs"]]
exprsmatrix_df_gse83456 <- as.data.frame(exprs_matrix_gse83456, row.names = rownames(exprs_matrix_gse83456))

ill_ids_exprs <- rownames(exprs_matrix_gse83456)
ill_ids_meta <- rownames(meta_GSE83456@featureData@data)
identical(ill_ids_exprs, ill_ids_meta)

exprs_df_gse83456_avereps <- avereps(exprsmatrix_df_gse83456, ID = meta_GSE83456@featureData@data$Symbol)


pheno_data_GSE83456 <- meta_GSE83456@phenoData@data
pheno_data_GSE83456 <- as.data.frame(pheno_data_GSE83456)
pheno_data_GSE83456$source_name_ch1 <- as.factor(pheno_data_GSE83456$source_name_ch1)
levels(pheno_data_GSE83456$source_name_ch1) 
# "Blood Human Control", "Blood Human EPTB",  "Blood Human PTB", "Blood Human Sarcoid"


# Remove rows with "Blood Human EPTB"
pheno_data_GSE83456 <- pheno_data_GSE83456[pheno_data_GSE83456$source_name_ch1 != "Blood Human EPTB", ]

# Drop unused factor levels
pheno_data_GSE83456$source_name_ch1 <- droplevels(pheno_data_GSE83456$source_name_ch1)


pheno_data_GSE83456 <- pheno_data_GSE83456 %>%
  mutate(Type = case_when(
    source_name_ch1 == "Blood Human Control" ~ "Control",
    source_name_ch1 == "Blood Human PTB" ~ "Active",
    source_name_ch1 == "Blood Human Sarcoid" ~ "Other"
  ))

control_Samples = row.names(pheno_data_GSE83456[pheno_data_GSE83456$Type=="Control",,drop=FALSE])
disease_Samples = row.names(pheno_data_GSE83456[pheno_data_GSE83456$Type=="Other",,drop=FALSE])

exprs_df_gse83456_avereps = exprs_df_gse83456_avereps[,c(control_Samples,disease_Samples)]

pheno_data_GSE83456 = pheno_data_GSE83456[c(control_Samples,disease_Samples),,drop=FALSE]

groups = pheno_data_GSE83456$Type
f <- factor(groups)
design <- model.matrix(~0+f)
colnames(design) <- levels(f)
cont.matrix <- makeContrasts(Other-Control,levels=design)
fit <- lmFit(exprs_df_gse83456_avereps, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit.b <- eBayes(fit2)
save(fit.b,file="GSE83456-DE-Results.RData")

#############
# GSE42830 #
#############

meta_GSE42830 = getGEO(filename = "/Users/dimplejanardhan/Downloads/LMU_Clinic_Student Assistant/GSE/Sarcoidosis_datasets/GSE42830/GSE42830_series_matrix.txt.gz")
# Op: Using locally cached version of GPL10558 found here: "/var/folders/2t/yrqd68g17sl423r_hbm32sz00000gn/T//Rtmphhl1rN/GPL10558.soft.gz"

exprs_matrix_gse42830 <- meta_GSE42830@assayData[["exprs"]]
exprsmatrix_df_gse42830 <- as.data.frame(exprs_matrix_gse42830, row.names = rownames(exprs_matrix_gse42830))
# write_xlsx(exprsmatrix_df_GSE42830, path = "/Users/dimplejanardhan/Downloads/LMU_Clinic_Student Assistant/GSE/Sarcoidosis_datasets/GSE42830/assaydata_GSE42830.xlsx")


ill_ids_exprs <- rownames(exprs_matrix_gse42830)
ill_ids_meta <- rownames(meta_GSE42830@featureData@data)
identical(ill_ids_exprs, ill_ids_meta)

exprs_df_gse42830_avereps <- avereps(exprsmatrix_df_gse42830, ID = meta_GSE42830@featureData@data$Symbol)

pheno_data_GSE42830 <- meta_GSE42830@phenoData@data
pheno_data_GSE42830 <- as.data.frame(pheno_data_GSE42830)
pheno_data_GSE42830$characteristics_ch1.2 <- as.factor(pheno_data_GSE42830$characteristics_ch1.2)
levels(pheno_data_GSE42830$characteristics_ch1.2) 
# [1] "disease state: Active Sarcoid"         "disease state: Baseline"              
# [3] "disease state: Control"                "disease state: lung cancer"           
# [5] "disease state: Non-active sarcoidosis" "disease state: TB" 
# Baseline -> Pneumonia sample

# Remove rows with "disease state: Baseline" & "disease state: lung cancer" & "disease state: Non-active sarcoidosis"    
pheno_data_GSE42830 <- pheno_data_GSE42830[pheno_data_GSE42830$characteristics_ch1.2 != "disease state: Baseline", ]
pheno_data_GSE42830 <- pheno_data_GSE42830[pheno_data_GSE42830$characteristics_ch1.2 != "disease state: lung cancer", ]
pheno_data_GSE42830 <- pheno_data_GSE42830[pheno_data_GSE42830$characteristics_ch1.2 != "disease state: Non-active sarcoidosis", ]

# Drop unused factor levels
pheno_data_GSE42830$characteristics_ch1.2 <- droplevels(pheno_data_GSE42830$characteristics_ch1.2)

pheno_data_GSE42830 <- pheno_data_GSE42830 %>%
  mutate(Type = case_when(
    characteristics_ch1.2 == "disease state: Control" ~ "Control",
    characteristics_ch1.2 == "disease state: TB" ~ "Active",
    characteristics_ch1.2 == "disease state: Active Sarcoid" ~ "Other"  ))

control_Samples = row.names(pheno_data_GSE42830[pheno_data_GSE42830$Type=="Control",,drop=FALSE])
disease_Samples = row.names(pheno_data_GSE42830[pheno_data_GSE42830$Type=="Other",,drop=FALSE])

exprs_df_gse42830_avereps = exprs_df_gse42830_avereps[,c(control_Samples,disease_Samples)]

pheno_data_GSE42830 = pheno_data_GSE42830[c(control_Samples,disease_Samples),,drop=FALSE]

groups = pheno_data_GSE42830$Type
f <- factor(groups)
design <- model.matrix(~0+f)
colnames(design) <- levels(f)
cont.matrix <- makeContrasts(Other-Control,levels=design)
fit <- lmFit(exprs_df_gse42830_avereps, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit.b <- eBayes(fit2)
save(fit.b,file="GSE42830-DE-Results.RData")

#############
# GSE42826 #
#############

meta_GSE42826 = getGEO(filename = "/Users/dimplejanardhan/Downloads/LMU_Clinic_Student Assistant/GSE/Sarcoidosis_datasets/GSE42826/GSE42826_series_matrix.txt.gz")
# Op: Using locally cached version of GPL10558 found here: "/var/folders/2t/yrqd68g17sl423r_hbm32sz00000gn/T//Rtmphhl1rN/GPL10558.soft.gz" 

exprs_matrix_gse42826 <- meta_GSE42826@assayData[["exprs"]]
exprsmatrix_df_gse42826 <- as.data.frame(exprs_matrix_gse42826, row.names = rownames(exprs_matrix_gse42826))

ill_ids_exprs <- rownames(exprs_matrix_gse42826)
ill_ids_meta <- rownames(meta_GSE42826@featureData@data)
identical(ill_ids_exprs, ill_ids_meta)

exprs_df_gse42826_avereps <- avereps(exprsmatrix_df_gse42826, ID = meta_GSE42826@featureData@data$Symbol)

pheno_data_GSE42826 <- meta_GSE42826@phenoData@data
pheno_data_GSE42826 <- as.data.frame(pheno_data_GSE42826)
pheno_data_GSE42826$characteristics_ch1.2 <- as.factor(pheno_data_GSE42826$characteristics_ch1.2)
levels(pheno_data_GSE42826$characteristics_ch1.2)
# [1] "disease state: Active Sarcoid"         "disease state: Control"               
# [3] "disease state: lung cancer"            "disease state: Non-active sarcoidosis"
# [5] "disease state: Pneumonia"              "disease state: TB" 

# Remove rows with "disease state: Pneumonia"  & "disease state: lung cancer" & "disease state: Non-active sarcoidosis"   
pheno_data_GSE42826 <- pheno_data_GSE42826[pheno_data_GSE42826$characteristics_ch1.2 != "disease state: Pneumonia", ]
pheno_data_GSE42826 <- pheno_data_GSE42826[pheno_data_GSE42826$characteristics_ch1.2 != "disease state: lung cancer", ]
pheno_data_GSE42826 <- pheno_data_GSE42826[pheno_data_GSE42826$characteristics_ch1.2 != "disease state: Non-active sarcoidosis", ]

# Drop unused factor levels
pheno_data_GSE42826$characteristics_ch1.2 <- droplevels(pheno_data_GSE42826$characteristics_ch1.2)

pheno_data_GSE42826 <- pheno_data_GSE42826 %>%
  mutate(Type = case_when(
    characteristics_ch1.2 == "disease state: Control" ~ "Control",
    characteristics_ch1.2 == "disease state: TB" ~ "Active",
    characteristics_ch1.2 == "disease state: Active Sarcoid" ~ "Other"  ))

control_Samples = row.names(pheno_data_GSE42826[pheno_data_GSE42826$Type=="Control",,drop=FALSE])
disease_Samples = row.names(pheno_data_GSE42826[pheno_data_GSE42826$Type=="Other",,drop=FALSE])

exprs_df_gse42826_avereps = exprs_df_gse42826_avereps[,c(control_Samples,disease_Samples)]

pheno_data_GSE42826 = pheno_data_GSE42826[c(control_Samples,disease_Samples),,drop=FALSE]

groups = pheno_data_GSE42826$Type
f <- factor(groups)
design <- model.matrix(~0+f)
colnames(design) <- levels(f)
cont.matrix <- makeContrasts(Other-Control,levels=design)
fit <- lmFit(exprs_df_gse42826_avereps, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit.b <- eBayes(fit2)
save(fit.b,file="GSE42826-DE-Results.RData")

################################################
# 3. Visualizing results by scatterplots #
################################################

DEPlot_func <- function(inputData, xaxis, yaxis,xpvalue, ypvalue,xlab,ylab) {
  inputData$TBsig <- ifelse(inputData[,xpvalue] <= 0.05&abs(inputData[,xaxis]) > 1, inputData$TBsig <- "Sig", inputData$TBsig <- "nSig")
  inputData$Othersig <- ifelse(inputData[,ypvalue] <= 0.05&abs(inputData[,yaxis]) > 1, inputData$Othersig <- "Sig", inputData$Othersig <- "nSig")
  ggplot(inputData, aes(x=.data[[xaxis]],y=.data[[yaxis]]))+
    geom_point(aes(color = Othersig,size=2)) +
    scale_x_continuous(limits=c(0,3)) +
    scale_y_continuous(limits=c(-0.5,2.5)) + 
    coord_fixed(ratio=1)  +
    scale_color_manual(values = c("grey", "red")) +
    theme_bw(base_size = 12) +
    geom_text_repel(
      data = subset(inputData, inputData[,ypvalue] <= 1 & abs(inputData[,yaxis]) > 0), 
      aes(label = Protein),
      size = 6,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(2, "lines"),
      segment.colour = "transparent") +
    labs(y=paste0("Log2( ",ylab," )"), x=paste0("Log2( ",xlab," )")) +
    theme(legend.position = "none",axis.title = element_blank(),axis.text = element_text(size=16))
}


################################################
# 4. Displaying all scatterplots #
################################################

# our olink data - UGHHHHHH - i mean our network proteins & their logFC & pvalues

network_correlations_with_symbols_andUniprot <- read.csv("/Users/dimplejanardhan/Downloads/LMU_Clinic_Student Assistant/R_Studio_LMU_Student Assistant_TB/network_correlations_with_symbols_andUniprot.csv")

listofUniProts_networkproteins <- unique(c(network_correlations_with_symbols_andUniprot$Protein1, network_correlations_with_symbols_andUniprot$Protein2))
listofUniProts_networkproteins <- as.data.frame(listofUniProts_networkproteins)
colnames(listofUniProts_networkproteins) <- "UniProt" 

replacement_symbols <- read_xlsx("/Users/dimplejanardhan/Downloads/LMU_Clinic_Student Assistant/R_Studio_LMU_Student Assistant_TB/replacedsymbolswiththenamessentbydrolga_MWtestresults.xlsx") # this has 126 proteins; there are only 43 proteins in our protein network

# removing the other 83 proteins now
replacement_symbols <- replacement_symbols %>%
  filter(UniProt %in% listofUniProts_networkproteins$UniProt)

load("~/Downloads/LMU_Clinic_Student Assistant/R_Studio_LMU_Student Assistant_TB/modules.RData")

# changing the names of proteins from UniProt ID to their respective symbols in modules.Rdata

# Convert mapping to a named vector for easy lookup
id_to_symbol <- setNames(replacement_symbols$New_Assay, replacement_symbols$UniProt)

# Function to replace UniProt IDs with Gene Symbols
replace_uniprot_with_symbol <- function(protein_list) {
  symbols <- id_to_symbol[protein_list]  # Lookup replacement
  symbols[is.na(symbols)] <- protein_list[is.na(symbols)]  # Keep original if no match
  return(symbols)
}

# Apply the replacement to all modules
modules <- lapply(modules, replace_uniprot_with_symbol)

# Save the updated modules list
save(modules, file = "modules_updated_symbols.RData")


load("/Users/dimplejanardhan/Downloads/LMU_Clinic_Student Assistant/R_Studio_LMU_Student Assistant_TB/modules_updated_symbols.RData")
load("~/Downloads/LMU_Clinic_Student Assistant/R_Studio_LMU_Student Assistant_TB/MWtest_results_unpaired.RData")

MWtest_results <- MWtest_results %>%
  filter(UniProt %in% listofUniProts_networkproteins$UniProt)

MWtest_results_updatedsymbols <- MWtest_results %>%
  left_join(replacement_symbols, by = "UniProt") %>% 
  mutate(Assay = ifelse(!is.na(New_Assay), New_Assay, Assay)) %>% # Replace only when New_Assay is not NA
  select(-New_Assay) 



MWtest_results_updatedsymbols <- as.data.frame(MWtest_results_updatedsymbols)

# Set the "Assay" column as row names
rownames(MWtest_results_updatedsymbols) <- MWtest_results_updatedsymbols$Assay

# Remove 'Assay' column
MWtest_results_updatedsymbols$Assay <- NULL


logFC <- as.data.frame(MWtest_results_updatedsymbols$estimate)
pValue <- as.data.frame(MWtest_results_updatedsymbols$Adjusted_pval)
DE_results <- cbind(logFC=logFC,pValue=pValue)
rownames(DE_results) <- rownames(MWtest_results_updatedsymbols)

# renaming the colnames because currently it's: "MWtest_results_updatedsymbols$estimate"; "MWtest_results_updatedsymbols$Adjusted_pval"
colnames(DE_results) <- c("logFC.Active.Cured", "pValue.Active.Cured")


# other studies data

load("~/Downloads/LMU_Clinic_Student Assistant/R_Studio_LMU_Student Assistant_TB/GSE42026-DE-Results.RData")
logFC <- as.data.frame(fit.b$coefficients)
pValue <- as.data.frame(fit.b$p.value)
new <- cbind(logFC=logFC,pValue=pValue)
new <- new[match(row.names(DE_results),row.names(new)),]
DE_results <- cbind(DE_results,new)


load("~/Downloads/LMU_Clinic_Student Assistant/R_Studio_LMU_Student Assistant_TB/GSE60244-DE-Results.RData")
logFC <- as.data.frame(fit.b$coefficients)
pValue <- as.data.frame(fit.b$p.value)
new <- cbind(logFC=logFC,pValue=pValue)
new <- new[match(row.names(DE_results),row.names(new)),]
DE_results <- cbind(DE_results,new)

load("~/Downloads/LMU_Clinic_Student Assistant/R_Studio_LMU_Student Assistant_TB/GSE40012-DE-Results.RData")
logFC <- as.data.frame(fit.b$coefficients)
pValue <- as.data.frame(fit.b$p.value)
new <- cbind(logFC=logFC,pValue=pValue)
new <- new[match(row.names(DE_results),row.names(new)),]
DE_results <- cbind(DE_results,new)


load("~/Downloads/LMU_Clinic_Student Assistant/R_Studio_LMU_Student Assistant_TB/GSE83456-DE-Results.RData")
logFC <- as.data.frame(fit.b$coefficients)
pValue <- as.data.frame(fit.b$p.value)
new <- cbind(logFC=logFC,pValue=pValue)
new <- new[match(row.names(DE_results),row.names(new)),]
names(new)=c(paste0("logFC.","GSE83456"),paste0("pValue.","GSE83456")) 
DE_results <- cbind(DE_results,new)


load("~/Downloads/LMU_Clinic_Student Assistant/R_Studio_LMU_Student Assistant_TB/GSE42830-DE-Results.RData")
logFC <- as.data.frame(fit.b$coefficients)
pValue <- as.data.frame(fit.b$p.value)
new <- cbind(logFC=logFC,pValue=pValue)
new <- new[match(row.names(DE_results),row.names(new)),]
names(new)=c(paste0("logFC.","GSE42830"),paste0("pValue.","GSE42830")) 
DE_results <- cbind(DE_results,new)


load("~/Downloads/LMU_Clinic_Student Assistant/R_Studio_LMU_Student Assistant_TB/GSE42826-DE-Results.RData")
logFC <- as.data.frame(fit.b$coefficients)
pValue <- as.data.frame(fit.b$p.value)
new <- cbind(logFC=logFC,pValue=pValue)
new <- new[match(row.names(DE_results),row.names(new)),]
names(new)=c(paste0("logFC.","GSE42826"),paste0("pValue.","GSE42826")) 
DE_results <- cbind(DE_results,new)



DE_results$Protein <- row.names(DE_results)
names(DE_results) <- lapply(names(DE_results), function(x) gsub(" - ",".",x))

colnames(DE_results) # cross-check correct column names
names(DE_results)[15]="logFC.Bacterial2.Control"
names(DE_results)[19]="pValue.Bacterial2.Control"
# Also, i don't see these 2^ in the plots in their code - it's ok - you make a 13th plot & use it anyway



pdf(file="/Users/dimplejanardhan/Downloads/LMU_Clinic_Student Assistant/R_Studio_LMU_Student Assistant_TB/Figure_Signature_Refinement.pdf",width = 40,height = 60)
p1=DEPlot_func(DE_results,"logFC.Active.Cured","logFC.Bacterial.Control","pValue.Active.Cured","pValue.Bacterial.Control","Active vs. Cured","Bacterial vs. Control")
p2=DEPlot_func(DE_results,"logFC.Active.Cured","logFC.H1N109.Control","pValue.Active.Cured","pValue.H1N109.Control","Active vs. Cured","H1N109 vs. Control")
p3=DEPlot_func(DE_results,"logFC.Active.Cured","logFC.RSV.Control","pValue.Active.Cured","pValue.RSV.Control","Active vs. Cured","RSV vs. Control")
p4=DEPlot_func(DE_results,"logFC.Active.Cured","logFC.Bacterial.Control","pValue.Active.Cured","pValue.Bacterial.Control","Active vs. Cured","Bacterial vs. Control")
p5=DEPlot_func(DE_results,"logFC.Active.Cured","logFC.InfluenzaA.Control","pValue.Active.Cured","pValue.InfluenzaA.Control","Active vs. Cured","InfluenzaA vs. Control")
p6=DEPlot_func(DE_results,"logFC.Active.Cured","logFC.SIR.Control","pValue.Active.Cured","pValue.SIR.Control","Active vs. Cured","SIR vs. Control")
p7=DEPlot_func(DE_results,"logFC.Active.Cured","logFC.Bacteria.Control","pValue.Active.Cured","pValue.Bacteria.Control","Active vs. Cured","Bacteria vs. Control")
p8=DEPlot_func(DE_results,"logFC.Active.Cured","logFC.Virus.Control","pValue.Active.Cured","pValue.Virus.Control","Active vs. Cured","Virus vs. Control")
p9=DEPlot_func(DE_results,"logFC.Active.Cured","logFC.Coinfection.Control","pValue.Active.Cured","pValue.Coinfection.Control","Active vs. Cured","Coinfection vs. Control")
p10=DEPlot_func(DE_results,"logFC.Active.Cured","logFC.Bacterial2.Control","pValue.Active.Cured","pValue.Bacterial2.Control","Active vs. Cured","Bacterial2 vs. Control")
p11=DEPlot_func(DE_results,"logFC.Active.Cured","logFC.GSE42826","pValue.Active.Cured","pValue.GSE42826","Active vs. Cured","Sarcidosis vs. Control")
p12=DEPlot_func(DE_results,"logFC.Active.Cured","logFC.GSE42830","pValue.Active.Cured","pValue.GSE42830","Active vs. Cured","Sarcoidosis vs. Control")
p13=DEPlot_func(DE_results,"logFC.Active.Cured","logFC.GSE83456","pValue.Active.Cured","pValue.GSE83456","Active vs. Cured","Sarcoidosis vs. Control")


cowplot::plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,ncol=3,nrow=5)
dev.off()
