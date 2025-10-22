########################
# Signature Refinement #
########################

library(tidyverse)
library(GEOquery)
library(ggplot2)
library(limma)
library(ggrepel)

basedir = '/home/obaranov/projects/TBSequel/DimplesProject/'

################################################
# 1. Lower respiratory infections datasets #
################################################

gse26 <- getGEO("GSE42026")
gse44 = getGEO("GSE60244")
gse12 = getGEO("GSE40012")
# gds42026 <- getGEO(filename='/var/folders/2t/yrqd68g17sl423r_hbm32sz00000gn/T//Rtmphhl1rN/GPL6947.soft.gz')
# not working
# write_xlsx(gds42026@dataTable@table, path = "table_gse42026.xlsx")

# using bgx file i got the file "table_gse42026_controls.xlsx"


# raw_GSE42026 <- readBGX("/Users/dimplejanardhan/Downloads/LMU_Clinic_Student Assistant/GSE/GPL6947_HumanHT-12_V3_0_R1_11283641_A.bgx")



# reading the series matrix file for all 3 studies - Metadata

# meta_GSE42026 = getGEO(filename="/Users/dimplejanardhan/Downloads/LMU_Clinic_Student Assistant/GSE/Lower_respiratory_infections_datasets/GSE42026/GSE42026_series_matrix.txt.gz")
# # output: Using locally cached version of GPL6947 found here: "/var/folders/2t/yrqd68g17sl423r_hbm32sz00000gn/T//Rtmphhl1rN/GPL6947.soft.gz"

# meta_GSE60244 = getGEO(filename = "/Users/dimplejanardhan/Downloads/LMU_Clinic_Student Assistant/GSE/Lower_respiratory_infections_datasets/GSE60244/GSE60244_series_matrix.txt.gz")

# meta_GSE40012 = getGEO(filename = "/Users/dimplejanardhan/Downloads/LMU_Clinic_Student Assistant/GSE/Lower_respiratory_infections_datasets/GSE40012/GSE40012_series_matrix.txt.gz")
# # output: Using locally cached version of GPL6947 found here: "/var/folders/2t/yrqd68g17sl423r_hbm32sz00000gn/T//Rtmphhl1rN/GPL6947.soft.gz" 


# meta_GSE1234 is the ExpressionSet object so:
exprs_matrix_gse42026 <-  exprs(gse26[[1]])
exprs_matrix_gse60244 <-  exprs(gse44[[1]])
exprs_matrix_gse40012 <-  exprs(gse12[[1]])

# limma needs log transformed data 
abslog = function(x){
  sig = sign(x)
  return(sig * log2(abs(x)))
}


# Convert the matrix to a DataFrame; why actually?
exprs_df_gse42026 <- as.data.frame(exprs_matrix_gse42026 %>% apply(c(1,2),abslog), 
                        row.names = rownames(exprs_matrix_gse42026))
exprs_df_gse60244 <- as.data.frame(exprs_matrix_gse60244 %>% apply(c(1,2),abslog), 
                        row.names = rownames(exprs_matrix_gse60244))
exprs_df_gse40012 <- as.data.frame(exprs_matrix_gse40012  %>% apply(c(1,2),abslog), 
                        row.names = rownames(exprs_matrix_gse40012))
#write.csv(exprs_df_gse40012, "assaydata_gse40012.csv", row.names = TRUE)
#write.csv(exprs_df_gse42026, "assaydata_gse42026.csv", row.names = TRUE)
#write.csv(exprs_df_gse60244, "assaydata_gse60244.csv", row.names = TRUE)


# including the symbol name in the expression matrix for the corresponding IluminIDs to collapse/average the same ones


# averaging using avereps()
exprs_df_gse42026_avereps <- avereps(exprs_df_gse42026, ID = gse26[[1]]@featureData@data$Symbol)
exprs_df_gse60244_avereps <- avereps(exprs_df_gse60244, ID = gse44[[1]]@featureData@data$Symbol)
exprs_df_gse40012_avereps <- avereps(exprs_df_gse40012, ID = gse12[[1]]@featureData@data$Symbol)



########################
########################
# GSE42026 #
########################
########################

# adding a new column "Type" - levels - Control, Bacterial, H1N109, RSV
df_meta_GSE42026 <- gse26[[1]]@phenoData@data
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
saveRDS(fit.b,file=paste0(basedir,"/FilteredData/GSE42026-DE-Results.RData"))


########################
########################
# GSE60244 #
########################
########################

ill_ids_exprs <- rownames(exprs_matrix_gse60244)
ill_ids_meta <- rownames(gse44[[1]]@featureData@data)
identical(ill_ids_exprs, ill_ids_meta)

exprs_df_gse60244_avereps <- avereps(exprs_df_gse60244, ID = gse44[[1]]@featureData@data$Symbol)

df_meta_GSE60244 <- gse44[[1]]@phenoData@data
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
saveRDS(fit.b,file=paste0(basedir,"/FilteredData/GSE60244-DE-Results.RData"))


########################
########################
# GSE40012 #
########################
########################

ill_ids_exprs <- rownames(exprs_matrix_gse40012)
ill_ids_meta <- rownames(gse12[[1]]@featureData@data)
identical(ill_ids_exprs, ill_ids_meta)

exprs_df_gse40012_avereps <- avereps(exprs_df_gse40012, ID = gse12[[1]]@featureData@data$Symbol)

df_meta_GSE40012 <- gse12[[1]]@phenoData@data
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
saveRDS(fit.b,file=paste0(basedir,"/FilteredData/GSE40012-DE-Results.RData"))



############################
# 2. Sarcoidosis datasets #
###########################

#############
# GSE83456 #
############

gse56 = getGEO("GSE83456")


exprs_matrix_gse83456 <- gse56[[1]]@assayData[["exprs"]]
exprsmatrix_df_gse83456 <- as.data.frame(exprs_matrix_gse83456, row.names = rownames(exprs_matrix_gse83456))

ill_ids_exprs <- rownames(exprs_matrix_gse83456)
ill_ids_meta <- rownames(gse56[[1]]@featureData@data)
identical(ill_ids_exprs, ill_ids_meta)

exprs_df_gse83456_avereps <- avereps(exprsmatrix_df_gse83456, ID = gse56[[1]]@featureData@data$Symbol)

pheno_data_GSE83456 <- gse56[[1]]@phenoData@data
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
saveRDS(fit.b,file=paste0(basedir,"/FilteredData/GSE83456-DE-Results.RData"))

#############
# GSE42830 #
#############

gse30 = getGEO("GSE42830")

exprs_matrix_gse42830 <- gse30[[1]]@assayData[["exprs"]]
exprsmatrix_df_gse42830 <- as.data.frame(exprs_matrix_gse42830, row.names = rownames(exprs_matrix_gse42830))
# write_xlsx(exprsmatrix_df_GSE42830, path = "/Users/dimplejanardhan/Downloads/LMU_Clinic_Student Assistant/GSE/Sarcoidosis_datasets/GSE42830/assaydata_GSE42830.xlsx")


ill_ids_exprs <- rownames(exprs_matrix_gse42830)
ill_ids_meta <- rownames(gse30[[1]]@featureData@data)
identical(ill_ids_exprs, ill_ids_meta)

exprs_df_gse42830_avereps <- avereps(exprsmatrix_df_gse42830, ID = gse30[[1]]@featureData@data$Symbol)

pheno_data_GSE42830 <- gse30[[1]]@phenoData@data
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
saveRDS(fit.b,file=paste0(basedir,"/FilteredData/GSE42830-DE-Results.RData"))

#############
# GSE42826 #
#############

gse26 = getGEO("GSE42826")

exprs_matrix_gse42826 <- gse26[[1]]@assayData[["exprs"]]
exprsmatrix_df_gse42826 <- as.data.frame(exprs_matrix_gse42826, row.names = rownames(exprs_matrix_gse42826))

ill_ids_exprs <- rownames(exprs_matrix_gse42826)
ill_ids_meta <- rownames(gse26[[1]]@featureData@data)
identical(ill_ids_exprs, ill_ids_meta)

exprs_df_gse42826_avereps <- avereps(exprsmatrix_df_gse42826, ID = gse26[[1]]@featureData@data$Symbol)

pheno_data_GSE42826 <- gse26[[1]]@phenoData@data
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
saveRDS(fit.b,file=paste0(basedir,"/FilteredData/GSE42826-DE-Results.RData"))

################################################
# 3. Visualizing results by scatterplots #
################################################

getAllSig = function(degtab, tbsdegs, setlabel){
  tmptab = NULL
  for(coefpair in degtab$coefficients %>% colnames){
    gsetmp = topTable(degtab, coef = coefpair, number = Inf) %>% rownames_to_column('Assay')

    ############# !!!!!!!!!!!!!!!!!!!! go on from here
    codeg12 = merge(gsetmp, tbsdegs, by = 'Assay')
    codeg12['abslog'] = abs(codeg12['logFC'])
    codeg12['absest'] = abs(codeg12['estimate'])
    codeg12['dolabel'] = abs(codeg12['logFC']) > 0.5
    codeg12['labtext'] = ""
    codeg12[codeg12['dolabel'] %>% pull(),'labtext'] = codeg12[codeg12['dolabel'] %>% pull(),'Assay']

    p = ggplot(codeg12, aes(x = absest, y = abslog)) + 
      geom_point() +
      geom_hline(yintercept = 0.5, col = 'red') + 
      geom_text_repel(aes(label = labtext))
    ggsave(paste0(basedir,'/OutputData/',setlabel,coefpair,'.pdf'), plot = p)



    tmp = codeg12[codeg12$dolabel,]
    tmp['dataContrast'] = paste0(setlabel,'_', coefpair)
    tmptab = rbind(tmptab, tmp)
  }
  return(tmptab)
}

################################################
# 4. Displaying all scatterplots #
################################################


tbsdegs = readRDS(paste0(basedir, '/FilteredData/MWtest_results_paired.RData'))
gse12degs = readRDS(paste0(basedir, '/FilteredData/GSE40012-DE-Results.RData')) 
gse44degs = readRDS(paste0(basedir, '/FilteredData/GSE60244-DE-Results.RData')) 
gse26degs = readRDS(paste0(basedir, '/FilteredData/GSE42026-DE-Results.RData')) 

gse12sig = getAllSig(gse12degs, tbsdegs, 'GSE40012')
gse44sig = getAllSig(gse44degs, tbsdegs, 'GSE60244')
gse26sig = getAllSig(gse26degs, tbsdegs, 'GSE42026')

gsesig = rbind(gse12sig, gse44sig, gse26sig)

tbssig = tbsdegs[tbsdegs$Adjusted_pval < 0.05,] 
candidates = tbssig[!(tbssig$Assay %in% (gsesig$Assay %>% unique())),]
candidates %>% write_delim(paste0(basedir,'/OutputData/candidates.csv'),delim = '\t')


###############################################################
### statistics & exclusion of genes that are DEGs elsewhere ###
###############################################################

