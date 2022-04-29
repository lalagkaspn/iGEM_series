## This script is used to download and pre-process series (GSE) from GEO for 
## PDAC patients (normal and tumor tissues) with samples taken
## from blood (whole blood samples or circulating tumor samples)

library(dplyr)
library(GEOquery)
library(org.Hs.eg.db)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(tidyr)
library(limma)
library(openxlsx)
library(EnhancedVolcano)

##### Downloading data #####
datasets = c("GSE125158", "GSE49641", "GSE74629", "GSE18670")
# the latter is the dataset for GSE18670 and contains the full gene expression matrix

# Run this before getGEO
# Probably, there is a bug with the newest version of readr and GEOquery
# https://github.com/seandavi/GEOquery/issues/114
readr::local_edition(1)

# Download data
GEOsets = list()
for (i in 1:length(datasets)){
  GEOsets[[i]] = getGEO(datasets[i])
}; rm(i)
GEOsets = unlist(GEOsets)
names(GEOsets) = datasets; rm(datasets)

###### pData #####
pdata = list()
for(i in 1:length(GEOsets)) {
  pdata[[i]] = pData(GEOsets[[i]])
}
names(pdata) = names(GEOsets); rm(i)

## Keep only necessary pdata
filt_pdata = list()

# GSE125158 
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125158
# Paper: https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC6447845/
# Comments:
#   - 30 samples (30 patients)
#     - 17 PDAC whole blood samples
#     - 13 normal whole blood samples

# Select informative data only
filt_pdata[["GSE125158"]] = pdata$GSE125158 %>%
  dplyr::select(GEO_accession = geo_accession,
         Patient_ID = title,
         Platform = platform_id,
         Gender = "gender:ch1",
         Age = "age:ch1",
         Tissue_type = "diagnosis:ch1")

# Change Tissue_type to tumor and non_tumor
for (i in 1:length(filt_pdata$GSE125158$Tissue_type)) {
  if (filt_pdata$GSE125158$Tissue_type[i] == "healthy") {
    filt_pdata$GSE125158$Tissue_type[i] = "non_tumor"
  } else {
    filt_pdata$GSE125158$Tissue_type[i] = "tumor"
  }
}; rm(i)

# Clear Patient_ID names
filt_pdata$GSE125158$Patient_ID = sub("\\:.*", "", filt_pdata$GSE125158$Patient_ID)

# Transform to factors with consistent universal levels
filt_pdata$GSE125158$Tissue_type = factor(x = filt_pdata$GSE125158$Tissue_type,
                                          levels = c("non_tumor","tumor"),
                                          labels = c("non_tumor","tumor"))
filt_pdata$GSE125158$Gender = factor(x = filt_pdata$GSE125158$Gender,
                                     levels = c("female","male"),
                                     labels = c("female","male"))
filt_pdata$GSE125158$Age = as.numeric(filt_pdata$GSE125158$Age)

# GSE49641
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49641
# Paper: https://pubmed.ncbi.nlm.nih.gov/25069573/
# Comments:
#   - 36 samples (36 patients)
#     - 18 unresectable PDAC
#     - 18 normal

# Select informative data only
filt_pdata[["GSE49641"]] = pdata$GSE49641 %>%
  dplyr::select(GEO_accession = geo_accession,
         Patient_ID = title,
         Platform = platform_id,
         Gender = "gender:ch1",
         Age = "age:ch1",
         Tissue_type = "disease state:ch1",
         AJCC_classification = "pathological staging:ch1")

# Change Tissue_type to tumor and non_tumor
for (i in 1:length(filt_pdata$GSE49641$Tissue_type)) {
  if (filt_pdata$GSE49641$Tissue_type[i] == "control") {
    filt_pdata$GSE49641$Tissue_type[i] = "non_tumor"
  } else {
    filt_pdata$GSE49641$Tissue_type[i] = "tumor"
  }
}; rm(i)

# Transform to factors with consistent universal levels
filt_pdata$GSE49641$Tissue_type = factor(x = filt_pdata$GSE49641$Tissue_type, 
                                         levels = c("non_tumor","tumor"), 
                                         labels = c("non_tumor","tumor"))
filt_pdata$GSE49641$Gender = factor(x = filt_pdata$GSE49641$Gender, 
                                    levels = c("Female","Male"), 
                                    labels = c("female","male"))
filt_pdata$GSE49641$Age = as.numeric(filt_pdata$GSE49641$Age)
filt_pdata$GSE49641$AJCC_classification = gsub("IV", "4", 
                                               filt_pdata$GSE49641$AJCC_classification)
filt_pdata$GSE49641$AJCC_classification = gsub("III", "3", 
                                               filt_pdata$GSE49641$AJCC_classification)
filt_pdata$GSE49641$AJCC_classification = gsub("NA", NA, 
                                               filt_pdata$GSE49641$AJCC_classification)
filt_pdata$GSE49641$AJCC_classification = factor(x = filt_pdata$GSE49641$AJCC_classification,
                                                 levels = c("1a","1b","2a","2b","3","4"),
                                                 labels = c("1a","1b","2a","2b","3","4"))

# GSE74629
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74629
# Paper: https://pubmed.ncbi.nlm.nih.gov/29617451/
# Comments:
#   - 50 samples (50 patients)
#     - 36 PDAC
#     - 14 normal

# Select informative data only
filt_pdata[["GSE74629"]] = pdata$GSE74629 %>%
  dplyr::select(GEO_accession = geo_accession,
         Patient_ID = title,
         Platform = platform_id,
         Gender = "gender:ch1",
         Age = "age:ch1",
         Tissue_type = "diagnosis:ch1",
         T_stage = "pathological staging:ch1",
         Chronic_pancreatitis = "chronic pancreatitis:ch1",
         Type_II_diabetes = "diabetes:ch1")
# Patient_ID are defined similarly to GSE49641

# Change Tissue_type to tumor and non_tumor
for (i in 1:length(filt_pdata$GSE74629$Tissue_type)) {
  if (filt_pdata$GSE74629$Tissue_type[i] == "healthy") {
    filt_pdata$GSE74629$Tissue_type[i] = "non_tumor"
  } else {
    filt_pdata$GSE74629$Tissue_type[i] = "tumor"
  }
}; rm(i)

# Transform to factors with consistent universal levels
filt_pdata$GSE74629$Tissue_type = factor(x = filt_pdata$GSE74629$Tissue_type,
                                         levels = c("non_tumor","tumor"),
                                         labels = c("non_tumor","tumor"))
filt_pdata$GSE74629$Gender = factor(x = filt_pdata$GSE74629$Gender,
                                    levels = c("female","male"),
                                    labels = c("female","male"))
filt_pdata$GSE74629$Age = as.numeric(filt_pdata$GSE74629$Age)
filt_pdata$GSE74629$Chronic_pancreatitis = as.factor(filt_pdata$GSE74629$Chronic_pancreatitis)
filt_pdata$GSE74629$Type_II_diabetes = as.factor(filt_pdata$GSE74629$Type_II_diabetes)
filt_pdata$GSE74629$T_stage = factor(x = filt_pdata$GSE74629$T_stage,
                                     levels = c("1","2","3","4"),
                                     labels = c("1","2","3","4"))

# GSE18670 (GDS4329)
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18670
# Paper: https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC3599097/
# Comments:
#   - 24 samples (6 patients)
#       - 6 circulating tumor (CTC)
#       - 6 haematological (G)
#       - 6 original tumor (T)
#       - 6 non-tumor pancreatic control (P)

# Select pdata
filt_pdata[["GSE18670"]] = pdata$GSE18670 %>%
  dplyr::select(GEO_accession = geo_accession,
         Patient_ID = title,
         Platform = platform_id,
         Tissue_type = "source_name_ch1",
         Age = "age:ch1",
         AJCC_classification = "ajcc classif. (2002):ch1",
         Gender = "gender:ch1",
         N_stage = "n-stage:ch1",
         T_stage = "t-stage:ch1",
         Tumor_grade = "tumor grade:ch1")

# Keep only circulating tumor (CTC)
# We do not keep the control samples here because they derive from pancreatic tissue
patients_keep = grep("_CTC_",x = filt_pdata$GSE18670$Tissue_type)

filt_pdata$GSE18670 = filt_pdata$GSE18670[patients_keep, ] ; rm(patients_keep)

# Add M_stage column
filt_pdata$GSE18670$M_stage = "0"

# Replace Tissue_type with tumor and non_tumor values
tumor = which(grepl("_CTC_", filt_pdata$GSE18670$Tissue_type))
for (i in tumor) {
  filt_pdata$GSE18670[i, "Tissue_type"] = "tumor"
}

non_tumor = which(grepl("_P_", filt_pdata$GSE18670$Tissue_type))
for (i in non_tumor) {
  filt_pdata$GSE18670[i, "Tissue_type"] = "non_tumor"
}
rm(tumor, non_tumor, i)

# Reorder columns
filt_pdata$GSE18670 = filt_pdata$GSE18670 %>%
  dplyr::select(GEO_accession, Patient_ID, Platform, Gender, Age, Tissue_type, 
         AJCC_classification, T_stage, N_stage, M_stage, Tumor_grade)

# Transform to factors with consistent universal levels
filt_pdata$GSE18670$Gender = factor(x = filt_pdata$GSE18670$Gender, 
                                    levels = c("Female", "Male"), 
                                    labels = c("female", "male"))
filt_pdata$GSE18670$Age = as.numeric(filt_pdata$GSE18670$Age)
filt_pdata$GSE18670$Tissue_type = factor(x = filt_pdata$GSE18670$Tissue_type, 
                                         levels = c("non_tumor", "tumor"), 
                                         labels = c("non_tumor", "tumor"))
filt_pdata$GSE18670$AJCC_classification = factor(x = filt_pdata$GSE18670$AJCC_classification, 
                                                 levels = c("1a","1b","2a","2b","3","4"),
                                                 labels = c("1a","1b","2a","2b","3","4"))
filt_pdata$GSE18670$T_stage = factor(x = filt_pdata$GSE18670$T_stage, 
                                     levels = c("1","2","3","4"), 
                                     labels = c("1","2","3","4"))
filt_pdata$GSE18670$N_stage = factor(x = filt_pdata$GSE18670$N_stage, 
                                     levels = c("0","1","2"), 
                                     labels = c("0","1","2"))
filt_pdata$GSE18670$M_stage = factor(x = filt_pdata$GSE18670$M_stage, 
                                     levels = c("0", "1"), 
                                     labels = c("non_metastatic", "metastatic"))
filt_pdata$GSE18670$Tumor_grade = factor(x = filt_pdata$GSE18670$Tumor_grade, 
                                         levels = c("1","2","3"), labels = c("1","2","3"))

# Add NAs to non_tumor samples clinical data
for(i in 1:length(filt_pdata$GSE18670$Tissue_type)){
  if(filt_pdata$GSE18670$Tissue_type[i] == "non_tumor") {
    filt_pdata$GSE18670$AJCC_classification[i] = NA
    filt_pdata$GSE18670$T_stage[i] = NA
    filt_pdata$GSE18670$N_stage[i] = NA
    filt_pdata$GSE18670$M_stage[i] = NA
    filt_pdata$GSE18670$Tumor_grade[i] = NA
  }
}; rm(i)

# full_pdata
# Keep only information for Study, GEO_accession and Tissue_type
# Useful for QC analysis
pdata125158 = filt_pdata$GSE125158 %>%
  dplyr::select(GEO_accession, Tissue_type) %>%
  dplyr::mutate(Study = "GSE125158")
pdata49641 = filt_pdata$GSE49641 %>%
  dplyr::select(GEO_accession, Tissue_type) %>%
  dplyr::mutate(Study = "GSE49641")
pdata74629 = filt_pdata$GSE74629 %>%
  dplyr::select(GEO_accession, Tissue_type) %>%
  dplyr::mutate(Study = "GSE74629")
pdata18670 = filt_pdata$GSE18670 %>%
  dplyr::select(GEO_accession, Tissue_type) %>%
  dplyr::mutate(Study = "GSE18670")

full_pdata = rbind(pdata125158, pdata49641, pdata74629, pdata18670)
rownames(full_pdata) = full_pdata$GSM
rm(pdata125158, pdata49641, pdata74629, pdata18670)

##### Expression data #####
esets = list()
for (i in 1:length(GEOsets)) {
  esets[[i]] = exprs(GEOsets[[i]])
  esets[[i]] = as.data.frame(esets[[i]])
}; rm(i)
names(esets) = names(GEOsets)

## Calculate NAs values
na_esets = c()
for(i in 1:length(esets)) {
  na_esets[i] = sum(is.na(esets[[i]]))
}; rm(i)
names(na_esets) = names(GEOsets)
na_esets # no missing values

##### Annotation esets with Entrez ID's #####
# GSE125158
fdata125158 = fData(GEOsets$GSE125158) %>%
  dplyr::select(ID, ENTREZ_GENE_ID = GENE) %>%
  dplyr::filter(!grepl("///", ENTREZ_GENE_ID)) %>%
  dplyr::filter(nchar(ENTREZ_GENE_ID)>0)
esets[["GSE125158"]]$ID = rownames(esets[["GSE125158"]])
esets[["GSE125158"]] = esets[["GSE125158"]] %>%
  inner_join(fdata125158) %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  summarise_all(mean, na.rm = TRUE)
rm(fdata125158)

# GSE49641
# Convert GB_LIST (RefSeq) to ENTREZ IDs
ref = org.Hs.egREFSEQ2EG
mapped_genes_official = mappedkeys(ref)
ref_df = as.data.frame(ref[mapped_genes_official])
ref_df = ref_df %>% dplyr::rename(ENTREZ_GENE_ID = gene_id, RefSeq = accession)
length(unique(ref_df$RefSeq)) == length(ref_df$RefSeq) # TRUE --> No duplicates in RefSeq

# The GB_LIST column contains RefSeq ID's which we aim to map to Entrez ID's
# It might be the case that a probe is matched to multiple RefSeq ID's which,
# however, correspond to the same Entrez. In order to account for that and keep
# these probes while discarding probes which match to multiple Entrez ID's we do
# the following

fdata49641 = fData(GEOsets$GSE49641) %>%
  dplyr::select(ID, RefSeq = GB_LIST) %>%
  dplyr::filter(is.na(RefSeq)==F) %>%
  dplyr::filter(nchar(RefSeq)>0)

names = c()
for (i in 1:202){
  names = c(names, paste0("RS", i))
}
fdata49641$ID = as.character(fdata49641$ID)
fdata49641_sep = melt(separate(fdata49641, col = RefSeq, into = names, sep = ","),
                      id.vars = "ID") %>%
  dplyr::select(-variable) %>%
  dplyr::rename(RefSeq = value) %>%
  dplyr::filter(is.na(RefSeq)==F) %>%
  dplyr::filter(nchar(RefSeq)>0) %>%
  inner_join(ref_df, by = "RefSeq") %>%
  distinct() %>%
  group_by(ID)
fdata49641_sep$Filter = NA

# Keep only the probes which match to a unique Entrez, even if these mapped to
# multiple RefSeq ID's
for (i in 1:nrow(fdata49641_sep)){
  if (length(unique(fdata49641_sep$ENTREZ_GENE_ID[fdata49641_sep$ID==fdata49641_sep$ID[i]]))==1){
    fdata49641_sep$Filter[i] = "keep"
  } else {
    fdata49641_sep$Filter[i] = "discard"
  }
}

final_fdata49641 = fdata49641_sep %>%
  dplyr::filter(Filter == "keep") %>%
  distinct() %>% 
  dplyr::select(-RefSeq, -Filter)
esets[["GSE49641"]]$ID = as.character(rownames(esets[["GSE49641"]]))
esets[["GSE49641"]] = esets[["GSE49641"]] %>%
  inner_join(final_fdata49641, by = "ID") %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  summarise_all(mean, na.rm = TRUE)
rm(fdata49641, fdata49641_sep, final_fdata49641, names)

# GSE74629
fdata74629 = fData(GEOsets$GSE74629) %>%
  dplyr::select(ID, ENTREZ_GENE_ID = Entrez_Gene_ID) %>%
  dplyr::filter(!grepl("///", ENTREZ_GENE_ID)) %>%
  dplyr::filter(nchar(ENTREZ_GENE_ID)>0)
esets[["GSE74629"]]$ID = rownames(esets[["GSE74629"]])
esets[["GSE74629"]] = esets[["GSE74629"]] %>%
  inner_join(fdata74629) %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  summarise_all(mean, na.rm = TRUE)
rm(fdata74629)

# GSE18670
fdata18670 = fData(GEOsets$GSE18670) %>%
  dplyr::select(ID, ENTREZ_GENE_ID) %>%
  dplyr::filter(!grepl("///", ENTREZ_GENE_ID)) %>%
  dplyr::filter(nchar(ENTREZ_GENE_ID)>0)
esets[["GSE18670"]]$ID = rownames(esets[["GSE18670"]])
esets[["GSE18670"]] = esets[["GSE18670"]] %>%
  inner_join(fdata18670) %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  summarise_all(mean, na.rm = TRUE)
rm(fdata18670)

##### z-score-transformation #####
# KBZ transformation method ( https:://www.biostars.org/p/283083/ )
z = list()
for(i in 1:length(esets)){
  df = as.data.frame(esets[[i]]) %>%
    dplyr::select(-ENTREZ_GENE_ID)
  t = as.data.frame(t(df))
  z_t = sapply(t, function(t) (t-mean(t, na.rm = T))/sd(t, na.rm = T))
  z[[i]] = as.matrix(t(z_t))
  rownames(z[[i]]) = esets[[i]]$ENTREZ_GENE_ID
  colnames(z[[i]]) = colnames(df)
  z[[i]] = as.data.frame(z[[i]])
  z[[i]]$EntrezGene.ID = as.character(esets[[i]]$ENTREZ_GENE_ID)
  rm(t, z_t, df)
}; rm(i)
names(z) = names(filt_pdata)

##### Quality Control #####
# Joining in one expression matrix: original version
for (i in 1:length(esets)){
  esets[[i]]$ENTREZ_GENE_ID = as.character(esets[[i]]$ENTREZ_GENE_ID)
  esets[[i]] = esets[[i]][,c("ENTREZ_GENE_ID", 
                             intersect(colnames(esets[[i]]),
                                       filt_pdata[[i]]$GEO_accession))]
}
original_exprs = esets[[1]] %>% inner_join(esets[[2]], by = "ENTREZ_GENE_ID") %>%
  inner_join(esets[[3]], by = "ENTREZ_GENE_ID") %>%
  inner_join(esets[[4]], by = "ENTREZ_GENE_ID") %>%
  dplyr::select(ENTREZ_GENE_ID, everything())

rows = original_exprs$ENTREZ_GENE_ID
original_exprs = as.matrix(original_exprs %>% dplyr::select(-ENTREZ_GENE_ID))
rownames(original_exprs) = rows; rm(rows) # 5016 x 122
# Making sure we do not have NAs in any row
original_exprs = original_exprs[rowSums(is.na(original_exprs)) != ncol(original_exprs), ]
original_exprs_nonas = na.omit(original_exprs) # 5016 x 122

for (i in 1:length(z)){
  z[[i]]$EntrezGene.ID = as.character(z[[i]]$EntrezGene.ID)
  z[[i]] = z[[i]][,c("EntrezGene.ID", 
                     intersect(colnames(z[[i]]),
                               filt_pdata[[i]]$GEO_accession))]
}
# Joining in one expression matrix: z-score normalised version
z_exprs = z[[1]] %>% inner_join(z[[2]], by = "EntrezGene.ID") %>%
  inner_join(z[[3]], by = "EntrezGene.ID") %>%
  inner_join(z[[4]], by = "EntrezGene.ID") %>%
  dplyr::select(EntrezGene.ID, everything())

rownames(z_exprs) = z_exprs$EntrezGene.ID
z_exprs = as.matrix(z_exprs %>% dplyr::select(-EntrezGene.ID)) # 5016 x 122
# Making sure we do not have NAs in any row
z_exprs = z_exprs[rowSums(is.na(z_exprs)) != ncol(z_exprs), ]
z_exprs_nonas = na.omit(z_exprs) # 5015 x 122
# ENTREZ_ID 2597 had some NAs. It was removed from z_exprs_nonas

# Multidimensional scaling plot: original matrix #####
original_mds = plotMDS(original_exprs_nonas)
original_pca = data.frame(cbind(original_mds$x, original_mds$y, 
                                as.character(full_pdata$Study), full_pdata$GEO_accession, 
                                as.character(full_pdata$Tissue_type)))
colnames(original_pca) = c("X1", "X2", "Study", "GSM", "Type")
original_pca$Study = factor(original_pca$Study)
original_pca$Type = factor(original_pca$Type)
original_pca$X1 = as.numeric(original_pca$X1)
original_pca$X2 = as.numeric(original_pca$X2)

tiff("QC/Blood_samples/Original_MDS.tif", width = 1920, height = 1080, res = 100)
ggplot(original_pca, aes(X1, X2, color = Study, shape = Type)) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Dark2") +
  theme(plot.title = element_text(face = "bold", size = 27, hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.text = element_text(angle = 0, hjust = 1, margin = margin(t = 1, unit = "cm"),
                                 size = 15),
        axis.title = element_text(angle = 0, hjust = 0.5, margin = margin(t = 3, unit = "cm"),
                                  size = 20),
        axis.line = element_line(),
        legend.position = "right",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 17),
        legend.key.size = unit(1, "cm"))+
  labs(title = "Multidimensional Scaling Plot",
       x = paste0("\nLeading logFC dimension 1 (", round(100*original_mds$var.explained[1],2), "% of variance)"),
       y = paste0("Leading logFC dimension 2 (", round(100*original_mds$var.explained[2],2), "% of variance)\n"))
dev.off()

# Multidimensional scaling plot: z-score normalised matrix
z_mds = plotMDS(z_exprs_nonas)
z_pca = data.frame(cbind(z_mds$x, z_mds$y, 
                         as.character(full_pdata$Study), full_pdata$GEO_accession, 
                         as.character(full_pdata$Tissue_type)))
colnames(z_pca) = c("X1", "X2", "Study", "GSM", "Type")
z_pca$Study = factor(z_pca$Study)
z_pca$Type = factor(z_pca$Type)
z_pca$X1 = as.numeric(z_pca$X1)
z_pca$X2 = as.numeric(z_pca$X2)

tiff("QC/Blood_samples/KBZ_MDS.tif", width = 1920, height = 1080, res = 100)
ggplot(z_pca, aes(X1, X2, color = Study, shape = Type)) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Dark2") +
  theme(plot.title = element_text(face = "bold", size = 27, hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.text = element_text(angle = 0, hjust = 1, margin = margin(t = 1, unit = "cm"),
                                 size = 15),
        axis.title = element_text(angle = 0, hjust = 0.5, margin = margin(t = 3, unit = "cm"),
                                  size = 20),
        axis.line = element_line(),
        legend.position = "right",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 17),
        legend.key.size = unit(1, "cm"))+
  labs(title = "Multidimensional Scaling Plot: z-score-normalised data",
       x = paste0("\nLeading logFC dimension 1 (", round(100*z_mds$var.explained[1],2), "% of variance)"),
       y = paste0("Leading logFC dimension 2 (", round(100*z_mds$var.explained[2],2), "% of variance)\n"))
dev.off()

# Global expression boxplot: original matrix
original_eset = as.data.frame(original_exprs_nonas)
tiff("QC/Blood_samples/Original_boxplot.tif", width = 1920, height = 1080, res = 100)
ggplot(melt(original_eset), aes(x=variable, y=value)) +
  geom_boxplot(outlier.size = 0.4, outlier.shape = 20,
               fill = c(rep("cyan", 30), rep("chartreuse", 36),
                        rep("orange", 50), rep("red", 6)), outlier.alpha = 0.1)+
  scale_y_continuous("Expression", limits = c(0,round(max(melt(original_eset)$value)+1)), 
                     breaks = seq(0,round(max(melt(original_eset)$value)+1), 1))+
  theme(plot.title = element_text(face = "bold", size = 27, hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.text.y = element_text(angle = 0, hjust = 1, margin = margin(t = 1, unit = "cm"),
                                   size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15, 
                                   margin = margin(t = .05, unit = "cm") ),
        axis.title = element_text(angle = 0, hjust = 0.5, margin = margin(t = 1, unit = "cm"),
                                  size = 25, face = "bold"),
        axis.line = element_line())+
  labs(title = "Boxplot of expression",
       x = "\nSamples",
       y = "Expression\n")
dev.off()

# Global expression boxplot: z-score normalised matrix
z_eset = as.data.frame(z_exprs_nonas)
tiff("QC/Blood_samples/KBZ_boxplot.tif", width = 1920, height = 1080, res = 100)
ggplot(melt(z_eset), aes(x=variable, y=value)) +
  geom_boxplot(outlier.size = 0.4, outlier.shape = 20,
               fill = c(rep("cyan", 30), rep("chartreuse", 36),
                        rep("orange", 50), rep("red", 6)), outlier.alpha = 0.1)+
  scale_y_continuous("Expression", limits = c(0,round(max(melt(z_eset)$value)+1)), 
                     breaks = seq(0,round(max(melt(z_eset)$value)+1), 1))+
  theme(plot.title = element_text(face = "bold", size = 27, hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.text.y = element_text(angle = 0, hjust = 1, margin = margin(t = 1, unit = "cm"),
                                   size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15, 
                                   margin = margin(t = .05, unit = "cm") ),
        axis.title = element_text(angle = 0, hjust = 0.5, margin = margin(t = 1, unit = "cm"),
                                  size = 25, face = "bold"),
        axis.line = element_line())+
  labs(title = "Boxplot of expression: z-score-normalised data",
       x = "\nSamples",
       y = "Expression\n")
dev.off()

# Heatmaps
save_pheatmap_png <- function(x, filename, width=2600, height=1800, res = 130) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_tiff <- function(x, filename, width=2600, height=1800, res = 130) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# Original
annotation_for_heatmap = full_pdata[, c("Study", "Tissue_type")]
rownames(annotation_for_heatmap) = full_pdata$GEO_accession

original_dists = as.matrix(dist(t(original_exprs_nonas), method = "manhattan"))

rownames(original_dists) = full_pdata$GEO_accession
hmcol = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255))
colnames(original_dists) <- NULL
diag(original_dists) <- NA

ann_colors <- list(
  Type = c(tumor = "deeppink4", non_tumor = "dodgerblue4"),
  Study = c(GSE125158 = "darkseagreen", GSE49641 = "darkorange",
            GSE74629 = "darkcyan", GSE18670 = "darkred")
)

original_heatmap = pheatmap(t(original_dists), col = hmcol,
                            annotation_col = annotation_for_heatmap,
                            annotation_colors = ann_colors,
                            legend = TRUE,
                            show_rownames = F,
                            show_colnames = F,
                            treeheight_col = 0,
                            legend_breaks = c(min(original_dists, na.rm = TRUE), 
                                              max(original_dists, na.rm = TRUE)), 
                            legend_labels = (c("small distance", "large distance")),
                            main = "Original heatmap")
save_pheatmap_tiff(original_heatmap, "QC/Blood_samples/original_heatmap.tif")

# Z-score version
annotation_for_heatmap = full_pdata[, c("Study", "Tissue_type")]
rownames(annotation_for_heatmap) = full_pdata$GEO_accession

z_dists = as.matrix(dist(t(z_exprs_nonas), method = "manhattan"))

rownames(z_dists) = full_pdata$GEO_accession
hmcol = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255))
colnames(z_dists) <- NULL
diag(z_dists) <- NA

ann_colors <- list(
  Type = c(tumor = "deeppink4", non_tumor = "dodgerblue4"),
  Study = c(GSE125158 = "darkseagreen", GSE49641 = "darkorange",
            GSE74629 = "darkcyan", GSE18670 = "darkred")
)

z_heatmap = pheatmap(t(z_dists), col = hmcol,
                     annotation_col = annotation_for_heatmap,
                     annotation_colors = ann_colors,
                     legend = TRUE,
                     show_rownames = F,
                     show_colnames = F,
                     treeheight_col = 0,
                     legend_breaks = c(min(z_dists, na.rm = TRUE), 
                                       max(z_dists, na.rm = TRUE)), 
                     legend_labels = (c("small distance", "large distance")),
                     main = "Z-score normalisation heatmap")
save_pheatmap_tiff(z_heatmap, "QC/Blood_samples/KBZ_heatmap.tif")

##### Differential Gene Expression (DGEA) - Union #####
# Joining in one expression matrix: z-score normalised version
z_exprs_blood = z[[1]] %>% full_join(z[[2]], by = "EntrezGene.ID") %>%
  full_join(z[[3]], by = "EntrezGene.ID") %>%
  full_join(z[[4]], by = "EntrezGene.ID") %>%
  dplyr::select(EntrezGene.ID, everything())

rownames(z_exprs_blood) = z_exprs_blood$EntrezGene.ID
z_exprs_blood = as.matrix(z_exprs_blood %>% dplyr::select(-EntrezGene.ID)) # 31581 x 128

full_pdata$Study = as.factor(full_pdata$Study)

# Design and contrast matrix
design = model.matrix(~0 + full_pdata$Tissue_type + full_pdata$Study)
colnames(design) = c("non_tumor", "tumor", "GSE18670", "GSE49641", "GSE74629")
rownames(design) = colnames(z_exprs_blood)
cont.matrix = makeContrasts(tumorvsnontumor=tumor-non_tumor, levels=design)

# DGEA
TN_z_fit = lmFit(z_exprs_blood, design)
TN_z_fit2 = contrasts.fit(TN_z_fit, cont.matrix)
TN_z_fit2 = eBayes(TN_z_fit2, robust = TRUE)
TN_z_results = decideTests(TN_z_fit2)
summary(TN_z_results)
TN_z_DE = as.data.frame(topTable(TN_z_fit2, adjust.method = "BH", number = Inf))
TN_z_DE$EntrezGene.ID = rownames(TN_z_DE)

# 7145 DEGs

# Annotation with official gene symbols
official = org.Hs.egSYMBOL
mapped_genes_official = mappedkeys(official)
official_df = as.data.frame(official[mapped_genes_official])
official_df = official_df %>% dplyr::rename(EntrezGene.ID = gene_id, Gene.Symbol = symbol)
official_df$HGNC_Official = "Yes"
official_df = official_df[-which(duplicated(official_df$Gene.Symbol)==T),]
official_df = distinct(official_df)

alias = org.Hs.egALIAS2EG
mapped_genes_alias = mappedkeys(alias)
alias_df = as.data.frame(alias[mapped_genes_alias])
alias_df = alias_df %>% dplyr::rename(EntrezGene.ID = gene_id, Gene.Symbol = alias_symbol)
alias_df = alias_df[-which(alias_df$Gene.Symbol %in% official_df$Gene.Symbol),]
alias_df$HGNC_Official = "No"

ID_Map = rbind(official_df, alias_df) %>% distinct()
ID_Map$EntrezGene.ID = as.numeric(ID_Map$EntrezGene.ID)
ID_Map = ID_Map[order(ID_Map$EntrezGene.ID),] %>%
  dplyr::rename(probe=Gene.Symbol) %>%
  dplyr::select(probe, EntrezGene.ID, HGNC_Official)

# Aliases
aliases_for_join = alias_df %>% dplyr::rename(Alias = Gene.Symbol)
Aliases = official_df %>% inner_join(aliases_for_join,
                                     by = "EntrezGene.ID") %>%
  dplyr::select(Alias, Gene.Symbol, EntrezGene.ID) %>%
  dplyr::rename(probe = Alias, HGNC_Symbol = Gene.Symbol,
                Entrez = EntrezGene.ID) %>%
  distinct()

rm(alias, alias_df, aliases_for_join, official,
   mapped_genes_alias, mapped_genes_official)

ID_Map = ID_Map %>% dplyr::rename(Gene.Symbol = probe)
ID_Map$EntrezGene.ID = as.character(ID_Map$EntrezGene.ID)
TN_z_DE_mapped = TN_z_DE %>% left_join(ID_Map, by = "EntrezGene.ID")
TN_z_DE_mapped$Filter = NA
unmapped = which(is.na(TN_z_DE_mapped$HGNC_Official))
TN_z_DE_mapped$HGNC_Official[unmapped] = "unmapped"
for(i in 1:nrow(TN_z_DE_mapped)){
  if(TN_z_DE_mapped$HGNC_Official[i] == "Yes"){
    TN_z_DE_mapped$Filter[i] = "Keep"
  } else if(length(unique(TN_z_DE_mapped$HGNC_Official[TN_z_DE_mapped$EntrezGene.ID ==
                                                       TN_z_DE_mapped$EntrezGene.ID[i]])) > 1 &&
            TN_z_DE_mapped$HGNC_Official[i] == "No"){
    TN_z_DE_mapped$Filter[i] = "Discard"
  } else if(unique(TN_z_DE_mapped$HGNC_Official[TN_z_DE_mapped$EntrezGene.ID ==
                                                TN_z_DE_mapped$EntrezGene.ID[i]]) == "No"){
    TN_z_DE_mapped$Filter[i] = "Keep"
    TN_z_DE_mapped$Gene.Symbol[i] = TN_z_DE_mapped$EntrezGene.ID[i]
  } else if(TN_z_DE_mapped$HGNC_Official[i] == "unmapped"){
    TN_z_DE_mapped$Gene.Symbol[i] = TN_z_DE_mapped$EntrezGene.ID[i]
    TN_z_DE_mapped$Filter[i] = "Keep"
  }
}

TN_z_DE_mapped = TN_z_DE_mapped %>% 
  dplyr::filter(Filter == "Keep") %>%
  dplyr::select(EntrezGene.ID, Gene.Symbol, everything()) %>%
  dplyr::select(-Filter) %>%
  distinct()
TN_z_DE_mapped = TN_z_DE_mapped[order(TN_z_DE_mapped$adj.P.Val),]
rownames(TN_z_DE_mapped) = TN_z_DE_mapped$EntrezGene.ID
TN_z_DE_mapped = TN_z_DE_mapped %>% dplyr::select(EntrezGene.ID, Gene.Symbol, everything())
write.xlsx(TN_z_DE_mapped, "DGEA/Union/Blood_samples_analysis/Blood_TN_z_DE_topTable.xlsx",
           overwrite = TRUE)

##### Volcano plots #####
# Tumor vs Normal
TN_z_volcano = EnhancedVolcano(TN_z_DE_mapped,
                               lab = TN_z_DE_mapped[, "Gene.Symbol"],
                               x = 'logFC',
                               y = 'adj.P.Val',
                               title = "Tumor vs. Non-tumor (Blood)",
                               pCutoff = 0.05,
                               FCcutoff = 1,
                               col=c('grey', 'pink', 'purple4', 'red4'),
                               colAlpha = 0.7)
tiff("DGEA/Union/Blood_samples_analysis/Blood_TN_z_Volcano.png", width = 1500, height = 1920, res = 100)
TN_z_volcano
dev.off()

##### Comparisons between tumor tissue identified DEGs and blood samples DEGs follow.
# Just for stage 1 comparisons we provide the overlap with the stage 1 vs stage 2 dichotomizers
# too to see if there are DEGs in the blood which could potentially separate stage 1 and stage 2
# tumors #####

##### Comparisons with Stage 1 vs normal analysis #####
# Intersect of DEGs with the tumor vs normal (stage 1 tumors) DEGs:
stage_1_z_results = read.xlsx("DGEA/Union/DGEA_results.xlsx", sheet = 2)
dichotomizers_1_2 = read.xlsx("DGEA/Union/Dichotomizers.xlsx", sheet = 2)
significants_blood = TN_z_DE_mapped$Gene.Symbol[TN_z_DE_mapped$adj.P.Val < 0.05]
significants_stage_1 = stage_1_z_results$Gene.Symbol[stage_1_z_results$adj.P.Val < 0.05]

DEG_overlap_stage_1 = intersect(significants_blood, significants_stage_1)
dichotomizers_1_2_overlap = intersect(DEG_overlap_stage_1, dichotomizers_1_2$Gene.Symbol)
# 2097 genes, 425 of which are included in the dichotomizers_1_2

# We also need to establish which of the overlapping genes are differentially expressed
# towards the same direction (up-/down-regulated):

subset1 = stage_1_z_results[stage_1_z_results$Gene.Symbol %in% DEG_overlap_stage_1, ] %>%
  dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
  dplyr::rename(logFC_stage_1 = logFC, adj_p_val_stage_1 = adj.P.Val)
subset2 = TN_z_DE_mapped[TN_z_DE_mapped$Gene.Symbol %in% DEG_overlap_stage_1, ] %>%
  dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
  dplyr::rename(logFC_blood = logFC, adj_p_val_blood = adj.P.Val)
common_set_stage_1 = inner_join(subset1, subset2, by = c("Gene.Symbol", "EntrezGene.ID"))
common_set_stage_1$concordance = ifelse(common_set_stage_1$logFC_stage_1*common_set_stage_1$logFC_blood > 0, 1, 0)
concordant_set_stage_1 = common_set_stage_1 %>%
  dplyr::filter(concordance == 1)
concordant_set_stage_1 = concordant_set_stage_1[order(concordant_set_stage_1$adj_p_val_blood,
                                                      concordant_set_stage_1$adj_p_val_stage_1), ]
write.xlsx(concordant_set_stage_1, "DGEA/Blood_DEG_overlap_stage_1.xlsx", overwrite = TRUE)

concordant_dichotomizers_1_2_overlap = intersect(concordant_set_stage_1$Gene.Symbol, 
                                              dichotomizers_1_2$Gene.Symbol)
# 1243 final genes, 218 of which are included in the dichotomizers_1_2

##### Comparisons with Stage 2 vs normal analysis #####
# Intersect of DEGs with the tumor vs normal (stage 2 tumors) DEGs:
stage_2_z_results = read.xlsx("DGEA/Union/DGEA_results.xlsx", sheet = 3)
significants_blood = TN_z_DE_mapped$Gene.Symbol[TN_z_DE_mapped$adj.P.Val < 0.05]
significants_stage_2 = stage_2_z_results$Gene.Symbol[stage_2_z_results$adj.P.Val < 0.05]

DEG_overlap_stage_2 = intersect(significants_blood, significants_stage_2)
# 3264 genes

# We also need to establish which of the overlapping genes are differentially expressed
# towards the same direction (up-/down-regulated):

subset1 = stage_2_z_results[stage_2_z_results$Gene.Symbol %in% DEG_overlap_stage_2, ] %>%
  dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
  dplyr::rename(logFC_stage_2 = logFC, adj_p_val_stage_2 = adj.P.Val)
subset2 = TN_z_DE_mapped[TN_z_DE_mapped$Gene.Symbol %in% DEG_overlap_stage_2, ] %>%
  dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
  dplyr::rename(logFC_blood = logFC, adj_p_val_blood = adj.P.Val)
common_set_stage_2 = inner_join(subset1, subset2, by = c("Gene.Symbol", "EntrezGene.ID"))
common_set_stage_2$concordance = ifelse(common_set_stage_2$logFC_stage_2*common_set_stage_2$logFC_blood > 0, 1, 0)
concordant_set_stage_2 = common_set_stage_2 %>%
  dplyr::filter(concordance == 1)
concordant_set_stage_2 = concordant_set_stage_2[order(concordant_set_stage_2$adj_p_val_blood,
                                                      concordant_set_stage_2$adj_p_val_stage_2), ]
write.xlsx(concordant_set_stage_2, "DGEA/Blood_DEG_overlap_stage_2.xlsx", overwrite = TRUE)

# 1825 final genes

##### Comparisons with Stage 3 vs normal analysis #####
# Intersect of DEGs with the tumor vs normal (stage 3 tumors) DEGs:
stage_3_z_results = read.xlsx("DGEA/Union/DGEA_results.xlsx", sheet = 4)
significants_blood = TN_z_DE_mapped$Gene.Symbol[TN_z_DE_mapped$adj.P.Val < 0.05]
significants_stage_3 = stage_3_z_results$Gene.Symbol[stage_3_z_results$adj.P.Val < 0.05]

DEG_overlap_stage_3 = intersect(significants_blood, significants_stage_3)
# 2252 genes

# We also need to establish which of the overlapping genes are differentially expressed
# towards the same direction (up-/down-regulated):

subset1 = stage_3_z_results[stage_3_z_results$Gene.Symbol %in% DEG_overlap_stage_3, ] %>%
  dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
  dplyr::rename(logFC_stage_3 = logFC, adj_p_val_stage_3 = adj.P.Val)
subset2 = TN_z_DE_mapped[TN_z_DE_mapped$Gene.Symbol %in% DEG_overlap_stage_3, ] %>%
  dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
  dplyr::rename(logFC_blood = logFC, adj_p_val_blood = adj.P.Val)
common_set_stage_3 = inner_join(subset1, subset2, by = c("Gene.Symbol", "EntrezGene.ID"))
common_set_stage_3$concordance = ifelse(common_set_stage_3$logFC_stage_3*common_set_stage_3$logFC_blood > 0, 1, 0)
concordant_set_stage_3 = common_set_stage_3 %>%
  dplyr::filter(concordance == 1)
concordant_set_stage_3 = concordant_set_stage_3[order(concordant_set_stage_3$adj_p_val_blood,
                                                      concordant_set_stage_3$adj_p_val_stage_3), ]
write.xlsx(concordant_set_stage_3, "DGEA/Blood_DEG_overlap_stage_3.xlsx", overwrite = TRUE)

# 1319 final genes

##### Comparisons with Stage 4 vs normal analysis #####
# Intersect of DEGs with the tumor vs normal (stage 4 tumors) DEGs:
stage_4_z_results = read.xlsx("DGEA/Union/DGEA_results.xlsx", sheet = 5)
significants_blood = TN_z_DE_mapped$Gene.Symbol[TN_z_DE_mapped$adj.P.Val < 0.05]
significants_stage_4 = stage_4_z_results$Gene.Symbol[stage_4_z_results$adj.P.Val < 0.05]

DEG_overlap_stage_4 = intersect(significants_blood, significants_stage_4)
# 1877 genes

# We also need to establish which of the overlapping genes are differentially expressed
# towards the same direction (up-/down-regulated):

subset1 = stage_4_z_results[stage_4_z_results$Gene.Symbol %in% DEG_overlap_stage_4, ] %>%
  dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
  dplyr::rename(logFC_stage_4 = logFC, adj_p_val_stage_4 = adj.P.Val)
subset2 = TN_z_DE_mapped[TN_z_DE_mapped$Gene.Symbol %in% DEG_overlap_stage_4, ] %>%
  dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
  dplyr::rename(logFC_blood = logFC, adj_p_val_blood = adj.P.Val)
common_set_stage_4 = inner_join(subset1, subset2, by = c("Gene.Symbol", "EntrezGene.ID"))
common_set_stage_4$concordance = ifelse(common_set_stage_4$logFC_stage_4*common_set_stage_4$logFC_blood > 0, 1, 0)
concordant_set_stage_4 = common_set_stage_4 %>%
  dplyr::filter(concordance == 1)
concordant_set_stage_4 = concordant_set_stage_4[order(concordant_set_stage_4$adj_p_val_blood,
                                                      concordant_set_stage_4$adj_p_val_stage_4), ]
write.xlsx(concordant_set_stage_4, "DGEA/Blood_DEG_overlap_stage_4.xlsx", overwrite = TRUE)

# 1037 final genes

rm(subset1, subset2)

##### Overlap of concordant overlaps #####
all_stage_blood_concordant_overlap_set = inner_join(concordant_set_stage_1,
                                                    concordant_set_stage_2,
                                                    by = c("EntrezGene.ID",
                                                           "Gene.Symbol")) %>%
  inner_join(concordant_set_stage_3, by = c("EntrezGene.ID",
                                            "Gene.Symbol")) %>%
  inner_join(concordant_set_stage_4, by = c("EntrezGene.ID",
                                            "Gene.Symbol"))

# 820 genes
consistent_genes_Entrez = all_stage_blood_concordant_overlap_set$EntrezGene.ID
consistent_genes_symbols = all_stage_blood_concordant_overlap_set$Gene.Symbol

write.xlsx(all_stage_blood_concordant_overlap_set, 
           "DGEA/all_stage_blood_concordant_overlap_set.xlsx",
           overwrite = TRUE)

# Save the volcano after loading the previous volcano plots. Use the same
# .RData file as before:

# rm(list=setdiff(ls(), c("union_four_normal_volcano", "union_three_normal_volcano", 
# "union_two_normal_volcano", "union_one_normal_volcano", "TN_z_volcano")))