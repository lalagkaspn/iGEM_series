## This script is used to download and pre-process series (GSE) from GEO for PDAC patients (normal and tumor tissues) with samples taken
## from blood (whole blood samples or circulating tumor samples)

library(dplyr)
library(GEOquery)
library(org.Hs.eg.db)
library(ggplot2)
library(reshape2)
library(pheatmap)

##### Downloading data #####
datasets = c("GSE125158", "GSE49641", "GSE74629", "GSE18670")

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
names(GEOsets) = datasets; rm(i)

###### pData #####
pdata = list()
for(i in 1:length(GEOsets)) {
  pdata[[i]] = pData(GEOsets[[i]])
}
names(pdata) = datasets; rm(i)

## Isolate only needed pdata
filt_pdata = list()

# GSE125158 
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125158
# Paper: https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC6447845/
# Comments/To discuss with Aris:
#   - 30 samples (30 patients)
#     - 17 PDAC whole blood samples
#     - 13 normal whole blood samples

# Select informative data only
filt_pdata[["GSE125158"]] = pdata$GSE125158 %>%
  select(GEO_accession = geo_accession,
         Patient_ID = title,
         Platform = platform_id,
         Gender = "gender:ch1",
         Age = "age:ch1",
         Tissue_type = "diagnosis:ch1")

# Change Tissue_type to tumor and non_tumor
for (i in 1:length(filt_pdata$GSE125158$Tissue_type)) {
  if(filt_pdata$GSE125158$Tissue_type[i] == "healthy") {
    filt_pdata$GSE125158$Tissue_type[i] = "non_tumor"
  } else {
    filt_pdata$GSE125158$Tissue_type[i] = "tumor"
  }
}; rm(i)

# Clear Patient_ID names
filt_pdata$GSE125158$Patient_ID = sub("\\:.*", "", filt_pdata$GSE125158$Patient_ID)

# Transform pdata to universal levels
filt_pdata$GSE125158$Tissue_type = factor(x = filt_pdata$GSE125158$Tissue_type, levels = c("non_tumor","tumor"), labels = c("non_tumor","tumor"))
filt_pdata$GSE125158$Gender = factor(x = filt_pdata$GSE125158$Gender, levels = c("female","male"), labels = c("female","male"))
filt_pdata$GSE125158$Age = as.numeric(filt_pdata$GSE125158$Age)

# GSE49641
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49641
# Paper: https://pubmed.ncbi.nlm.nih.gov/25069573/
# Comments/To discuss with Aris:
#   - 36 samples (36 patients)
#     - 18 unresectable PDAC
#     - 18 normal

# Select informative data only
filt_pdata[["GSE49641"]] = pdata$GSE49641 %>%
  select(GEO_accession = geo_accession,
         Patient_ID = title,
         Platform = platform_id,
         Gender = "gender:ch1",
         Age = "age:ch1",
         Tissue_type = "disease state:ch1",
         AJCC_classification = "pathological staging:ch1")

# Change Tissue_type to tumor and non_tumor
for (i in 1:length(filt_pdata$GSE49641$Tissue_type)) {
  if(filt_pdata$GSE49641$Tissue_type[i] == "control") {
    filt_pdata$GSE49641$Tissue_type[i] = "non_tumor"
  } else {
    filt_pdata$GSE49641$Tissue_type[i] = "tumor"
  }
}; rm(i)

# Transform pdata to universal levels
filt_pdata$GSE49641$Tissue_type = factor(x = filt_pdata$GSE49641$Tissue_type, levels = c("non_tumor","tumor"), labels = c("non_tumor","tumor"))
filt_pdata$GSE49641$Gender = factor(x = filt_pdata$GSE49641$Gender, levels = c("Female","Male"), labels = c("female","male"))
filt_pdata$GSE49641$Age = as.numeric(filt_pdata$GSE49641$Age)
filt_pdata$GSE49641$AJCC_classification = gsub("IV", "4", filt_pdata$GSE49641$AJCC_classification)
filt_pdata$GSE49641$AJCC_classification = gsub("III", "3", filt_pdata$GSE49641$AJCC_classification)
filt_pdata$GSE49641$AJCC_classification = gsub("NA", NA, filt_pdata$GSE49641$AJCC_classification)
filt_pdata$GSE49641$AJCC_classification = factor(x = filt_pdata$GSE49641$AJCC_classification, levels = c("1a","1b","2a","2b","3","4"), labels = c("1a","1b","2a","2b","3","4"))

# GSE74629
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74629
# Paper: https://pubmed.ncbi.nlm.nih.gov/29617451/
# Comments/To discuss with Aris:
#   - 50 samples (50 patients)
#     - 36 PDAC
#     - 14 normal

# Select informative data only
filt_pdata[["GSE74629"]] = pdata$GSE74629 %>%
  select(GEO_accession = geo_accession,
         Patient_ID = title,
         Platform = platform_id,
         Gender = "gender:ch1",
         Age = "age:ch1",
         Tissue_type = "diagnosis:ch1",
         T_stage = "pathological staging:ch1",
         Chronic_pancreatitis = "chronic pancreatitis:ch1",
         Type_II_diabetes = "diabetes:ch1")
# Patient_ID are named in similar way with GSE49641

# Change Tissue_type to tumor and non_tumor
for (i in 1:length(filt_pdata$GSE74629$Tissue_type)) {
  if(filt_pdata$GSE74629$Tissue_type[i] == "healthy") {
    filt_pdata$GSE74629$Tissue_type[i] = "non_tumor"
  } else {
    filt_pdata$GSE74629$Tissue_type[i] = "tumor"
  }
}; rm(i)

# Transform pdata to universal levels
filt_pdata$GSE74629$Tissue_type = factor(x = filt_pdata$GSE74629$Tissue_type, levels = c("non_tumor","tumor"), labels = c("non_tumor","tumor"))
filt_pdata$GSE74629$Gender = factor(x = filt_pdata$GSE74629$Gender, levels = c("female","male"), labels = c("female","male"))
filt_pdata$GSE74629$Age = as.numeric(filt_pdata$GSE74629$Age)
filt_pdata$GSE74629$Chronic_pancreatitis = as.factor(filt_pdata$GSE74629$Chronic_pancreatitis)
filt_pdata$GSE74629$Type_II_diabetes = as.factor(filt_pdata$GSE74629$Type_II_diabetes)
filt_pdata$GSE74629$T_stage = factor(x = filt_pdata$GSE74629$T_stage, levels = c("1","2","3","4"), labels = c("1","2","3","4"))

# GSE18670
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
  select(GEO_accession = geo_accession,
         Patient_ID = title,
         Platform = platform_id,
         Tissue_type = "source_name_ch1",
         Age = "age:ch1",
         AJCC_classification = "ajcc classif. (2002):ch1",
         Gender = "gender:ch1",
         N_stage = "n-stage:ch1",
         T_stage = "t-stage:ch1",
         Tumor_grade = "tumor grade:ch1")

# Keep only circulating tumor (CTC) and control (P) samples
patients_keep = grep("_CTC_",x = filt_pdata$GSE18670$Tissue_type)
patients_keep = c(patients_keep, grep("_P_",x = filt_pdata$GSE18670$Tissue_type))

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
  select(GEO_accession, Patient_ID, Platform, Gender, Age, Tissue_type, AJCC_classification, T_stage, N_stage, M_stage, Tumor_grade)

# Transform pdata to universal levels
filt_pdata$GSE18670$Gender = factor(x = filt_pdata$GSE18670$Gender, levels = c("Female", "Male"), labels = c("female", "male"))
filt_pdata$GSE18670$Age = as.numeric(filt_pdata$GSE18670$Age)
filt_pdata$GSE18670$Tissue_type = factor(x = filt_pdata$GSE18670$Tissue_type, levels = c("non_tumor", "tumor"), labels = c("non_tumor", "tumor"))
filt_pdata$GSE18670$AJCC_classification = factor(x = filt_pdata$GSE18670$AJCC_classification, levels = c("1a","1b","2a","2b","3","4"), labels = c("1a","1b","2a","2b","3","4"))
filt_pdata$GSE18670$T_stage = factor(x = filt_pdata$GSE18670$T_stage, levels = c("1","2","3","4"), labels = c("1","2","3","4"))
filt_pdata$GSE18670$N_stage = factor(x = filt_pdata$GSE18670$N_stage, levels = c("0","1","2"), labels = c("0","1","2"))
filt_pdata$GSE18670$M_stage = factor(x = filt_pdata$GSE18670$M_stage, levels = c("0", "1"), labels = c("non_metastatic", "metastatic"))
filt_pdata$GSE18670$Tumor_grade = factor(x = filt_pdata$GSE18670$Tumor_grade, levels = c("1","2","3"), labels = c("1","2","3"))

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
pdata74629 = filt_pdata$GSE74629 %>%
  dplyr::select(GEO_accession, Tissue_type) %>%
  dplyr::mutate(Study = "GSE74629")
pdata18670 = filt_pdata$GSE18670 %>%
  dplyr::select(GEO_accession, Tissue_type) %>%
  dplyr::mutate(Study = "GSE18670")

full_pdata = rbind(pdata125158, pdata74629, pdata18670)
rownames(full_pdata) = full_pdata$GSM
rm(pdata125158, pdata74629, pdata18670)

##### Expression data #####
esets = list()
for (i in 1:length(GEOsets)) {
  esets[[i]] = exprs(GEOsets[[i]])
  esets[[i]] = as.data.frame(esets[[i]])
}; rm(i)
names(esets) = datasets

## Calculate NAs values
na_esets = c()
for(i in 1:length(esets)) {
  na_esets[i] = sum(is.na(esets[[i]]))
}; rm(i)
names(na_esets) = datasets
na_esets

## Remove unneeded/inappropriate samples
# GSE18670 
# Keep only samples for tumor (CTC) and non_tumor (P)
esets[["GSE18670"]] = esets$GSE18670[, filt_pdata$GSE18670$GEO_accession]

##### Annotation esets with Entrez ID's #####
# GSE125158
fdata125158 = fData(GEOsets$GSE125158) %>%
  dplyr::select(ID, ENTREZ_GENE_ID = GENE) %>%
  dplyr::filter(!grepl("///", ENTREZ_GENE_ID)) %>%
  dplyr::filter(nchar(ENTREZ_GENE_ID)>0)
esets[[1]]$ID = rownames(esets[[1]])
esets[[1]] = esets[[1]] %>%
  inner_join(fdata125158) %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  summarise_all(mean, na.rm = TRUE)
esets[[1]]$ENTREZ_GENE_ID = as.numeric(esets[[1]]$ENTREZ_GENE_ID)
rm(fdata125158)

# GSE49641
# Convert GB_LIST (RefSeq) to ENTREZ IDs
ref = org.Hs.egREFSEQ2EG
mapped_genes_official = mappedkeys(ref)
ref_df = as.data.frame(ref[mapped_genes_official])
ref_df = ref_df %>% dplyr::rename(EntrezGene.ID = gene_id, RefSeq = accession)
length(unique(ref_df$RefSeq)) == length(ref_df$RefSeq) # TRUE --> No duplicates in RefSeq

# From GB_LIST remove:
#   - Values that have comma (,) (multiple genes to a probe)
#   - Everything after the bullet (.)
fdata49641 = fData(GEOsets$GSE49641) %>%
  dplyr::select(ID, RefSeq = GB_LIST) %>%
  dplyr::filter(!grepl(",", RefSeq)) %>%
  dplyr::filter(nchar(RefSeq)>0) %>%
  mutate(RefSeq = gsub("\\..*","", RefSeq)) %>%
  inner_join(ref_df, by = "RefSeq") %>%
  dplyr::select(ID, ENTREZ_GENE_ID = EntrezGene.ID)
esets[[2]]$ID = as.numeric(rownames(esets[[2]]))
esets[[2]] = esets[[2]] %>%
  inner_join(fdata49641) %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  summarise_all(mean, na.rm = TRUE)
esets[[2]]$ENTREZ_GENE_ID = as.numeric(esets[[2]]$ENTREZ_GENE_ID)
rm(fdata49641)

# GSE74629
fdata74629 = fData(GEOsets$GSE74629) %>%
  dplyr::select(ID, ENTREZ_GENE_ID = Entrez_Gene_ID) %>%
  dplyr::filter(!grepl("///", ENTREZ_GENE_ID)) %>%
  dplyr::filter(nchar(ENTREZ_GENE_ID)>0)
esets[[3]]$ID = rownames(esets[[3]])
esets[[3]] = esets[[3]] %>%
  inner_join(fdata74629) %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  summarise_all(mean, na.rm = TRUE)
esets[[3]]$ENTREZ_GENE_ID = as.numeric(esets[[3]]$ENTREZ_GENE_ID)
rm(fdata74629)

# GSE18670
fdata18670 = fData(GEOsets$GSE18670) %>%
  dplyr::select(ID, ENTREZ_GENE_ID) %>%
  dplyr::filter(!grepl("///", ENTREZ_GENE_ID)) %>%
  dplyr::filter(nchar(ENTREZ_GENE_ID)>0)
esets[[4]]$ID = rownames(esets[[4]])
esets[[4]] = esets[[4]] %>%
  inner_join(fdata18670) %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  summarise_all(mean, na.rm = TRUE)
esets[[4]]$ENTREZ_GENE_ID = as.numeric(esets[[4]]$ENTREZ_GENE_ID)
rm(fdata18670)

##### z-score-transformation #####
# KBZ transformation method ( https:://www.biostars.org/p/283083/ )
z = list()
for(i in 1:length(esets)){
  df = as.data.frame(esets[[i]][,colnames(esets[[i]]) %in% filt_pdata[[i]]$GEO_accession])
  t = as.data.frame(t(df))
  z_t = sapply(t, function(t) (t-mean(t, na.rm = T))/sd(t, na.rm = T))
  z[[i]] = as.matrix(t(z_t))
  rownames(z[[i]]) = esets[[i]]$ENTREZ_GENE_ID
  colnames(z[[i]]) = colnames(df)
  z[[i]] = as.data.frame(z[[i]])
  z[[i]]$EntrezGene.ID = esets[[i]]$ENTREZ_GENE_ID
  rm(t, z_t, df)
}; rm(i)
names(z) = names(filt_pdata)

##### Quality Control #####
# Joining in one expression matrix
# z_exprs = z[[1]] %>% inner_join(z[[2]], by = "EntrezGene.ID") %>%
#   inner_join(., z[[3]], by = "EntrezGene.ID") %>%
#   inner_join(., z[[4]], by = "EntrezGene.ID") %>%
#   dplyr::select(EntrezGene.ID, everything())
# 73 genes

# Exclude z[[2]] (GSE49641)
z_exprs = z[[1]] %>% inner_join(z[[3]], by = "EntrezGene.ID") %>%
  inner_join(., z[[4]], by = "EntrezGene.ID") %>%
  dplyr::select(EntrezGene.ID, everything())

rownames(z_exprs) = z_exprs$EntrezGene.ID
z_exprs = as.matrix(z_exprs %>% dplyr::select(-EntrezGene.ID))
# Making sure we do not have NAs in any row
z_exprs = z_exprs[rowSums(is.na(z_exprs)) != ncol(z_exprs), ]
z_exprs_nonas = na.omit(z_exprs)
# ENTREZ_ID 2597 had some NAs. It was removed from z_ezprs_nonas

# Multidimensional plotting
z_mds = plotMDS(z_exprs_nonas)
z_pca = data.frame(cbind(z_mds$x, z_mds$y, 
                         as.character(full_pdata$Study), full_pdata$GEO_accession, 
                         as.character(full_pdata$Tissue_type)))
colnames(z_pca) = c("X1", "X2", "Study", "GSM", "Type")
z_pca$Study = factor(z_pca$Study, labels = c("GSE125158", "GSE74629", "GSE18670"))
z_pca$Type = factor(z_pca$Type, labels = c("tumor", "non_tumor"))
z_pca$X1 = as.numeric(z_pca$X1)
z_pca$X2 = as.numeric(z_pca$X2)

png("Plots/QC/KBZ/Blood_samples/MDS.png", width = 1024, height = 768)
ggplot(z_pca, aes(X1, X2, color = Type, shape = Study)) +
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
       x = paste0("\nLeading logFC dimension 1 (", round(100*z_mds$var.explained[1],2), "% of variance)"),
       y = paste0("Leading logFC dimension 2 (", round(100*z_mds$var.explained[2],2), "% of variance)\n"))
dev.off()

# Density boxplot - leave for now
# library(oligo)
# o = oligo::boxplot(z_exprs_nonas, main = "Density boxplot",
#                target = "core")

# Global expression boxplot
z_eset = as.data.frame(z_exprs_nonas)
png("Plots/QC/KBZ/Blood_samples/Boxplot.png", width = 1920, height = 1080)
ggplot(melt(z_eset), aes(x=variable, y=value)) +
  geom_boxplot(outlier.size = 0.4, outlier.shape = 20,
               fill = c(rep("cyan", 30), rep("orange", 12), rep("red", 50)), outlier.alpha = 0.1)+
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
  labs(title = "Boxplot of expression",
       x = "\nSamples",
       y = "Expression\n")
dev.off()

# Heatmap
save_pheatmap_png <- function(x, filename, width=2600, height=1800, res = 130) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

annotation_for_heatmap = full_pdata[, c("Study", "Tissue_type")]
rownames(annotation_for_heatmap) = full_pdata$GEO_accession

z_dists = as.matrix(dist(t(z_exprs_nonas), method = "manhattan"))

rownames(z_dists) = full_pdata$GEO_accession
hmcol = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255))
colnames(z_dists) <- NULL
diag(z_dists) <- NA

ann_colors <- list(
  Type = c(tumor = "deeppink4", non_tumor = "dodgerblue4"),
  Study = c(GSE125158 = "darkseagreen", GSE74629 = "darkcyan", GSE18670 = "darkred")
)

z_heatmap = pheatmap(t(z_dists), col = hmcol,
                     annotation_col = annotation_for_heatmap,
                     annotation_colors = ann_colors,
                     legend = TRUE,
                     treeheight_col = 0,
                     legend_breaks = c(min(z_dists, na.rm = TRUE), 
                                       max(z_dists, na.rm = TRUE)), 
                     legend_labels = (c("small distance", "large distance")),
                     main = "Clustering heatmap")
save_pheatmap_png(z_heatmap, "Plots/QC/KBZ/Blood_samples/Heatmap.png")
