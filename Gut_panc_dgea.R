### Gut journal PDAC/pancreatitis exosome RNA paper external validation ###

# Libraries #####
library(dplyr)
library(DESeq2)
library(edgeR)
library(data.table)

# Preprocessing #####

# Count expression data:  not available in GEO
# TPM expression data:  available as .txt file here:
# https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE133684&format=file&file=GSE133684%5Fexp%5FTPM%2Dall%2Etxt%2Egz
# The file is gitignored in the repo due to its size
# Pheno data: available

# However, only useful pheno data is disease category (Normal, PDAC, chr. pancr.)
# which can be inferred from the sample names in the TPM matrix

# Load TPM data into R Studio
TPM = fread("GSE133684_exp_TPM-all.txt")
pdata = as.data.frame(list(Sample.ID = colnames(TPM[,-1])))
cp_greps = grepl("CP_", pdata$Sample.ID)
pdac_greps = grepl("PDAC_", pdata$Sample.ID)
norm_greps = grepl("N_", pdata$Sample.ID)
pdata$type = NA
pdata$description = NA
pdata$type[cp_greps] = "CP"
pdata$description[cp_greps] = "Chronic pancreatitis"
pdata$type[pdac_greps] = "PDAC"
pdata$description[pdac_greps] = "Pancreatic ductal adenocarcinoma"
pdata$type[norm_greps] = "Healthy"
pdata$description[norm_greps] = "Healthy samples"
rm(cp_greps, pdac_greps, norm_greps); gc()

# Convert TPM to matrix format and use first column as row names
genes = TPM$V1
TPM = as.matrix(TPM[, -1])
class(TPM) = "numeric"
rownames(TPM) = genes
log2tpm = log2(TPM + 1)
class(log2tpm) = "numeric"
rm(genes)

# Convert the type column into factor type with base level: Healthy
pdata$type = factor(pdata$type)
pdata$type = relevel(pdata$type, ref = "Healthy")

# Map the ENSEMBL ids to Entrez IDs and Gene Symbols
library(org.Hs.eg.db)

# ENSEMBL and Gene Symbol objects
official = org.Hs.egSYMBOL
ensembl = org.Hs.egENSEMBL2EG

# Mapping keys
mapped_genes_official = mappedkeys(official)
mapped_genes_ensembl = mappedkeys(ensembl)

# Annotation df: ensembl - entrez
ensembl_df = as.data.frame(ensembl[mapped_genes_ensembl])
ensembl_df = ensembl_df %>% dplyr::rename(EntrezGene.ID = gene_id)
ensembl_df = ensembl_df[-which(duplicated(ensembl_df$ensembl_id)), ]

# Annotation df: symbols - entrez
official_df = as.data.frame(official[mapped_genes_official])
official_df = official_df %>% dplyr::rename(EntrezGene.ID = gene_id, Gene.Symbol = symbol)
official_df$HGNC_Official = "Yes"
official_df = official_df[-which(duplicated(official_df$Gene.Symbol)), ]

rm(official, mapped_genes_official, ensembl, mapped_genes_ensembl)

# Exploratory analysis #####

# How many genes in this dataset are all zeros?
table(rowSums(log2tpm==0)==501) # 501: number of samples

# FALSE  TRUE 
# 41606  12542

print(paste0(round(100*12542/54148, 2), "% of genes only have zero counts."))
# 23.16% of genes only have zero counts.

# We need to filter out uninformative genes
# We define a function which is similar to filterByExprs from edgeR

# Define the filter_TPM_matrix function
filter_TPM_matrix <- function(tpm_matrix, group_df, min_log2tpm = 1, min_fraction = 0.5) {
  # Check if the input is a matrix
  if (!is.matrix(tpm_matrix)) {
    stop("Input should be a matrix.")
  }
  
  # Check if the group_df is a data frame with the required columns
  if (!is.data.frame(group_df) || !all(c("Sample.ID", "group") %in% colnames(group_df))) {
    stop("The group argument should be a data frame with columns 'Sample.ID' and 'group'.")
  }
  
  # Check if the colnames of the matrix match the Sample.ID column in the group_df
  if (!all(colnames(tpm_matrix) %in% group_df$Sample.ID)) {
    stop("The colnames of the matrix should match the Sample.ID column in the group data frame.")
  }
  
  # Get unique groups
  unique_groups <- unique(group_df$group)
  
  # Calculate the number of samples required to pass the filter within each group
  min_samples_by_group <- round(tapply(group_df$Sample.ID, group_df$group, length) * min_fraction)
  
  # Define a filtering function for log2(TPM + 1) data
  filter_fun <- function(x, min_samples) {
    sum(x >= min_log2tpm) >= min_samples
  }
  
  # Apply the filter to the input matrix, considering each group
  filter_result_list <- lapply(unique_groups, function(g) {
    group_samples <- group_df$Sample.ID[group_df$group == g]
    matrix_columns <- match(group_samples, colnames(tpm_matrix))
    apply(tpm_matrix[, matrix_columns], 1, function(x) filter_fun(x, min_samples_by_group[g]))
  })
  
  # Combine the filter results using the logical OR operation
  combined_filter_result <- Reduce("|", filter_result_list)
  
  # Apply the combined filter to the input matrix
  filtered_data <- tpm_matrix[combined_filter_result, ]
  
  # Return the filtered matrix
  return(filtered_data)
}

group_df = pdata[, c(1,2)] %>% dplyr::rename(group = type)
group_df$group = as.character(group_df$group)

# min_fraction = 0.1 means that at least 10% of the samples of each group
# have a minimum log2tpm value of 1. The smallest group here is the chronic
# pancreatitis group which has 100 samples. PDAC has 284 and healthy has 117.
# Therefore, genes retained for DGEA have at least 10 CP samples, 28 PDAC samples
# or 11 healthy samples with expression of >= 1 log2(TPM + 1).

log2tpm_filt = filter_TPM_matrix(log2tpm, group_df = group_df,
                                 min_log2tpm = 1, min_fraction = 0.5)
dim(log2tpm_filt)
# 12756   501

class(log2tpm_filt) = "numeric"
rm(group_df)

# Density curve function
library(ggplot2)
library(tidyr)
create_density_curve <- function(matrix) {
  as.data.frame(matrix) %>%
    gather(.) %>%
    ggplot(aes(x = value, fill = "Density Curve")) + 
    geom_density(alpha = 0.5, position = "identity") +
    xlim(-1, NA) +
    labs(title = "Distribution of values by sample",
         x = "Value", y = "Density") +
    theme_classic()+
    theme(legend.position = "none")
}

# Per column density plot
create_density_plot_color <- function(matrix) {
  as.data.frame(matrix) %>%
    gather(., key = "sample", value = "value") %>%
    ggplot(aes(x = value, fill = sample)) + 
    geom_density(alpha = 0.5, position = "identity") +
    xlim(-1, NA) +
    labs(title = "Distribution of values by sample",
         x = "Value", y = "Density") +
    theme_classic()+
    theme(legend.position = "none")
}

dir.create("Exosome_PDAC_plots")

# Density plots/curves for our data
create_density_curve(log2tpm_filt)
ggsave(filename = "ExosomePDAC_filtered_log2(TPM+1)_density_curves.tiff",
       path = "Exosome_PDAC_plots",
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 150, compression = "lzw")
dev.off()

# Color density plots for counts
create_density_plot_color(log2tpm_filt)
ggsave(filename = "ExosomePDAC_filtered_log2(TPM+1)_color_densities.tiff",
       path = "Exosome_PDAC_plots",
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 150, compression = "lzw")
dev.off()

# Produce mean - sd plots 
library(vsn)
library(ggpubr)
library(cowplot)

plot1 = meanSdPlot(log2tpm_filt)$gg + 
  ggtitle(bquote(bold("log2(TPM + 1)"))) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5))
dev.off()
tiff("Exosome_PDAC_plots/meanSdPlot.tiff", res = 300, 
     width = 3200, height = 2000, compression = "lzw")
plot1
dev.off()

# Heatmap of the count matrix - top 100 genes
library(pheatmap)
library(rcartocolor)
select = order(rowMeans(log2tpm_filt),
               decreasing=TRUE)[1:50]
df = as.data.frame(pdata[, c(1,2)])
colnames(df) = c("Sample.ID", "type")
df = df[order(df$type), ]
rownames(df) = df$Sample.ID
df = df %>% dplyr::select(-Sample.ID)
log2tpm_filt = log2tpm_filt[, rownames(df)]

# Annotation colors
ann_colors = list(
  type = c(CP = carto_pal(n = 7, "Emrld")[5], 
                PDAC = "deeppink4", 
                Healthy = "dodgerblue2")
)

# log2(TPM + 1)
tiff("Exosome_PDAC_plots/top50_heatmap_log2(counts+1).tiff",
     res = 200, width = 1920, 
     height = 1080, compression = "lzw")
pheatmap(log2tpm_filt[select, ], cluster_rows=FALSE, show_rownames=FALSE,
         show_colnames = FALSE,
         cluster_cols=FALSE, annotation_col=df, main = "log2(TPM + 1)",
         annotation_colors = ann_colors, fontsize = 6)
dev.off()

# Heatmap of sample distances
library(RColorBrewer)
library(tibble)

sampledists = dist(t(log2tpm_filt))
sampleDistMatrix = as.matrix(sampledists)
#colnames(sampleDistMatrix) = NULL
colors = colorRampPalette(viridisLite::magma(10))(255)

tiff("Exosome_PDAC_plots/sample_distance_heatmap.tiff", 
     res = 700, width = 6760, 
     height = 4760, compression = "lzw")
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampledists,
         clustering_distance_cols = sampledists,
         col = colors,
         show_rownames = FALSE,
         show_colnames = FALSE,
         cutree_rows = 3,
         cutree_cols = 3,
         treeheight_row = 3,
         treeheight_col = 3,
         fontsize = 8,
         main = "Sample distance heatmap",
         annotation_col = df,
         annotation_colors = ann_colors)
dev.off()

# Principal Component Analysis
library(M3C)
type = as.character(df$type)
PCA = pca(log2tpm_filt,
    labels = type,
    legendtextsize = 2, legendtitle = "Sample type", axistextsize = 4,
    dotsize = 0.1)+
  aes(shape = type, size = type,
      alpha = type, color = type) +
  scale_size_manual(name = "Sample type", values = c(0.1, 0.1, 0.1), 
                    labels = c("CP", "Healthy", "PDAC")) +
  scale_alpha_manual(name = "Sample type", values = c(0.7, 0.7, 0.7), 
                     labels = c("CP", "Healthy", "PDAC")) +
  scale_shape_manual(name = "Sample type", values = c(16, 17, 18), 
                     labels = c("CP", "Healthy", "PDAC"))+
  scale_color_manual(name = "Sample type", values = c(carto_pal(n = 7, "Emrld")[5],
                                                      "dodgerblue2",
                                                      "deeppink4"), 
                     labels = c("CP", "Healthy", "PDAC"))+
  theme(panel.border = element_rect(linewidth = 0.2),
        plot.title = element_text(size = 5, face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(size = 4, face = "bold"),
        axis.title.y = element_text(size = 4, face = "bold"),
        axis.ticks = element_line(linewidth = 0.15),
        legend.background = element_rect(fill = "white", linetype = "solid"),
        legend.key.size = unit(0.3, "cm"),
        legend.spacing.y = unit(0.5, units = "mm"))+
  labs(title = "PCA plot: sample types, log2(TPM + 1)") +
  guides(size = "none", alpha = "none")
PCA
ggsave(filename = "PCA.tiff",
       path = "Exosome_PDAC_plots", 
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 700, compression = "lzw")
dev.off()

# Multidimensional scaling plots
mds <- as.data.frame(df[rownames(sampleDistMatrix),])  %>%
  cbind(cmdscale(sampleDistMatrix))
colnames(mds)[1] = "Sample type"
mdsplot = ggplot(mds, aes(x = `1`, y = `2`, color = `Sample type`, 
                          shape = `Sample type`)) +
  geom_point(size = 0.3, alpha = 0.65) +
  ggtitle("Multidimensional scaling: log2(TPM + 1)") +
  theme_classic() +
  scale_color_manual(name = "Sample type", values = c(carto_pal(n = 7, "Emrld")[5],
                                                      "dodgerblue2",
                                                      "deeppink4"), 
                     labels = c("CP", "Healthy", "PDAC"))+
  theme(plot.title = element_text(size = 4, face = "bold", vjust = 0.5, hjust = 0.5),
        axis.text = element_text(size = 3, hjust = 0.5, vjust = 0.5, 
                                 color = "black"),
        axis.title = element_text(size = 4, face = "bold"),
        axis.ticks = element_line(linewidth = 0.1),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "right",
        legend.key.size = unit(2, units = "mm"),
        legend.text = element_text(size = 3),
        legend.title = element_text(face = "bold", size = 3.5),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        legend.spacing.y = unit(0.5, units = "mm"),
        legend.spacing.x = unit(0.5, units = "mm"),
        legend.background = element_blank())+
  labs(x = "MDS1", y = "MDS2") +
  guides(alpha = "none", size = "none", shape = "none")
mdsplot
ggsave(filename = "MDS.tiff",
       path = "Exosome_PDAC_plots", 
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 700, compression = "lzw")
dev.off()

# Differential Gene Expression
design = model.matrix(~0 + df$type)
colnames(design) = c("Healthy", "CP", "PDAC")
rownames(design) = rownames(df)
cm = makeContrasts(CPvsNormal = CP - Healthy,
                   PDACvsNormal = PDAC - Healthy,
                   PDACvsCP = PDAC - CP,
                   PDACvsAll = PDAC - (CP + Healthy)/2,
                   levels = design)
fit = lmFit(log2tpm_filt, design = design)
fit2 = contrasts.fit(fit, contrasts = cm)
fit2 = eBayes(fit2, robust = TRUE)
results = summary(decideTests(fit2))

DE_CPvsNormal = as.data.frame(topTable(fit2, adjust = "BH", number = Inf,
                            coef = "CPvsNormal"))
DE_CPvsNormal$ensembl_id = as.character(rownames(DE_CPvsNormal))
DE_CPvsNormal_mapped = DE_CPvsNormal %>% left_join(ensembl_df, by = "ensembl_id") %>%
  left_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(EntrezGene.ID, ensembl_id, Gene.Symbol, everything())
rownames(DE_CPvsNormal_mapped) = DE_CPvsNormal_mapped$ensembl_id

DE_PDACvsNormal = as.data.frame(topTable(fit2, adjust = "BH", number = Inf,
                                         coef = "PDACvsNormal"))
DE_PDACvsNormal$ensembl_id = as.character(rownames(DE_PDACvsNormal))
DE_PDACvsNormal_mapped = DE_PDACvsNormal %>% left_join(ensembl_df, by = "ensembl_id") %>%
  left_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(EntrezGene.ID, ensembl_id, Gene.Symbol, everything())
rownames(DE_PDACvsNormal_mapped) = DE_PDACvsNormal_mapped$ensembl_id

DE_PDACvsCP = as.data.frame(topTable(fit2, adjust = "BH", number = Inf,
                                     coef = "PDACvsCP"))
DE_PDACvsCP$ensembl_id = as.character(rownames(DE_PDACvsCP))
DE_PDACvsCP_mapped = DE_PDACvsCP %>% left_join(ensembl_df, by = "ensembl_id") %>%
  left_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(EntrezGene.ID, ensembl_id, Gene.Symbol, everything())
rownames(DE_PDACvsCP_mapped) = DE_PDACvsCP_mapped$ensembl_id

DE_PDACvsAll = as.data.frame(topTable(fit2, adjust = "BH", number = Inf,
                                      coef = "PDACvsAll"))
DE_PDACvsAll$ensembl_id = as.character(rownames(DE_PDACvsAll))
DE_PDACvsAll_mapped = DE_PDACvsAll %>% left_join(ensembl_df, by = "ensembl_id") %>%
  left_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(EntrezGene.ID, ensembl_id, Gene.Symbol, everything())
rownames(DE_PDACvsAll_mapped) = DE_PDACvsAll_mapped$ensembl_id

library(openxlsx)
wb = createWorkbook()
addWorksheet(wb, "DE_CPvsNormal")
writeData(wb, "DE_CPvsNormal", DE_CPvsNormal_mapped)
addWorksheet(wb, "DE_PDACvsNormal")
writeData(wb, "DE_PDACvsNormal", DE_PDACvsNormal_mapped)
addWorksheet(wb, "DE_PDACvsCP")
writeData(wb, "DE_PDACvsCP", DE_PDACvsCP_mapped)
addWorksheet(wb, "DE_PDACvsAll")
writeData(wb, "DE_PDACvsAll", DE_PDACvsAll_mapped)

# Comparisons between DEGs #####

# Isolate lists of gene symbols with adj.p.val < 0.05
# Focus on what is specific for PDAC

CPvsNormal_gs_list = DE_CPvsNormal_mapped %>% dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::filter(!is.na(Gene.Symbol))
PDACvsNormal_gs_list = DE_PDACvsNormal_mapped %>% dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::filter(!is.na(Gene.Symbol))
PDACvsCP_gs_list = DE_PDACvsCP_mapped %>% dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::filter(!is.na(Gene.Symbol))
PDACvsAll_gs_list = DE_PDACvsAll_mapped %>% dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::filter(!is.na(Gene.Symbol))

# How many CPvsNormal DEGs overlap with PDACvsNormal DEGs?
length(intersect(CPvsNormal_gs_list$Gene.Symbol, PDACvsNormal_gs_list$Gene.Symbol))
# 749

# The intersect with different directions of logFC could be used for differential diagnosis
diff_diagnosis = CPvsNormal_gs_list %>% dplyr::select(Gene.Symbol, CP_logFC = logFC,
                                                      CP_adj.P.Val = adj.P.Val) %>%
  inner_join(PDACvsNormal_gs_list %>%
               dplyr::select(Gene.Symbol, PDAC_logFC = logFC,
                             PDAC_adj.P.Val = adj.P.Val), by = "Gene.Symbol")
diff_diagnosis$logFC_product = diff_diagnosis$CP_logFC*diff_diagnosis$PDAC_logFC
diff_diagnosis$diagnosis_potential = ifelse(diff_diagnosis$logFC_product < 0,
                                    "diff_diagnosis", "inconclusive_diagnosis")
diff_diagnosis_filt = diff_diagnosis %>%
  dplyr::filter(diagnosis_potential == "diff_diagnosis") %>%
  dplyr::arrange(logFC_product)

# First subset: first ten rows with positive values in PDAC_logFC
# PDAC_top10_upregulated <- diff_diagnosis_filt %>%
#  filter(PDAC_logFC > 0) %>%
#  slice_head(n = 10)

# Second subset: first ten rows with positive values in CP_logFC
# CP_top10_upregulated <- diff_diagnosis_filt %>%
#  filter(CP_logFC > 0) %>%
#  slice_head(n = 10)

# Final diff_diagnosis set
# final_diff_diagnosis = rbind(PDAC_top10_upregulated, CP_top10_upregulated) %>%
#  inner_join(DE_CPvsNormal_mapped %>% dplyr::select(ensembl_id, Gene.Symbol),
#             by = "Gene.Symbol")

final_diff_diagnosis = diff_diagnosis_filt %>%
  inner_join(DE_CPvsNormal_mapped %>% 
               dplyr::select(ensembl_id, Gene.Symbol, EntrezGene.ID), by = "Gene.Symbol",
             multiple = "all")

addWorksheet(wb, "diff_diagnosis")
writeData(wb, "diff_diagnosis", final_diff_diagnosis)
saveWorkbook(wb, "Exosome_PDAC_plots/DGEA.xlsx", overwrite = TRUE)
gc()

# Diff diagnosis dimensionality reduction plots
# PCA
PCA_dd = pca(log2tpm_filt[final_diff_diagnosis$ensembl_id,],
          labels = type,
          legendtextsize = 2, legendtitle = "Sample type", axistextsize = 4,
          dotsize = 0.1)+
  aes(shape = type, size = type,
      alpha = type, color = type) +
  scale_size_manual(name = "Sample type", values = c(0.1, 0.1, 0.1), 
                    labels = c("CP", "Healthy", "PDAC")) +
  scale_alpha_manual(name = "Sample type", values = c(0.7, 0.7, 0.7), 
                     labels = c("CP", "Healthy", "PDAC")) +
  scale_shape_manual(name = "Sample type", values = c(16, 17, 18), 
                     labels = c("CP", "Healthy", "PDAC"))+
  scale_color_manual(name = "Sample type", values = c(carto_pal(n = 7, "Emrld")[5],
                                                      "dodgerblue2",
                                                      "deeppink4"), 
                     labels = c("CP", "Healthy", "PDAC"))+
  theme(panel.border = element_rect(linewidth = 0.2),
        plot.title = element_text(size = 5, face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(size = 4, face = "bold"),
        axis.title.y = element_text(size = 4, face = "bold"),
        axis.ticks = element_line(linewidth = 0.15),
        legend.background = element_rect(fill = "white", linetype = "solid"),
        legend.key.size = unit(0.3, "cm"),
        legend.spacing.y = unit(0.5, units = "mm"))+
  labs(title = "PCA plot: sample types, log2(TPM + 1)") +
  guides(size = "none", alpha = "none")
PCA_dd
ggsave(filename = "PCA_diff_diag.tiff",
       path = "Exosome_PDAC_plots", 
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 700, compression = "lzw")
dev.off()

# MDS
sampledists_dd = dist(t(log2tpm_filt[final_diff_diagnosis$ensembl_id, ]))
sampleDistMatrix_dd = as.matrix(sampledists_dd)
mds_dd <- as.data.frame(df[rownames(sampleDistMatrix_dd),])  %>%
  cbind(cmdscale(sampleDistMatrix_dd))
colnames(mds_dd)[1] = "Sample type"
mdsplot_dd = ggplot(mds_dd, aes(x = `1`, y = `2`, color = `Sample type`, 
                          shape = `Sample type`)) +
  geom_point(size = 0.3, alpha = 0.65) +
  ggtitle("Multidimensional scaling: log2(TPM + 1)") +
  theme_classic() +
  scale_color_manual(name = "Sample type", values = c(carto_pal(n = 7, "Emrld")[5],
                                                      "dodgerblue2",
                                                      "deeppink4"), 
                     labels = c("CP", "Healthy", "PDAC"))+
  theme(plot.title = element_text(size = 4, face = "bold", vjust = 0.5, hjust = 0.5),
        axis.text = element_text(size = 3, hjust = 0.5, vjust = 0.5, 
                                 color = "black"),
        axis.title = element_text(size = 4, face = "bold"),
        axis.ticks = element_line(linewidth = 0.1),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "right",
        legend.key.size = unit(2, units = "mm"),
        legend.text = element_text(size = 3),
        legend.title = element_text(face = "bold", size = 3.5),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        legend.spacing.y = unit(0.5, units = "mm"),
        legend.spacing.x = unit(0.5, units = "mm"),
        legend.background = element_blank())+
  labs(x = "MDS1", y = "MDS2") +
  guides(alpha = "none", size = "none", shape = "none")
mdsplot_dd
ggsave(filename = "MDS_diff_diag.tiff",
       path = "Exosome_PDAC_plots", 
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 700, compression = "lzw")
dev.off()

# GSEA with the diff_diagnosis set #####

# PDAC input
gsea_input_pdac = final_diff_diagnosis %>% dplyr::select(EntrezGene.ID, Gene.Symbol,
                                                         PDAC_logFC, PDAC_adj.P.Val)
log2fc_pdac = gsea_input_pdac$PDAC_logFC
genes = gsea_input_pdac$EntrezGene.ID
gene_list_pdac = log2fc_pdac
names(gene_list_pdac) = genes
rm(log2fc_pdac, genes)

# CP input
gsea_input_cp = final_diff_diagnosis %>% dplyr::select(EntrezGene.ID, Gene.Symbol,
                                                       CP_logFC, CP_adj.P.Val)
log2fc_cp = gsea_input_cp$CP_logFC
genes = gsea_input_cp$EntrezGene.ID
gene_list_cp = log2fc_cp
names(gene_list_cp) = genes
rm(log2fc_cp, genes)

# We need to account for duplicate Entrez IDs

# Define a function that keeps the names for which the absolute logFC is highest, retaining their original sign
get_unique_max_abs <- function(named_list) {
  unique_names <- unique(names(named_list))
  max_abs_values_with_sign <- sapply(unique_names, function(name) {
    values <- named_list[names(named_list) == name]
    max_abs_value_with_sign <- values[which.max(abs(values))]
    return(max_abs_value_with_sign)
  })
  return(setNames(max_abs_values_with_sign, unique_names))
}

gene_list_cp = get_unique_max_abs(gene_list_cp)
gene_list_cp = sort(gene_list_cp, decreasing = TRUE)
gene_list_pdac = get_unique_max_abs(gene_list_pdac)
gene_list_pdac = sort(gene_list_pdac, decreasing = TRUE)


# Prepare output workbook
gsea_wb = createWorkbook()

# Perform GO gene set enrichment analysis using clusterProfiler
library(clusterProfiler)

# PDAC
RNGversion("4.2.2")
set.seed(123)
gseaGO_pdac = gseGO(geneList = gene_list_pdac,
               ont = "ALL",
               OrgDb = "org.Hs.eg.db",
               keyType = "ENTREZID",
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               minGSSize = 10,
               maxGSSize = 800,
               seed = TRUE)

gseaGO_pdac = clusterProfiler::setReadable(gseaGO_pdac, 'org.Hs.eg.db')

addWorksheet(gsea_wb, "PDAC_diff_GSEA_GO")
writeData(gsea_wb, "PDAC_diff_GSEA_GO", as.data.frame(gseaGO_pdac@result))

# KEGG
options(download.file.method = "auto")
RNGversion("4.2.2")
set.seed(123)
gseaKEGG_pdac = gseKEGG(geneList = gene_list_pdac,
                   organism = "hsa",
                   keyType = "ncbi-geneid",
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   minGSSize = 10,
                   maxGSSize = 800,
                   seed = TRUE)

gseaKEGG_pdac = clusterProfiler::setReadable(gseaKEGG_pdac, 'org.Hs.eg.db',
                                             keyType = "ENTREZID")
addWorksheet(gsea_wb, "PDAC_diff_GSEA_KEGG")
writeData(gsea_wb, "PDAC_diff_GSEA_KEGG", as.data.frame(gseaKEGG_pdac@result))

# Reactome
library(ReactomePA)
RNGversion("4.2.2")
set.seed(123)
gseaReactome_pdac = gsePathway(geneList = gene_list_pdac,
                          organism = "human",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          minGSSize = 10,
                          maxGSSize = 800,
                          seed = TRUE)

gseaReactome_pdac = clusterProfiler::setReadable(gseaReactome_pdac, 'org.Hs.eg.db')
addWorksheet(gsea_wb, "PDAC_diff_GSEA_Reactome")
writeData(gsea_wb, "PDAC_diff_GSEA_Reactome", as.data.frame(gseaReactome_pdac@result))

# CP
RNGversion("4.2.2")
set.seed(123)
gseaGO_cp = gseGO(geneList = gene_list_cp,
                  ont = "ALL",
                  OrgDb = "org.Hs.eg.db",
                  keyType = "ENTREZID",
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",
                  minGSSize = 10,
                  maxGSSize = 800,
                  seed = TRUE)

gseaGO_cp = clusterProfiler::setReadable(gseaGO_cp, 'org.Hs.eg.db')

addWorksheet(gsea_wb, "CP_diff_GSEA_GO")
writeData(gsea_wb, "CP_diff_GSEA_GO", as.data.frame(gseaGO_cp@result))

# KEGG
# options(download.file.method = "auto")
RNGversion("4.2.2")
set.seed(123)
gseaKEGG_cp = gseKEGG(geneList = gene_list_cp,
                      organism = "hsa",
                      keyType = "ncbi-geneid",
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      minGSSize = 10,
                      maxGSSize = 800,
                      seed = TRUE)

gseaKEGG_cp = clusterProfiler::setReadable(gseaKEGG_cp, 'org.Hs.eg.db',
                                           keyType = "ENTREZID")
addWorksheet(gsea_wb, "CP_diff_GSEA_KEGG")
writeData(gsea_wb, "CP_diff_GSEA_KEGG", as.data.frame(gseaKEGG_cp@result))

# Reactome
RNGversion("4.2.2")
set.seed(123)
gseaReactome_cp = gsePathway(geneList = gene_list_cp,
                             organism = "human",
                             pvalueCutoff = 0.05,
                             pAdjustMethod = "BH",
                             minGSSize = 10,
                             maxGSSize = 800,
                             seed = TRUE)

gseaReactome_cp = clusterProfiler::setReadable(gseaReactome_cp, 'org.Hs.eg.db')
addWorksheet(gsea_wb, "CP_diff_GSEA_Reactome")
writeData(gsea_wb, "CP_diff_GSEA_Reactome", as.data.frame(gseaReactome_cp@result))
saveWorkbook(gsea_wb, "Exosome_PDAC_plots/GSEA.xlsx",
             overwrite = TRUE)
gc()

# Comparison of the 346 genes with the 820-gene signature from our PDAC work #####
pdac_sig = read.xlsx("DGEA/all_stage_blood_concordant_overlap_set.xlsx")
length(intersect(pdac_sig$EntrezGene.ID, gsea_input_pdac$EntrezGene.ID))
# 27 common Entrez IDs
length(intersect(pdac_sig$Gene.Symbol, gsea_input_pdac$Gene.Symbol))
# 27 common gene symbols

overlap_820_subset = pdac_sig[pdac_sig$Gene.Symbol %in%
                                intersect(pdac_sig$Gene.Symbol, gsea_input_pdac$Gene.Symbol), ] %>%
  left_join(final_diff_diagnosis, by = "Gene.Symbol")

write.xlsx(overlap_820_subset, 
           "Exosome_PDAC_plots/overlap_with_820sig.xlsx")
