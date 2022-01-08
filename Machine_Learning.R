# Machine Learning on Prostate Cancer samples
# Libraries and importing #####

library(dplyr)
library(org.Hs.eg.db)
library(ggplot2)
library(reshape2)
library(tidyr)
library(openxlsx)
library(C50)
library(caret)
library(irr)
library(ipred)
library(adabag)
library(vcd)
library(randomForest)
library(pROC)
library(stringr)
library(kernlab)
library(glmnet)
library(MASS)
library(LiblineaR)

# Importing the two expression data frames
blood_frame = read.xlsx("Blood_samples_z_expression_matrix.xlsx")
tumor_frame = read.xlsx("Tumor_samples_z_expression_matrix.xlsx")

# Importing the blood - tumor DEG concordance set
concordance_set = read.xlsx("DGEA/Blood_Tumor_DEG_overlap.xlsx")

# Subsetting for common DEGs
blood_frame = blood_frame[, c("GEO_accession", "Tissue_type",
                              "Study", concordance_set$Gene.Symbol)] %>%
  mutate(Source = "Blood") %>%
  dplyr::select(GEO_accession, Tissue_type, Study, Source, everything())
tumor_frame = tumor_frame[, c("GEO_accession", "Tissue_type",
                              "Study", concordance_set$Gene.Symbol)] %>%
  mutate(Source = "Tumor") %>%
  dplyr::select(GEO_accession, Tissue_type, Study, Source, everything())

# We will remove the normal samples from the tumor_frame for study GSE18670
# because these are duplicated in blood_frame:
removals = tumor_frame$GEO_accession[tumor_frame$Tissue_type == "non_tumor" &
                                       tumor_frame$Study == "GSE18670"]
tumor_frame = tumor_frame %>%
  dplyr::filter(!GEO_accession %in% removals)
rm(removals)

# Binding into a full expression frame
global_expression = rbind(blood_frame, tumor_frame) # 571 samples x 533 genes
global_expression[colnames(global_expression[,5:ncol(global_expression)])] = 
  lapply(global_expression[colnames(global_expression[,5:ncol(global_expression)])], 
         function(x) as.numeric(x))
global_expression[colnames(global_expression[,2:4])] = 
  lapply(global_expression[colnames(global_expression[,2:4])], 
         function(x) as.factor(x))
global_expression$Tissue_type = relevel(global_expression$Tissue_type,
                                        ref = "non_tumor")
global_expression$Source = relevel(global_expression$Source,
                                        ref = "Tumor")

# Check if gene symbols are syntactically valid for R:
isValidName <- function(string) {
  grepl("^((([[:alpha:]]|[.][._[:alpha:]])[._[:alnum:]]*)|[.])$", string)
}

isValidAndUnreservedName <- function(string) {
  make.names(string) == string
}

testValidity <- function(string) {
  valid <- isValidName(string)
  unreserved <- isValidAndUnreservedName(string)
  reserved <- (valid & ! unreserved)
  list("Valid"=valid,
       "Unreserved"=unreserved,
       "Reserved"=reserved)
}

testNames <- colnames(global_expression)
res = t(sapply(testNames, testValidity)) # All are ok

# Consensus clustering #####
# Preparing data
library(M3C)
rownames(global_expression) = global_expression$GEO_accession
data = global_expression %>%
  dplyr::select(-GEO_accession, -Tissue_type, -Study, -Source)
data = as.data.frame(t(data))
RNGversion("4.0.2")
mock = M3C(mydata = data, des = global_expression, iters = 100, repsref = 250, 
           repsreal = 250, seed = 123, fsize = 20, lthick = 2, dotsize = 1.25)

# optimal K: 2
cluster_number = c(2,3,4,5,6,7,8,9,10)

chifit_resp_p = list()
for (k in seq(2,10)){
  myresults <- mock$realdataresults[[k]]$ordered_annotation
  chifit_resp <- suppressWarnings(chisq.test(table(myresults[c('consensuscluster','class')])))
  chifit_resp_p[[k-1]] = chifit_resp$p.value
}

chifit_time_p = list()
for (k in seq(2,10)){
  myresults <- mock$realdataresults[[k]]$ordered_annotation
  chifit_time <- suppressWarnings(chisq.test(table(myresults[c('consensuscluster','Timepoint_coded')])))
  chifit_time_p[[k-1]] = chifit_time$p.value
}

chifit_pam50_p = list()
for (k in seq(2,10)){
  myresults <- mock$realdataresults[[k]]$ordered_annotation
  chifit_pam50 <- suppressWarnings(chisq.test(table(myresults[c('consensuscluster','pam50')])))
  chifit_pam50_p[[k-1]] = chifit_pam50$p.value
}

chifit_treat_p = list()
for (k in seq(2,10)){
  myresults <- mock$realdataresults[[k]]$ordered_annotation
  chifit_treat <- suppressWarnings(chisq.test(table(myresults[c('consensuscluster','Treatment')])))
  chifit_treat_p[[k-1]] = chifit_treat$p.value
}

chifit_mammaprint_p = list()
for (k in seq(2,10)){
  myresults <- mock$realdataresults[[k]]$ordered_annotation
  chifit_mammaprint <- suppressWarnings(chisq.test(table(myresults[c('consensuscluster','Mammaprint_risk')])))
  chifit_mammaprint_p[[k-1]] = chifit_mammaprint$p.value
}

chifit_rorS_p = list()
for (k in seq(2,10)){
  myresults <- mock$realdataresults[[k]]$ordered_annotation
  chifit_rorS <- suppressWarnings(chisq.test(table(myresults[c('consensuscluster','rorS_risk')])))
  chifit_rorS_p[[k-1]] = chifit_rorS$p.value
}

chifit_scmod1_p = list()
for (k in seq(2,10)){
  myresults <- mock$realdataresults[[k]]$ordered_annotation
  chifit_scmod1 <- suppressWarnings(chisq.test(table(myresults[c('consensuscluster','scmod1')])))
  chifit_scmod1_p[[k-1]] = chifit_scmod1$p.value
}

chifit_IC10_p = list()
for (k in seq(2,10)){
  myresults <- mock$realdataresults[[k]]$ordered_annotation
  chifit_IC10 <- suppressWarnings(chisq.test(table(myresults[c('consensuscluster','IC10')])))
  chifit_IC10_p[[k-1]] = chifit_IC10$p.value
}

cluster_corr = as.data.frame(cbind(cluster_number, chifit_resp_p, chifit_treat_p,
                                   chifit_time_p, chifit_pam50_p, chifit_mammaprint_p,
                                   chifit_rorS_p, chifit_scmod1_p, chifit_IC10_p))
colnames(cluster_corr) = c("No. of clusters", "Response", "Treatment", "Timepoint", "pam50",
                           "Mammaprint_risk", "rorS_risk", "scmod1", "IC10")

# PCA clustering
png("data/Output sets/DGEA_GSEA_ML/kNN/Cluster_PCA.png",
    width = 1920, height = 1080)
pca(mock$realdataresults[[2]]$ordered_data,
    labels = mock$realdataresults[[2]]$ordered_annotation$consensuscluster,
    legendtextsize = 20, axistextsize = 20, dotsize = 5)
dev.off()

# Consensus index plot
png("data/Output sets/DGEA_GSEA_ML/kNN/Consensus_index.png",
    width = 1920, height = 1080)
mock[["plots"]][[1]]
dev.off()

# Entropy plot
png("data/Output sets/DGEA_GSEA_ML/kNN/Entropy.png",
    width = 1920, height = 1080)
mock[["plots"]][[2]]
dev.off()

# Statistical significance of clusters
png("data/Output sets/DGEA_GSEA_ML/kNN/Stat_Sig.png",
    width = 1920, height = 1080)
mock[["plots"]][[3]]
dev.off()

# RCSI plot
png("data/Output sets/DGEA_GSEA_ML/kNN/RCSI.png",
    width = 1920, height = 1080)
mock[["plots"]][[4]]
dev.off()

kmwb = createWorkbook()
addWorksheet(kmwb, "Results")
writeData(kmwb, "Results", mock[["realdataresults"]][[2]][["ordered_annotation"]])
addWorksheet(kmwb, "Chisq.test")
writeData(kmwb, "Chisq.test", cluster_corr)
addWorksheet(kmwb, "Performance")
writeData(kmwb, "Performance", mock[["scores"]])
saveWorkbook(kmwb, file = "data/Output sets/DGEA_GSEA_ML/kNN/Clustering.xlsx",
             overwrite = TRUE); rm(kmwb); gc(verbose = FALSE)

rm(KBZ_Entrez_list, data, genes, full_data); gc(verbose = FALSE)

cluster_results = data.frame(names(mock[["realdataresults"]][[2]][["assignments"]]), 
                             mock[["realdataresults"]][[2]][["assignments"]])
# 2: optimal number of clusters
colnames(cluster_results) = c("Sample.ID", "Cluster")

# Training, Validation, Test set splitting #####
# The sets are generated once and then filtered for variables of interest in
# the ML training, validating, testing process

KBZ_Entrez_list = KBZ_DE_mapped %>% dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::select(EntrezGene.ID, Gene.Symbol)
data = as.data.frame(t(z_exprs))
genes = intersect(colnames(data), KBZ_Entrez_list$EntrezGene.ID)
data = data %>% dplyr::select(all_of(genes))
data$Sample.ID = rownames(data)
data = data %>%
  inner_join(Pheno_no_nit %>% 
               dplyr::select(Sample.ID, Response, pam50, Timepoint_coded,
                             Treatment, Dataset, rorS_risk, scmod1, IC10,
                             Mammaprint_risk), 
             by = "Sample.ID") %>%
  inner_join(cluster_results, by = "Sample.ID")
data$Dataset = as.factor(data$Dataset)
data$Cluster = as.factor(data$Cluster)
data$Timepoint_coded = as.factor(data$Timepoint_coded)
data$Treatment = as.factor(data$Treatment)
str_sub(colnames(data)[1:length(genes)],0,0) = "X_"
data = na.omit(data)
data$Response = as.factor(data$Response)

# Train, Validation, Test splits
RNGversion("4.0.2")
set.seed(123)
tobesplit = data %>%
  mutate(nr = row_number()) %>%
  dplyr::select(nr, everything()) %>%
  as.data.frame()
RNGversion("4.0.2")
set.seed(123)
train_set = tobesplit %>%
  group_by(Treatment, Response, Timepoint_coded, pam50, Dataset, Cluster) %>%
  sample_frac(0.7) %>%
  as.data.frame()
RNGversion("4.0.2")
set.seed(123)
validation_set = anti_join(tobesplit, train_set) %>%
  group_by(Treatment, Response, pam50, Timepoint_coded, Dataset, Cluster) %>%
  sample_frac(0.6) %>%
  as.data.frame()
test_set = anti_join(tobesplit, as.data.frame(rbind(train_set, validation_set)))

train_set = train_set %>% dplyr::ungroup()
validation_set = validation_set %>% dplyr::ungroup()
test_set = test_set %>% dplyr::ungroup()


# Lists of training and validation sets:
# (Because the number of clusters is 2, no models will be trained using just the cluster
# as a predictor. So, correlation with variables of interest will suffice in this case.)
# Training
train_gs = train_set %>% 
  dplyr::select(-nr, -pam50, -Sample.ID, -Timepoint_coded, -Treatment,
                -Dataset, -Cluster, -rorS_risk, -scmod1, -IC10,
                -Mammaprint_risk)
unchanged = c(colnames(train_gs)[c(1:62)], "Response")

train_cluster_meta = train_set %>% 
  dplyr::select(Response, pam50, Timepoint_coded, Treatment, Cluster, IC10,
                Mammaprint_risk, rorS_risk, scmod1)
train_cluster_meta = cbind(train_cluster_meta[,"Response"], as.data.frame(model.matrix(~0 +
                                                                                         train_cluster_meta$pam50 +
                                                                                         train_cluster_meta$Timepoint_coded +
                                                                                         train_cluster_meta$Treatment +
                                                                                         train_cluster_meta$Cluster +
                                                                                         train_cluster_meta$IC10 +
                                                                                         train_cluster_meta$Mammaprint_risk +
                                                                                         train_cluster_meta$rorS_risk +
                                                                                         train_cluster_meta$scmod1)))
colnames(train_cluster_meta) = c("Response", "Basal", "HER2", "LumB", "LumA",
                                 "Normal", "T2", "Endo", "Cluster_2", "IC2", "IC3", "IC4",
                                 "IC5", "IC6", "IC7", "IC8", "IC9", "IC10", "Mammaprint_risk_yes",
                                 "rorS_risk_interm", "rorS_risk_high", "ER_hp", "ER_lp", "HER2_scmod1")
train_cluster_meta = train_cluster_meta %>% dplyr::select(-Basal)
train_cluster_meta[colnames(train_cluster_meta[,2:ncol(train_cluster_meta)])] = 
  lapply(train_cluster_meta[colnames(train_cluster_meta[,2:ncol(train_cluster_meta)])], 
         function(x) factor(x, levels = c(0,1), labels = c(0,1)))

train_gs_meta = train_set %>% 
  dplyr::select(-nr, -Sample.ID, -Cluster)
train_gs_meta = cbind(train_gs_meta[,unchanged], as.data.frame(model.matrix(~0 +
                                                                              train_gs_meta$pam50 +
                                                                              train_gs_meta$Timepoint_coded +
                                                                              train_gs_meta$Treatment +
                                                                              train_gs_meta$IC10 +
                                                                              train_gs_meta$Mammaprint_risk +
                                                                              train_gs_meta$rorS_risk +
                                                                              train_gs_meta$scmod1)))
colnames(train_gs_meta) = c(unchanged, "Basal", "HER2", "LumB", "LumA",
                            "Normal", "T2", "Endo", "IC2", "IC3", "IC4",
                            "IC5", "IC6", "IC7", "IC8", "IC9", "IC10", "Mammaprint_risk_yes",
                            "rorS_risk_interm", "rorS_risk_high", "ER_hp", "ER_lp", "HER2_scmod1")
train_gs_meta = train_gs_meta %>% dplyr::select(-Basal)
train_gs_meta[c("HER2", "LumB", "LumA",
                "Normal", "T2", "Endo", "IC2", "IC3", "IC4",
                "IC5", "IC6", "IC7", "IC8", "IC9", "IC10", "Mammaprint_risk_yes",
                "rorS_risk_interm", "rorS_risk_high", "ER_hp", "ER_lp", "HER2_scmod1")] = 
  lapply(train_gs_meta[c("HER2", "LumB", "LumA", "Normal", "T2", "Endo",
                         "IC2", "IC3", "IC4",
                         "IC5", "IC6", "IC7", "IC8", "IC9", "IC10", "Mammaprint_risk_yes",
                         "rorS_risk_interm", "rorS_risk_high", "ER_hp", "ER_lp", "HER2_scmod1")],
         function(x) factor(x, levels = c(0,1), 
                            labels = c(0,1)))

train_all = train_set %>% 
  dplyr::select(-nr, -Sample.ID)
train_all = cbind(train_all[,unchanged], 
                  as.data.frame(model.matrix(~0 + train_all$pam50 +
                                               train_all$Timepoint_coded +
                                               train_all$Treatment +
                                               train_all$Cluster +
                                               train_all$IC10 +
                                               train_all$Mammaprint_risk +
                                               train_all$rorS_risk +
                                               train_all$scmod1)))

colnames(train_all) = c(unchanged, "Basal", "HER2", "LumB", "LumA",
                        "Normal", "T2", "Endo", "Cluster_2", "IC2", "IC3", "IC4",
                        "IC5", "IC6", "IC7", "IC8", "IC9", "IC10", "Mammaprint_risk_yes",
                        "rorS_risk_interm", "rorS_risk_high", "ER_hp", "ER_lp", "HER2_scmod1")
train_all = train_all %>% dplyr::select(-Basal)
train_all[c("HER2", "LumB", "LumA",
            "Normal", "T2", "Endo", "Cluster_2", "IC2", "IC3", "IC4",
            "IC5", "IC6", "IC7", "IC8", "IC9", "IC10", "Mammaprint_risk_yes",
            "rorS_risk_interm", "rorS_risk_high", "ER_hp", "ER_lp", "HER2_scmod1")] = 
  lapply(train_all[c("HER2", "LumB", "LumA", "Normal", "T2", "Endo",
                     "Cluster_2", "IC2", "IC3", "IC4",
                     "IC5", "IC6", "IC7", "IC8", "IC9", "IC10", "Mammaprint_risk_yes",
                     "rorS_risk_interm", "rorS_risk_high", "ER_hp", "ER_lp", "HER2_scmod1")],
         function(x) factor(x, levels = c(0,1), 
                            labels = c(0,1)))

Training = list(train_gs, train_cluster_meta, train_gs_meta, train_all)
names(Training) = c("GS", "Cluster_meta", "GS_meta", "All_together")
rm(train_gs, train_cluster_meta, train_gs_meta, train_all)

# Validation
validation_gs = validation_set %>% 
  dplyr::select(-nr, -pam50, -Sample.ID, -Timepoint_coded, -Treatment,
                -Dataset, -Cluster, -rorS_risk, -scmod1, -IC10,
                -Mammaprint_risk)
unchanged = c(colnames(validation_gs)[c(1:62)], "Response")

validation_cluster_meta = validation_set %>% 
  dplyr::select(Response, pam50, Timepoint_coded, Treatment, Cluster, IC10,
                Mammaprint_risk, rorS_risk, scmod1)
validation_cluster_meta = cbind(validation_cluster_meta[,"Response"], 
                                as.data.frame(model.matrix(~0 + validation_cluster_meta$pam50 +
                                                             validation_cluster_meta$Timepoint_coded +
                                                             validation_cluster_meta$Treatment +
                                                             validation_cluster_meta$Cluster +
                                                             validation_cluster_meta$IC10 +
                                                             validation_cluster_meta$Mammaprint_risk +
                                                             validation_cluster_meta$rorS_risk +
                                                             validation_cluster_meta$scmod1)))
colnames(validation_cluster_meta) = c("Response", "Basal", "HER2", "LumB", "LumA",
                                      
                                      "Normal", "T2", "Endo", "Cluster_2", "IC2", "IC3", "IC4",
                                      "IC5", "IC6", "IC7", "IC8", "IC9", "IC10", "Mammaprint_risk_yes",
                                      "rorS_risk_interm", "rorS_risk_high", "ER_hp", "ER_lp", "HER2_scmod1")
validation_cluster_meta = validation_cluster_meta %>% dplyr::select(-Basal)
validation_cluster_meta[colnames(validation_cluster_meta[,2:ncol(validation_cluster_meta)])] = 
  lapply(validation_cluster_meta[colnames(validation_cluster_meta[,2:ncol(validation_cluster_meta)])], 
         function(x) factor(x, levels = c(0,1), labels = c(0,1)))

validation_gs_meta = validation_set %>% 
  dplyr::select(-nr, -Sample.ID, -Cluster)
validation_gs_meta = cbind(validation_gs_meta[,unchanged], 
                           as.data.frame(model.matrix(~0 + validation_gs_meta$pam50 +
                                                        validation_gs_meta$Timepoint_coded +
                                                        validation_gs_meta$Treatment +
                                                        validation_gs_meta$IC10 +
                                                        validation_gs_meta$Mammaprint_risk +
                                                        validation_gs_meta$rorS_risk +
                                                        validation_gs_meta$scmod1)))
colnames(validation_gs_meta) = c(unchanged, "Basal", "HER2", "LumB", "LumA",
                                 "Normal", "T2", "Endo", "IC2", "IC3", "IC4",
                                 "IC5", "IC6", "IC7", "IC8", "IC9", "IC10", "Mammaprint_risk_yes",
                                 "rorS_risk_interm", "rorS_risk_high", "ER_hp", "ER_lp", "HER2_scmod1")
validation_gs_meta = validation_gs_meta %>% dplyr::select(-Basal)
validation_gs_meta[c("HER2", "LumB", "LumA",
                     "Normal", "T2", "Endo", "IC2", "IC3", "IC4",
                     "IC5", "IC6", "IC7", "IC8", "IC9", "IC10", "Mammaprint_risk_yes",
                     "rorS_risk_interm", "rorS_risk_high", "ER_hp", "ER_lp", "HER2_scmod1")] = 
  lapply(validation_gs_meta[c("HER2", "LumB", "LumA", "Normal", "T2", "Endo",
                              "IC2", "IC3", "IC4",
                              "IC5", "IC6", "IC7", "IC8", "IC9", "IC10", "Mammaprint_risk_yes",
                              "rorS_risk_interm", "rorS_risk_high", "ER_hp", "ER_lp", "HER2_scmod1")],
         function(x) factor(x, levels = c(0,1), 
                            labels = c(0,1)))

validation_all = validation_set %>% 
  dplyr::select(-nr, -Sample.ID)
validation_all = cbind(validation_all[,unchanged], 
                       as.data.frame(model.matrix(~0 + validation_all$pam50 +
                                                    validation_all$Timepoint_coded +
                                                    validation_all$Treatment +
                                                    validation_all$Cluster +
                                                    validation_all$IC10 +
                                                    validation_all$Mammaprint_risk +
                                                    validation_all$rorS_risk +
                                                    validation_all$scmod1)))

colnames(validation_all) = c(unchanged, "Basal", "HER2", "LumB", "LumA",
                             "Normal", "T2", "Endo", "Cluster_2", "IC2", "IC3", "IC4",
                             "IC5", "IC6", "IC7", "IC8", "IC9", "IC10", "Mammaprint_risk_yes",
                             "rorS_risk_interm", "rorS_risk_high", "ER_hp", "ER_lp", "HER2_scmod1")
validation_all = validation_all %>% dplyr::select(-Basal)
validation_all[c("HER2", "LumB", "LumA",
                 "Normal", "T2", "Endo", "Cluster_2", "IC2", "IC3", "IC4",
                 "IC5", "IC6", "IC7", "IC8", "IC9", "IC10", "Mammaprint_risk_yes",
                 "rorS_risk_interm", "rorS_risk_high", "ER_hp", "ER_lp", "HER2_scmod1")] = 
  lapply(validation_all[c("HER2", "LumB", "LumA", "Normal", "T2", "Endo",
                          "Cluster_2", "IC2", "IC3", "IC4",
                          "IC5", "IC6", "IC7", "IC8", "IC9", "IC10", "Mammaprint_risk_yes",
                          "rorS_risk_interm", "rorS_risk_high", "ER_hp", "ER_lp", "HER2_scmod1")],
         function(x) factor(x, levels = c(0,1), 
                            labels = c(0,1)))

Validation = list(validation_gs, validation_cluster_meta, validation_gs_meta, validation_all)
names(Validation) = c("GS", "Cluster_meta", "GS_meta", "All_together")
rm(validation_gs, validation_cluster_meta, validation_gs_meta, validation_all)

# Test set with dummy variables to be used later:
test_all = test_set %>% 
  dplyr::select(-nr, -Sample.ID)
test_all = cbind(test_all[,unchanged], 
                 as.data.frame(model.matrix(~0 + test_all$pam50 +
                                              test_all$Timepoint_coded +
                                              test_all$Treatment +
                                              test_all$Cluster +
                                              test_all$IC10 +
                                              test_all$Mammaprint_risk +
                                              test_all$rorS_risk +
                                              test_all$scmod1)))

colnames(test_all) = c(unchanged, "Basal", "HER2", "LumB", "LumA",
                       "Normal", "T2", "Endo", "Cluster_2", "IC2", "IC3", "IC4",
                       "IC5", "IC6", "IC7", "IC8", "IC9", "IC10", "Mammaprint_risk_yes",
                       "rorS_risk_interm", "rorS_risk_high", "ER_hp", "ER_lp", "HER2_scmod1")
test_all = test_all %>% dplyr::select(-Basal)
test_all[c("HER2", "LumB", "LumA",
           "Normal", "T2", "Endo", "Cluster_2", "IC2", "IC3", "IC4",
           "IC5", "IC6", "IC7", "IC8", "IC9", "IC10", "Mammaprint_risk_yes",
           "rorS_risk_interm", "rorS_risk_high", "ER_hp", "ER_lp", "HER2_scmod1")] = 
  lapply(test_all[c("HER2", "LumB", "LumA", "Normal", "T2", "Endo",
                    "Cluster_2", "IC2", "IC3", "IC4",
                    "IC5", "IC6", "IC7", "IC8", "IC9", "IC10", "Mammaprint_risk_yes",
                    "rorS_risk_interm", "rorS_risk_high", "ER_hp", "ER_lp", "HER2_scmod1")],
         function(x) factor(x, levels = c(0,1), 
                            labels = c(0,1)))

# LogR #####
if(length(KBZ_DE_mapped$EntrezGene.ID[KBZ_DE_mapped$adj.P.Val < 0.05])>1){
  
  # Preparing output
  logrwb = createWorkbook()
  LogR_models_back = list()
  LogR_models_reg  = list()
  
  # Model tuning and building
  trControl = trainControl(method = "repeatedcv",
                           repeats = 10,
                           number = 10)
  
  for(i in 1:length(Training)){
    LogR_accuracy = data.frame(matrix(ncol=3,0))
    colnames(LogR_accuracy) = c("Method", "Accuracy_train", "Accuracy_validation")
    
    RNGversion("4.0.2")
    set.seed(123)
    LogR_model_back = train(Response ~., 
                            data = Training[[i]], 
                            method = "glmStepAIC", 
                            family = "binomial",
                            trControl = trControl,
                            na.action = "na.omit",
                            direction = "backward")
    
    train_pred = predict(LogR_model_back, Training[[i]])
    agree_train = ifelse(train_pred == Training[[i]]$Response, 1, 0)
    accuracy_train = sum(agree_train) / nrow(Training[[i]])
    
    validation_pred = predict(LogR_model_back, Validation[[i]])
    agree_validation = ifelse(validation_pred == Validation[[i]]$Response, 1, 0)
    accuracy_validation = sum(agree_validation) / nrow(Validation[[i]])
    
    LogR_models_back[[i]] = LogR_model_back
    LogR_accuracy = rbind(LogR_accuracy, c("Backward LogR", accuracy_train, accuracy_validation))
    rm(train_pred, agree_train, accuracy_train, validation_pred, 
       agree_validation, accuracy_validation)
    
    RNGversion("4.0.2")
    set.seed(123)
    LogR_model_reg = train(Response ~., 
                           data = Training[[i]], 
                           method = "glmnet", 
                           family = "binomial",
                           trControl = trControl,
                           na.action = "na.omit")
    
    train_pred = predict(LogR_model_reg, Training[[i]])
    agree_train = ifelse(train_pred == Training[[i]]$Response, 1, 0)
    accuracy_train = sum(agree_train) / nrow(Training[[i]])
    
    validation_pred = predict(LogR_model_reg, Validation[[i]])
    agree_validation = ifelse(validation_pred == Validation[[i]]$Response, 1, 0)
    accuracy_validation = sum(agree_validation) / nrow(Validation[[i]])
    
    LogR_models_reg[[i]] = LogR_model_reg
    LogR_accuracy = rbind(LogR_accuracy, c("Regularised LogR", accuracy_train, accuracy_validation))
    LogR_accuracy = LogR_accuracy[2:nrow(LogR_accuracy),]
    
    addWorksheet(logrwb, names(Training)[i])
    writeData(logrwb, names(Training)[i], LogR_accuracy)
    
    rm(train_pred, agree_train, accuracy_train, validation_pred, 
       agree_validation, accuracy_validation, LogR_model_reg, LogR_model_back)
  }
  
  saveWorkbook(logrwb, file = "data/Output sets/DGEA_GSEA_ML/kNN/LogR.xlsx",
               overwrite = TRUE); rm(logrwb, trControl)
  
  names(LogR_models_reg) = names(Training)
  names(LogR_models_back) = names(Training)
  LogR = list(LogR_models_reg, LogR_models_back)
  names(LogR) = c("Regularised models", "Backward models")
  rm(LogR_models_back, LogR_models_reg)
} 

# Decision trees #####
if(length(KBZ_DE_mapped$EntrezGene.ID[KBZ_DE_mapped$adj.P.Val < 0.05])>0){
  
  # Preparing output
  grid1 = expand.grid(model = c("tree", "rule"),
                      trials = c(1,5,10,15,20,25,30,35),
                      winnow = c(TRUE, FALSE))
  ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 10,
                      selectionFunction = "best")
  grid_c50 = expand.grid(model="tree", trials=c(10,15,25,50,100), winnow=T)
  ctrls = trainControl(method = "repeatedcv", number = 10, repeats = 10,
                       selectionFunction = "best", savePredictions = TRUE,
                       classProbs = TRUE, summaryFunction = twoClassSummary)
  grid2 = expand.grid(mtry=c(2,4,8,10,12,14,16,17, 25, 35, 45, 55))
  dtwb = createWorkbook()
  DT = list()
  
  for(i in 1:length(Training)){
    DT_Accuracy = data.frame(matrix(ncol=3,0))
    colnames(DT_Accuracy) = c("Model", "Accuracy_train", "Accuracy_validation")
    
    # C5.0 - k
    RNGversion("4.0.2")
    set.seed(123)
    C50_model = train(Response ~., data = Training[[i]],
                      method = "C5.0",
                      trControl = ctrl,
                      metric = "Kappa",
                      tuneGrid = grid1,
                      na.action = "na.omit")
    
    model_preds_validation = predict(C50_model, Validation[[i]])
    model_table_validation = table(model_preds_validation, Validation[[i]]$Response)
    model_accuracy_validation = (model_table_validation[1] + model_table_validation[4])/
      sum(model_table_validation)
    
    DT_Accuracy = rbind(DT_Accuracy, c("C5.0 - k", max(C50_model$results$Accuracy),
                                       model_accuracy_validation))
    
    rm(model_preds_validation, model_table_validation, model_accuracy_validation)
    
    # Bagging
    RNGversion("4.0.2")
    set.seed(123)
    model_bag = bagging(Response ~., data = Training[[i]],
                        nbagg = 100)
    
    bag_pred_train = predict(model_bag, Training[[i]])
    bag_train_t = table(bag_pred_train$class, Training[[i]]$Response)
    bagging_train_accuracy = (bag_train_t[1] + bag_train_t[4])/
      sum(bag_train_t)
    
    bag_pred_validation = predict(model_bag, Validation[[i]])
    bag_validation_t = table(bag_pred_validation$class, Validation[[i]]$Response)
    bagging_accuracy_validation = (bag_validation_t[1] + bag_validation_t[4])/
      sum(bag_validation_t)
    
    DT_Accuracy = rbind(DT_Accuracy, c("100X Bagging", bagging_train_accuracy, 
                                       bagging_accuracy_validation))
    
    rm(bag_pred_train, bag_train_t, bagging_train_accuracy,
       bag_pred_validation, bag_validation_t, bagging_accuracy_validation)
    
    # Boosting
    RNGversion("4.0.2")
    set.seed(123)
    boost = boosting(Response ~., data = Training[[i]])
    
    boost_train_pred = predict(boost, Training[[i]])
    cm_boost_train = boost_train_pred$confusion
    boost_accuracy_train = (cm_boost_train[1] + cm_boost_train[4])/
      sum(cm_boost_train)
    
    boost_validation_pred = predict(boost, Validation[[i]])
    cm_boost_validation = boost_validation_pred$confusion
    boost_accuracy_validation = (cm_boost_validation[1] + cm_boost_validation[4])/
      sum(cm_boost_validation)
    
    DT_Accuracy = rbind(DT_Accuracy, c("Boosting", boost_accuracy_train,
                                       boost_accuracy_validation))
    
    rm(boost_train_pred, cm_boost_train, boost_accuracy_train,
       boost_validation_pred, cm_boost_validation, boost_accuracy_validation)
    
    # Boosting with cross-validation on the whole dataset
    # train
    RNGversion("4.0.2")
    set.seed(123)
    adaboost_cv = boosting.cv(Response ~., data = cbind(rbind(Training[[i]],
                                                              Validation[[i]])),
                              mfinal = 100)
    
    cm_adaboost = adaboost_cv$confusion
    adaboost_accuracy = (cm_adaboost[1] + cm_adaboost[4])/sum(cm_adaboost)
    
    DT_Accuracy = rbind(DT_Accuracy, c("Adaboost", adaboost_accuracy, adaboost_accuracy))
    
    rm(cm_adaboost, adaboost_accuracy)
    
    # Random Forest
    RNGversion("4.0.2")
    set.seed(123)
    rf = randomForest(Response ~.,data = Training[[i]],
                      na.action = "na.omit")
    rf
    rf_accuracy_train = (rf$confusion[1] + rf$confusion[4])/
      sum(rf$confusion)
    
    rf_preds_validation = predict(rf, Validation[[i]])
    rf_table_validation = table(rf_preds_validation, Validation[[i]]$Response)
    rf_validation_accuracy = (rf_table_validation[1] + rf_table_validation[4])/
      sum(rf_table_validation)
    
    DT_Accuracy = rbind(DT_Accuracy, c("Random Forest - k", rf_accuracy_train,
                                       rf_validation_accuracy))
    
    rm(rf_accuracy_train, rf_preds_validation, rf_table_validation, rf_validation_accuracy)
    
    # RF with grid2 and ctrls, ROC
    RNGversion("4.0.2")
    set.seed(123)
    comp_rf = train(Response ~., data = Training[[i]],
                    method = "rf",
                    metric = "ROC", trControl = ctrls, tuneGrid = grid2,
                    na.action = "na.omit")
    
    comp_rf_preds_validation = predict(comp_rf, Validation[[i]])
    comp_rf_table_validation = table(comp_rf_preds_validation, Validation[[i]]$Response)
    comp_rf_accuracy_validation = (comp_rf_table_validation[1] + comp_rf_table_validation[4])/
      sum(comp_rf_table_validation)
    
    DT_Accuracy = rbind(DT_Accuracy, c("Random Forest - ROC", max(comp_rf$results$ROC),
                                       comp_rf_accuracy_validation))
    
    rm(comp_rf_preds_validation, comp_rf_table_validation, comp_rf_accuracy_validation)
    
    # C5.0 with grid_c50
    RNGversion("4.0.2")
    set.seed(123)
    m_c50 = train(Response ~., data = Training[[i]],
                  method="C5.0", metric = "ROC",
                  trControl = ctrls, tuneGrid = grid_c50, 
                  na.action = "na.omit")
    
    m_c50_preds_validation = predict(m_c50, Validation[[i]])
    m_c50_table_validation = table(m_c50_preds_validation, Validation[[i]]$Response)
    m_c50_accuracy_validation = (m_c50_table_validation[1] + m_c50_table_validation[4])/
      sum(m_c50_table_validation)
    
    DT_Accuracy = rbind(DT_Accuracy, c("C5.0 - ROC", max(m_c50$results$ROC),
                                       m_c50_accuracy_validation))
    
    rm(m_c50_preds_validation, m_c50_table_validation, m_c50_accuracy_validation)
    
    # Finalising
    DT_models = list(C50_model, model_bag, boost, adaboost_cv, 
                     rf, comp_rf, m_c50)
    names(DT_models) = c("C5.0 - k", "100X Bagging", "Boosting", "Boosting.cv",
                         "RForest - k", "Rforest - ROC", 
                         "C5.0 - ROC")
    DT[[i]] = DT_models
    DT_Accuracy = DT_Accuracy[2:nrow(DT_Accuracy),]
    addWorksheet(dtwb, names(Training)[i])
    writeData(dtwb, names(Training)[i], DT_Accuracy)
    rm(C50_model, model_bag, boost, adaboost_cv, rf, comp_rf, m_c50, 
       DT_Accuracy, DT_models) 
    gc()
    cat(paste0("Done with ", names(Training)[i], "\n"))
  }
  
  saveWorkbook(dtwb, file = "data/Output sets/DGEA_GSEA_ML/kNN/DT.xlsx",
               overwrite = TRUE); rm(dtwb)
  names(DT) = names(Training)
}

# SVMs #####
if(length(KBZ_DE_mapped$EntrezGene.ID[KBZ_DE_mapped$adj.P.Val < 0.05])>0){
  
  # Preparing output
  SVM = list()
  svmwb = createWorkbook()
  cost_values = seq(from = 0.5, to = 20, by = 0.5)
  linear_svm_ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 10,
                                 selectionFunction = "best")
  linear_svm_tune = expand.grid(C = cost_values)
  L2_linear_svm_ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 10,
                                    selectionFunction = "best")
  L2_linear_svm_tune = expand.grid(cost = cost_values,
                                   Loss = "L2")
  rbf_svm_ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 10,
                              selectionFunction = "best")
  rbf_svm_tune = expand.grid(C = cost_values,
                             sigma = c(0, 0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9))
  
  for(i in 1:length(Training)){
    SVM_Accuracy = data.frame(matrix(ncol=5,0))
    colnames(SVM_Accuracy) = c("Kernel", "Cost_value_C", "Sigma", 
                               "Accuracy_train", "Accuracy_validation")
    
    # Linear kernel - simple
    RNGversion("4.0.2")
    set.seed(123)
    linear_svm_model = train(Response ~., data = Training[[i]], 
                             method = "svmLinear", 
                             preProcess = NULL,
                             trControl = linear_svm_ctrl,
                             tuneGrid = linear_svm_tune,
                             na.action = "na.omit")
    
    train_pred = predict(linear_svm_model, Training[[i]])
    agree_train = ifelse(train_pred == Training[[i]]$Response, 1, 0)
    accuracy_train = sum(agree_train) / nrow(Training[[i]])
    
    validation_pred = predict(linear_svm_model, Validation[[i]])
    agree_validation = ifelse(validation_pred == Validation[[i]]$Response, 1, 0)
    accuracy_validation = sum(agree_validation) / nrow(Validation[[i]])
    
    SVM_Accuracy = rbind(SVM_Accuracy, c("Linear", linear_svm_model$bestTune[["C"]], "NA",
                                         accuracy_train, accuracy_validation))
    
    rm(train_pred, agree_train, accuracy_train, validation_pred, agree_validation,
       accuracy_validation)
    
    # Linear kernel - L2 regularised
    RNGversion("4.0.2")
    set.seed(123)
    L2_linear_svm_model = train(Response ~., data = Training[[i]], 
                                method = "svmLinear3", 
                                preProcess = NULL,
                                trControl = L2_linear_svm_ctrl,
                                tuneGrid = L2_linear_svm_tune,
                                na.action = "na.omit")
    
    train_pred = predict(L2_linear_svm_model, Training[[i]])
    agree_train = ifelse(train_pred == Training[[i]]$Response, 1, 0)
    accuracy_train = sum(agree_train) / nrow(Training[[i]])
    
    validation_pred = predict(L2_linear_svm_model, Validation[[i]])
    agree_validation = ifelse(validation_pred == Validation[[i]]$Response, 1, 0)
    accuracy_validation = sum(agree_validation) / nrow(Validation[[i]])
    
    SVM_Accuracy = rbind(SVM_Accuracy, c("L2_Linear", L2_linear_svm_model$bestTune[["cost"]], 
                                         "NA", accuracy_train, accuracy_validation))
    
    rm(train_pred, agree_train, accuracy_train, validation_pred, agree_validation,
       accuracy_validation)
    
    # Radial Basis Function
    RNGversion("4.0.2")
    set.seed(123)
    rbf_svm_model = train(Response ~., data = Training[[i]], 
                          method = "svmRadial", 
                          preProcess = NULL,
                          trControl = rbf_svm_ctrl,
                          tuneGrid = rbf_svm_tune,
                          na.action = "na.omit")
    
    
    train_pred = predict(rbf_svm_model, Training[[i]])
    agree_train = ifelse(train_pred == Training[[i]]$Response, 1, 0)
    accuracy_train = sum(agree_train) / nrow(Training[[i]])
    
    validation_pred = predict(rbf_svm_model, Validation[[i]])
    agree_validation = ifelse(validation_pred == Validation[[i]]$Response, 1, 0)
    accuracy_validation = sum(agree_validation) / nrow(Validation[[i]])
    
    SVM_Accuracy = rbind(SVM_Accuracy, c("RBF", rbf_svm_model$bestTune[["C"]], 
                                         rbf_svm_model$bestTune[["sigma"]],
                                         accuracy_train, accuracy_validation))
    
    rm(train_pred, agree_train, accuracy_train, validation_pred, agree_validation,
       accuracy_validation)
    
    # Saving
    SVM_Accuracy = SVM_Accuracy[2:nrow(SVM_Accuracy),]
    SVM_models = list(linear_svm_model, L2_linear_svm_model, rbf_svm_model)
    names(SVM_models) = c("Linear", "L2 Linear", "RBF")
    SVM[[i]] = SVM_models
    addWorksheet(svmwb, names(Training)[i])
    writeData(svmwb, names(Training)[i], SVM_Accuracy)
    rm(linear_svm_model, L2_linear_svm_model,
       rbf_svm_model, SVM_Accuracy, SVM_models); gc()
    cat(paste0("Done with ", names(Training)[i], "\n"))
  }
  
  saveWorkbook(svmwb, file = "data/Output sets/DGEA_GSEA_ML/kNN/SVM.xlsx",
               overwrite = TRUE); rm(svmwb)
  names(SVM) = names(Training)
}
