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

# Enable parallel programming
library(doParallel)
system <- Sys.info()['sysname']
cores <- makeCluster(detectCores(), type='PSOCK')
cl <- NULL
if (system == 'Windows') {
  cl <- makeCluster(getOption('cl.cores', cores))
  registerDoParallel(cl)
  registerDoSEQ()
  on.exit(stopCluster(cl))
} else {
  options('mc.cores' = cores)
  registerDoParallel(cores)
}
# A function that outputs performance measures
err_metric=function(CM)
{
  TN =CM[1,1]
  TP =CM[2,2]
  FP =CM[1,2]
  FN =CM[2,1]
  precision =(TP)/(TP+FP)
  recall_score =(FP)/(FP+TN)
  
  f1_score=2*((precision*recall_score)/(precision+recall_score))
  accuracy_model  =(TP+TN)/(TP+TN+FP+FN)
  False_positive_rate =(FP)/(FP+TN)
  False_negative_rate =(FN)/(FN+TP)
  
  print(paste("Precision value of the model: ",round(precision,2)))
  print(paste("Accuracy of the model: ",round(accuracy_model,2)))
  print(paste("Recall value of the model: ",round(recall_score,2)))
  print(paste("False Positive rate of the model: ",round(False_positive_rate,2)))
  
  print(paste("False Negative rate of the model: ",round(False_negative_rate,2)))
  
  print(paste("f1 score of the model: ",round(f1_score,2)))
}

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
global_expression = global_expression %>%
  dplyr::rename(ID = GEO_accession)

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
rownames(global_expression) = global_expression$ID
data = global_expression %>%
  dplyr::select(-ID, -Tissue_type, -Study, -Source)
data = as.data.frame(t(data))
RNGversion("4.0.2")
mock = M3C(mydata = data, des = global_expression, iters = 100, repsref = 250, 
           repsreal = 250, seed = 123, fsize = 20, lthick = 2, dotsize = 1.25)

# optimal K: 2
cluster_number = c(2,3,4,5,6,7,8,9,10)

chifit_type_p = list()
for (k in seq(2,10)){
  myresults <- mock$realdataresults[[k]]$ordered_annotation
  chifit_type <- suppressWarnings(chisq.test(table(myresults[c('consensuscluster','Tissue_type')])))
  chifit_type_p[[k-1]] = chifit_type$p.value
}

chifit_source_p = list()
for (k in seq(2,10)){
  myresults <- mock$realdataresults[[k]]$ordered_annotation
  chifit_source <- suppressWarnings(chisq.test(table(myresults[c('consensuscluster','Source')])))
  chifit_source_p[[k-1]] = chifit_source$p.value
}

cluster_corr = as.data.frame(cbind(cluster_number, chifit_type_p, chifit_source_p))
colnames(cluster_corr) = c("No. of clusters", "Tissue_type", "Source")

# PCA clustering
png("ML/Cluster_PCA.png",
    width = 1920, height = 1080)
pca(mock$realdataresults[[2]]$ordered_data,
    labels = mock$realdataresults[[2]]$ordered_annotation$consensuscluster,
    legendtextsize = 20, axistextsize = 20, dotsize = 5)
dev.off()

# Consensus index plot
png("ML/Consensus_index.png",
    width = 1920, height = 1080)
mock[["plots"]][[1]]
dev.off()

# Entropy plot
png("ML/Entropy.png",
    width = 1920, height = 1080)
mock[["plots"]][[2]]
dev.off()

# Statistical significance of clusters
png("ML/Stat_Sig.png",
    width = 1920, height = 1080)
mock[["plots"]][[3]]
dev.off()

# RCSI plot
png("ML/RCSI.png",
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
saveWorkbook(kmwb, file = "ML/Clustering.xlsx",
             overwrite = TRUE); rm(kmwb); gc()

rm(data, res); gc()

cluster_results = data.frame(names(mock[["realdataresults"]][[2]][["assignments"]]), 
                             mock[["realdataresults"]][[2]][["assignments"]])
# 2: optimal number of clusters
colnames(cluster_results) = c("ID", "Cluster")

# Training, Validation, Test set splitting #####
# The sets are generated once and then filtered for variables of interest in
# the ML training, validating, testing process
data = global_expression %>%
  inner_join(cluster_results, by = "ID") %>%
  dplyr::select(ID, Tissue_type, Study, Source, Cluster, everything())
data$Cluster = as.factor(data$Cluster)

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
  group_by(Tissue_type, Source, Study, Cluster) %>%
  sample_frac(0.65) %>%
  as.data.frame()
RNGversion("4.0.2")
set.seed(123)
validation_set = anti_join(tobesplit, train_set) %>%
  group_by(Tissue_type, Source, Study, Cluster) %>%
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
  dplyr::select(-nr, -Cluster, -ID, -Source, -Study)
unchanged = colnames(train_gs)

train_cluster_meta = train_set %>% 
  dplyr::select(Cluster, Source, Tissue_type)
train_cluster_meta = cbind(train_cluster_meta[,"Tissue_type"], as.data.frame(model.matrix(~0 +
                                                                                         train_cluster_meta$Cluster +
                                                                                         train_cluster_meta$Source)))
colnames(train_cluster_meta) = c("Tissue_type", "Cluster_1", "Cluster_2", "Blood")
train_cluster_meta = train_cluster_meta %>% dplyr::select(-Cluster_1)
train_cluster_meta[colnames(train_cluster_meta[,2:ncol(train_cluster_meta)])] = 
  lapply(train_cluster_meta[colnames(train_cluster_meta[,2:ncol(train_cluster_meta)])], 
         function(x) factor(x, levels = c(0,1), labels = c(0,1)))

train_gs_meta = train_set %>% 
  dplyr::select(-nr, -ID, -Cluster, -Study)
train_gs_meta = cbind(train_gs_meta[,unchanged], as.data.frame(model.matrix(~0 +
                                                                              train_gs_meta$Source)))
colnames(train_gs_meta) = c(unchanged, "Tumor", "Blood")
train_gs_meta = train_gs_meta %>% dplyr::select(-Tumor)
train_gs_meta["Blood"] = 
  lapply(train_gs_meta[c("Blood")],
         function(x) factor(x, levels = c(0,1), 
                            labels = c(0,1)))

train_all = train_set %>% 
  dplyr::select(-nr, -ID, -Study)
train_all = cbind(train_all[,unchanged], 
                  as.data.frame(model.matrix(~0 + train_all$Cluster +
                                               train_all$Source)))

colnames(train_all) = c(unchanged, "Cluster_1", "Cluster_2", "Blood")
train_all = train_all %>% dplyr::select(-Cluster_1)
train_all[c("Cluster_2", "Blood")] = 
  lapply(train_all[c("Cluster_2", "Blood")],
         function(x) factor(x, levels = c(0,1), 
                            labels = c(0,1)))

Training = list(train_gs, train_cluster_meta, train_gs_meta, train_all)
names(Training) = c("GS", "Cluster_meta", "GS_meta", "All_together")
rm(train_gs, train_cluster_meta, train_gs_meta, train_all)

# Validation
validation_gs = validation_set %>% 
  dplyr::select(-nr, -Cluster, -ID, -Source, -Study)
unchanged = colnames(validation_gs)

validation_cluster_meta = validation_set %>% 
  dplyr::select(Cluster, Source, Tissue_type)
validation_cluster_meta = cbind(validation_cluster_meta[,"Tissue_type"], as.data.frame(model.matrix(~0 +
                                                                                                      validation_cluster_meta$Cluster +
                                                                                                      validation_cluster_meta$Source)))
colnames(validation_cluster_meta) = c("Tissue_type", "Cluster_1", "Cluster_2", "Blood")
validation_cluster_meta = validation_cluster_meta %>% dplyr::select(-Cluster_1)
validation_cluster_meta[colnames(validation_cluster_meta[,2:ncol(validation_cluster_meta)])] = 
  lapply(validation_cluster_meta[colnames(validation_cluster_meta[,2:ncol(validation_cluster_meta)])], 
         function(x) factor(x, levels = c(0,1), labels = c(0,1)))

validation_gs_meta = validation_set %>% 
  dplyr::select(-nr, -ID, -Cluster, -Study)
validation_gs_meta = cbind(validation_gs_meta[,unchanged], as.data.frame(model.matrix(~0 +
                                                                                        validation_gs_meta$Source)))
colnames(validation_gs_meta) = c(unchanged, "Tumor", "Blood")
validation_gs_meta = validation_gs_meta %>% dplyr::select(-Tumor)
validation_gs_meta["Blood"] = 
  lapply(validation_gs_meta[c("Blood")],
         function(x) factor(x, levels = c(0,1), 
                            labels = c(0,1)))

validation_all = validation_set %>% 
  dplyr::select(-nr, -ID, -Study)
validation_all = cbind(validation_all[,unchanged], 
                       as.data.frame(model.matrix(~0 + validation_all$Cluster +
                                                    validation_all$Source)))

colnames(validation_all) = c(unchanged, "Cluster_1", "Cluster_2", "Blood")
validation_all = validation_all %>% dplyr::select(-Cluster_1)
validation_all[c("Cluster_2", "Blood")] = 
  lapply(validation_all[c("Cluster_2", "Blood")],
         function(x) factor(x, levels = c(0,1), 
                            labels = c(0,1)))

Validation = list(validation_gs, validation_cluster_meta, validation_gs_meta, validation_all)
names(Validation) = c("GS", "Cluster_meta", "GS_meta", "All_together")
rm(validation_gs, validation_cluster_meta, validation_gs_meta, validation_all)

# Test set with dummy variables to be used later:
test_all = test_set %>% 
  dplyr::select(-nr, -ID, -Study)
test_all = cbind(test_all[,unchanged], 
                 as.data.frame(model.matrix(~0 + test_all$Cluster +
                                              test_all$Source)))

colnames(test_all) = c(unchanged, "Cluster_1", "Cluster_2", "Blood")
test_all = test_all %>% dplyr::select(-Cluster_1)
test_all[c("Cluster_2", "Blood")] = 
  lapply(test_all[c("Cluster_2", "Blood")],
         function(x) factor(x, levels = c(0,1), 
                            labels = c(0,1)))

# LogR #####
# Preparing output
logrwb = createWorkbook()
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
  LogR_model_reg = train(Tissue_type ~., 
                         data = Training[[i]], 
                         method = "glmnet", 
                         family = "binomial",
                         trControl = trControl,
                         na.action = "na.omit")
  
  train_pred = predict(LogR_model_reg, Training[[i]])
  agree_train = ifelse(train_pred == Training[[i]]$Tissue_type, 1, 0)
  accuracy_train = sum(agree_train) / nrow(Training[[i]])
  
  validation_pred = predict(LogR_model_reg, Validation[[i]])
  agree_validation = ifelse(validation_pred == Validation[[i]]$Tissue_type, 1, 0)
  accuracy_validation = sum(agree_validation) / nrow(Validation[[i]])
  
  LogR_models_reg[[i]] = LogR_model_reg
  LogR_accuracy = rbind(LogR_accuracy, c("Regularised LogR", accuracy_train, accuracy_validation))
  LogR_accuracy = LogR_accuracy[2:nrow(LogR_accuracy),]
  
  addWorksheet(logrwb, names(Training)[i])
  writeData(logrwb, names(Training)[i], LogR_accuracy)
  
  rm(train_pred, agree_train, accuracy_train, validation_pred, 
     agree_validation, accuracy_validation, LogR_model_reg)
}

saveWorkbook(logrwb, file = "ML/LogR.xlsx",
             overwrite = TRUE); rm(logrwb, trControl)

names(LogR_models_reg) = names(Training)

# Decision trees #####
# Preparing output
grid1 = expand.grid(model = c("tree", "rule"),
                    trials = c(10,15,25,50,100),
                    winnow = c(TRUE, FALSE))
ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 10,
                    selectionFunction = "best", allowParallel = TRUE)
grid_c50 = expand.grid(model = c("tree", "rule"),
                       trials = c(10,15,25,50,100),
                       winnow = c(TRUE, FALSE))
ctrls = trainControl(method = "repeatedcv", number = 10, repeats = 10,
                     selectionFunction = "best", savePredictions = TRUE,
                     classProbs = TRUE, summaryFunction = twoClassSummary,
                     allowParallel = TRUE)
grid2 = expand.grid(mtry=c(2, 10, 15, 25, 50, 100, 200, 300, 400, 500))
dtwb = createWorkbook()
DT = list()

for(i in 1:length(Training)){
  DT_Accuracy = data.frame(matrix(ncol=3,0))
  colnames(DT_Accuracy) = c("Model", "Accuracy_train", "Accuracy_validation")
  
  # C5.0 - k
  RNGversion("4.0.2")
  set.seed(123)
  C50_model = train(Tissue_type ~., data = Training[[i]],
                    method = "C5.0",
                    trControl = ctrl,
                    metric = "Kappa",
                    tuneGrid = grid1,
                    na.action = "na.omit")
  
  model_preds_validation = predict(C50_model, Validation[[i]])
  model_table_validation = table(model_preds_validation, Validation[[i]]$Tissue_type)
  model_accuracy_validation = (model_table_validation[1] + model_table_validation[4])/
    sum(model_table_validation)
  
  DT_Accuracy = rbind(DT_Accuracy, c("C5.0 - k", max(C50_model$results$Accuracy),
                                     model_accuracy_validation))
  
  rm(model_preds_validation, model_table_validation, model_accuracy_validation)
  
  # Bagging
  RNGversion("4.0.2")
  set.seed(123)
  model_bag = bagging(Tissue_type ~., data = Training[[i]],
                      nbagg = 100)
  
  bag_pred_train = predict(model_bag, Training[[i]])
  bag_train_t = table(bag_pred_train$class, Training[[i]]$Tissue_type)
  bagging_train_accuracy = (bag_train_t[1] + bag_train_t[4])/
    sum(bag_train_t)
  
  bag_pred_validation = predict(model_bag, Validation[[i]])
  bag_validation_t = table(bag_pred_validation$class, Validation[[i]]$Tissue_type)
  bagging_accuracy_validation = (bag_validation_t[1] + bag_validation_t[4])/
    sum(bag_validation_t)
  
  DT_Accuracy = rbind(DT_Accuracy, c("100X Bagging", bagging_train_accuracy, 
                                     bagging_accuracy_validation))
  
  rm(bag_pred_train, bag_train_t, bagging_train_accuracy,
     bag_pred_validation, bag_validation_t, bagging_accuracy_validation)
  
  # Boosting
  RNGversion("4.0.2")
  set.seed(123)
  boost = boosting(Tissue_type ~., data = Training[[i]])
  
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
  adaboost_cv = boosting.cv(Tissue_type ~., data = cbind(rbind(Training[[i]],
                                                               Validation[[i]])),
                            mfinal = 100)
  
  cm_adaboost = adaboost_cv$confusion
  adaboost_accuracy = (cm_adaboost[1] + cm_adaboost[4])/sum(cm_adaboost)
  
  DT_Accuracy = rbind(DT_Accuracy, c("Adaboost", adaboost_accuracy, adaboost_accuracy))
  
  rm(cm_adaboost, adaboost_accuracy)
  
  # Random Forest
  RNGversion("4.0.2")
  set.seed(123)
  rf = train(Tissue_type ~., data = Training[[i]],
             method = "rf",
             metric = "Kappa", trControl = ctrl, tuneGrid = grid2,
             na.action = "na.omit")
  
  rf_preds_validation = predict(rf, Validation[[i]])
  rf_table_validation = table(rf_preds_validation, Validation[[i]]$Tissue_type)
  rf_validation_accuracy = (rf_table_validation[1] + rf_table_validation[4])/
    sum(rf_table_validation)
  
  DT_Accuracy = rbind(DT_Accuracy, c("Random Forest - k", max(rf$results$Accuracy),
                                     rf_validation_accuracy))
  
  rm(rf_preds_validation, rf_table_validation, rf_validation_accuracy)
  
  # RF with grid2 and ctrls, ROC
  RNGversion("4.0.2")
  set.seed(123)
  comp_rf = train(Tissue_type ~., data = Training[[i]],
                  method = "rf",
                  metric = "ROC", trControl = ctrls, tuneGrid = grid2,
                  na.action = "na.omit")
  
  comp_rf_preds_validation = predict(comp_rf, Validation[[i]])
  comp_rf_table_validation = table(comp_rf_preds_validation, Validation[[i]]$Tissue_type)
  comp_rf_accuracy_validation = (comp_rf_table_validation[1] + comp_rf_table_validation[4])/
    sum(comp_rf_table_validation)
  
  DT_Accuracy = rbind(DT_Accuracy, c("Random Forest - ROC", max(comp_rf$results$ROC),
                                     comp_rf_accuracy_validation))
  
  rm(comp_rf_preds_validation, comp_rf_table_validation, comp_rf_accuracy_validation)
  
  # C5.0 with grid_c50
  RNGversion("4.0.2")
  set.seed(123)
  m_c50 = train(Tissue_type ~., data = Training[[i]],
                method="C5.0", metric = "ROC",
                trControl = ctrls, tuneGrid = grid_c50, 
                na.action = "na.omit")
  
  m_c50_preds_validation = predict(m_c50, Validation[[i]])
  m_c50_table_validation = table(m_c50_preds_validation, Validation[[i]]$Tissue_type)
  m_c50_accuracy_validation = (m_c50_table_validation[1] + m_c50_table_validation[4])/
    sum(m_c50_table_validation)
  
  DT_Accuracy = rbind(DT_Accuracy, c("C5.0 - ROC", max(m_c50$results$ROC),
                                     m_c50_accuracy_validation))
  
  rm(m_c50_preds_validation, m_c50_table_validation, m_c50_accuracy_validation)
  
  # Finalising
  DT_models = list(C50_model, model_bag, boost, adaboost_cv, 
                   rf, comp_rf, m_c50)
  names(DT_models) = c("C5.0 - k", "100X Bagging", "Boosting", "Boosting.cv",
                       "RForest - k", "RForest - ROC", 
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

saveWorkbook(dtwb, file = "ML/DT.xlsx",
             overwrite = TRUE); rm(dtwb)
names(DT) = names(Training)

# SVMs #####
# Preparing output
SVM = list()
svmwb = createWorkbook()
cost_values = seq(from = 0.5, to = 20, by = 0.5)
linear_svm_ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 10,
                               selectionFunction = "best", allowParallel = TRUE)
linear_svm_tune = expand.grid(C = cost_values)
L2_linear_svm_ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 10,
                                  selectionFunction = "best", allowParallel = TRUE)
L2_linear_svm_tune = expand.grid(cost = cost_values,
                                 Loss = "L2")
rbf_svm_ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 10,
                            selectionFunction = "best", allowParallel = TRUE)
rbf_svm_tune = expand.grid(C = cost_values,
                           sigma = c(0, 0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9))

for(i in 1:length(Training)){
  SVM_Accuracy = data.frame(matrix(ncol=5,0))
  colnames(SVM_Accuracy) = c("Kernel", "Cost_value_C", "Sigma", 
                             "Accuracy_train", "Accuracy_validation")
  
  # Linear kernel - simple
  RNGversion("4.0.2")
  set.seed(123)
  linear_svm_model = train(Tissue_type ~., data = Training[[i]], 
                           method = "svmLinear", 
                           preProcess = NULL,
                           trControl = linear_svm_ctrl,
                           tuneGrid = linear_svm_tune,
                           na.action = "na.omit")
  
  train_pred = predict(linear_svm_model, Training[[i]])
  agree_train = ifelse(train_pred == Training[[i]]$Tissue_type, 1, 0)
  accuracy_train = sum(agree_train) / nrow(Training[[i]])
  
  validation_pred = predict(linear_svm_model, Validation[[i]])
  agree_validation = ifelse(validation_pred == Validation[[i]]$Tissue_type, 1, 0)
  accuracy_validation = sum(agree_validation) / nrow(Validation[[i]])
  
  SVM_Accuracy = rbind(SVM_Accuracy, c("Linear", linear_svm_model$bestTune[["C"]], "NA",
                                       accuracy_train, accuracy_validation))
  
  rm(train_pred, agree_train, accuracy_train, validation_pred, agree_validation,
     accuracy_validation)
  
  # Linear kernel - L2 regularised
  RNGversion("4.0.2")
  set.seed(123)
  L2_linear_svm_model = train(Tissue_type ~., data = Training[[i]], 
                              method = "svmLinear3", 
                              preProcess = NULL,
                              trControl = L2_linear_svm_ctrl,
                              tuneGrid = L2_linear_svm_tune,
                              na.action = "na.omit")
  
  train_pred = predict(L2_linear_svm_model, Training[[i]])
  agree_train = ifelse(train_pred == Training[[i]]$Tissue_type, 1, 0)
  accuracy_train = sum(agree_train) / nrow(Training[[i]])
  
  validation_pred = predict(L2_linear_svm_model, Validation[[i]])
  agree_validation = ifelse(validation_pred == Validation[[i]]$Tissue_type, 1, 0)
  accuracy_validation = sum(agree_validation) / nrow(Validation[[i]])
  
  SVM_Accuracy = rbind(SVM_Accuracy, c("L2_Linear", L2_linear_svm_model$bestTune[["cost"]], 
                                       "NA", accuracy_train, accuracy_validation))
  
  rm(train_pred, agree_train, accuracy_train, validation_pred, agree_validation,
     accuracy_validation)
  
  # Radial Basis Function
  RNGversion("4.0.2")
  set.seed(123)
  rbf_svm_model = train(Tissue_type ~., data = Training[[i]], 
                        method = "svmRadial", 
                        preProcess = NULL,
                        trControl = rbf_svm_ctrl,
                        tuneGrid = rbf_svm_tune,
                        na.action = "na.omit")
  
  
  train_pred = predict(rbf_svm_model, Training[[i]])
  agree_train = ifelse(train_pred == Training[[i]]$Tissue_type, 1, 0)
  accuracy_train = sum(agree_train) / nrow(Training[[i]])
  
  validation_pred = predict(rbf_svm_model, Validation[[i]])
  agree_validation = ifelse(validation_pred == Validation[[i]]$Tissue_type, 1, 0)
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

saveWorkbook(svmwb, file = "ML/SVM.xlsx",
             overwrite = TRUE); rm(svmwb)
names(SVM) = names(Training)

# Evaluating the best model on the test set #####
final_model_test_pred = predict(DT[["GS"]][["RForest - k"]], test_all)
cm_final_model_test = table(final_model_test_pred, test_all$Tissue_type)
final_model_test_accuracy = (cm_final_model_test[1] + cm_final_model_test[4])/
  sum(cm_final_model_test)

err_metric(cm_final_model_test)

# "Precision value of the model:  0.93"
# "Accuracy of the model:  0.85"
# "Recall value of the model:  0.2"
# "False Positive rate of the model:  0.2"
# "False Negative rate of the model:  0.14"
# "f1 score of the model:  0.33"