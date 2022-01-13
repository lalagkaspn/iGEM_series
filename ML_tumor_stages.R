# Machine Learning on Prostate Cancer samples: tumor stages
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

# Importing the tumor expression data frame
tumor_frame = read.xlsx("Tumor_samples_z_expression_matrix.xlsx") %>%
  dplyr::select(GEO_accession, Tissue_type, AJCC_classification, Study, everything()) %>%
  dplyr::filter(Tissue_type == "tumor") %>%
  dplyr::filter(is.na(AJCC_classification) == FALSE)

# Binding into a full expression frame
tumor_frame$Stage = ifelse(tumor_frame$AJCC_classification == 3 |
                             tumor_frame$AJCC_classification == 4, "Late", "Early")
tumor_frame = tumor_frame %>%
  dplyr::select(GEO_accession, Tissue_type, AJCC_classification, Stage, Study, everything()) %>%
  dplyr::rename(ID = GEO_accession)
tumor_frame[colnames(tumor_frame[,6:ncol(tumor_frame)])] = 
  lapply(tumor_frame[colnames(tumor_frame[,6:ncol(tumor_frame)])], 
         function(x) as.numeric(x))
tumor_frame[colnames(tumor_frame[,2:5])] = 
  lapply(tumor_frame[colnames(tumor_frame[,2:5])], 
         function(x) as.factor(x))

cols = colnames(tumor_frame)
cols_new = gsub("-", "_", cols)
colnames(tumor_frame) = cols_new
rm(cols, cols_new)

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

testNames <- colnames(tumor_frame)
res = t(sapply(testNames, testValidity)) # All are ok

# Consensus clustering #####
# Preparing data
library(M3C)
rownames(tumor_frame) = tumor_frame$ID
data = tumor_frame %>%
  dplyr::select(-ID, -Tissue_type, -Study, -AJCC_classification, -Stage)
data = as.data.frame(t(data))
RNGversion("4.0.2")
mock = M3C(mydata = data, des = tumor_frame, iters = 100, repsref = 250, 
           repsreal = 250, seed = 123, fsize = 20, lthick = 2, dotsize = 1.25)

# optimal K: 2
cluster_number = c(2,3,4,5,6,7,8,9,10)

chifit_stage_p = list()
for (k in seq(2,10)){
  myresults <- mock$realdataresults[[k]]$ordered_annotation
  chifit_stage <- suppressWarnings(chisq.test(table(myresults[c('consensuscluster','Stage')])))
  chifit_stage_p[[k-1]] = chifit_stage$p.value
}

chifit_AJCC_p = list()
for (k in seq(2,10)){
  myresults <- mock$realdataresults[[k]]$ordered_annotation
  chifit_AJCC <- suppressWarnings(chisq.test(table(myresults[c('consensuscluster',
                                                               'AJCC_classification')])))
  chifit_AJCC_p[[k-1]] = chifit_AJCC$p.value
}

cluster_corr = as.data.frame(cbind(cluster_number, chifit_stage_p, chifit_AJCC_p))
colnames(cluster_corr) = c("No. of clusters", "Tissue_type", "AJCC_classification")

# PCA clustering
png("ML/Tumor_stage/Cluster_PCA.png",
    width = 1920, height = 1080)
pca(mock$realdataresults[[2]]$ordered_data,
    labels = mock$realdataresults[[2]]$ordered_annotation$consensuscluster,
    legendtextsize = 20, axistextsize = 20, dotsize = 5)
dev.off()

# Consensus index plot
png("ML/Tumor_stage/Consensus_index.png",
    width = 1920, height = 1080)
mock[["plots"]][[1]]
dev.off()

# Entropy plot
png("ML/Tumor_stage/Entropy.png",
    width = 1920, height = 1080)
mock[["plots"]][[2]]
dev.off()

# Statistical significance of clusters
png("ML/Tumor_stage/Stat_Sig.png",
    width = 1920, height = 1080)
mock[["plots"]][[3]]
dev.off()

# RCSI plot
png("ML/Tumor_stage/RCSI.png",
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
saveWorkbook(kmwb, file = "ML/Tumor_stage/Clustering.xlsx",
             overwrite = TRUE); rm(kmwb); gc()

rm(data, res); gc()

cluster_results = data.frame(names(mock[["realdataresults"]][[2]][["assignments"]]), 
                             mock[["realdataresults"]][[2]][["assignments"]])
# 2: optimal number of clusters
colnames(cluster_results) = c("ID", "Cluster")

# Results are not significant, indicating that there are not clusters of patients
# in the data, i.e. the samples are quite similar

# Training, Validation, Test set splitting #####
# The sets are generated once and then filtered for variables of interest in
# the ML training, validating, testing process
data = tumor_frame %>%
  dplyr::select(ID, Tissue_type, Study, AJCC_classification, Stage, everything())

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
  group_by(Stage, Study) %>%
  sample_frac(0.6) %>%
  as.data.frame()
RNGversion("4.0.2")
set.seed(123)
validation_set = anti_join(tobesplit, train_set) %>%
  group_by(Stage, Study) %>%
  sample_frac(0.6) %>%
  as.data.frame()
test_set = anti_join(tobesplit, as.data.frame(rbind(train_set, validation_set)))

train_set = train_set %>% dplyr::ungroup() %>%
  dplyr::select(-nr, -ID, -Tissue_type, -Study, -AJCC_classification)
validation_set = validation_set %>% dplyr::ungroup() %>%
  dplyr::select(-nr, -ID, -Tissue_type, -Study, -AJCC_classification)
test_set = test_set %>% dplyr::ungroup() %>%
  dplyr::select(-nr, -ID, -Tissue_type, -Study, -AJCC_classification)

# LogR #####
# Preparing output
logrwb = createWorkbook()

# Model tuning and building
trControl = trainControl(method = "repeatedcv",
                         repeats = 10,
                         number = 10, allowParallel = TRUE)

LogR_accuracy = data.frame(matrix(ncol=3,0))
colnames(LogR_accuracy) = c("Method", "Accuracy_train", "Accuracy_validation")

RNGversion("4.0.2")
set.seed(123)
LogR_model_reg = train(Stage ~., 
                       data = train_set, 
                       method = "glmnet", 
                       family = "binomial",
                       trControl = trControl,
                       na.action = "na.omit")

train_pred = predict(LogR_model_reg, train_set)
agree_train = ifelse(train_pred == train_set$Stage, 1, 0)
accuracy_train = sum(agree_train) / nrow(train_set)

validation_pred = predict(LogR_model_reg, validation_set)
agree_validation = ifelse(validation_pred == validation_set$Stage, 1, 0)
accuracy_validation = sum(agree_validation) / nrow(validation_set)

LogR_accuracy = rbind(LogR_accuracy, c("Regularised LogR", accuracy_train, accuracy_validation))
LogR_accuracy = LogR_accuracy[2:nrow(LogR_accuracy),]

addWorksheet(logrwb, "LogR_performance")
writeData(logrwb, "LogR_performance", LogR_accuracy)

rm(train_pred, agree_train, accuracy_train, validation_pred, 
   agree_validation, accuracy_validation)

saveWorkbook(logrwb, file = "ML/Tumor_stage/LogR.xlsx",
             overwrite = TRUE); rm(logrwb, trControl)

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

DT_Accuracy = data.frame(matrix(ncol=3,0))
colnames(DT_Accuracy) = c("Model", "Accuracy_train", "Accuracy_validation")

# C5.0 - k
RNGversion("4.0.2")
set.seed(123)
C50_model = train(Stage ~., data = train_set,
                  method = "C5.0",
                  trControl = ctrl,
                  metric = "Kappa",
                  tuneGrid = grid1,
                  na.action = "na.omit")

model_preds_validation = predict(C50_model, validation_set)
model_table_validation = table(model_preds_validation, validation_set$Stage)
model_accuracy_validation = (model_table_validation[1] + model_table_validation[4])/
  sum(model_table_validation)

DT_Accuracy = rbind(DT_Accuracy, c("C5.0 - k", max(C50_model$results$Accuracy),
                                   model_accuracy_validation))

rm(model_preds_validation, model_table_validation, model_accuracy_validation)

# Bagging
RNGversion("4.0.2")
set.seed(123)
model_bag = bagging(Stage ~., data = train_set,
                    nbagg = 100)

bag_pred_train = predict(model_bag, train_set)
bag_train_t = table(bag_pred_train$class, train_set$Stage)
bagging_train_accuracy = (bag_train_t[1] + bag_train_t[4])/
  sum(bag_train_t)

bag_pred_validation = predict(model_bag, validation_set)
bag_validation_t = table(bag_pred_validation$class, validation_set$Stage)
bagging_accuracy_validation = (bag_validation_t[1] + bag_validation_t[4])/
  sum(bag_validation_t)

DT_Accuracy = rbind(DT_Accuracy, c("100X Bagging", bagging_train_accuracy, 
                                   bagging_accuracy_validation))

rm(bag_pred_train, bag_train_t, bagging_train_accuracy,
   bag_pred_validation, bag_validation_t, bagging_accuracy_validation)

# Boosting
RNGversion("4.0.2")
set.seed(123)
boost = boosting(Stage ~., data = train_set)

boost_train_pred = predict(boost, train_set)
cm_boost_train = boost_train_pred$confusion
boost_accuracy_train = (cm_boost_train[1] + cm_boost_train[4])/
  sum(cm_boost_train)

boost_validation_pred = predict(boost, validation_set)
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
adaboost_cv = boosting.cv(Stage ~., data = cbind(rbind(train_set,
                                                       validation_set)),
                          mfinal = 100)

cm_adaboost = adaboost_cv$confusion
adaboost_accuracy = (cm_adaboost[1] + cm_adaboost[4])/sum(cm_adaboost)

DT_Accuracy = rbind(DT_Accuracy, c("Adaboost", adaboost_accuracy, adaboost_accuracy))

rm(cm_adaboost, adaboost_accuracy)

# Random Forest
RNGversion("4.0.2")
set.seed(123)
rf = train(Stage ~., data = train_set,
           method = "rf",
           metric = "Kappa", trControl = ctrl, tuneGrid = grid2,
           na.action = "na.omit")

rf_preds_validation = predict(rf, validation_set)
rf_table_validation = table(rf_preds_validation, validation_set$Stage)
rf_validation_accuracy = (rf_table_validation[1] + rf_table_validation[4])/
  sum(rf_table_validation)

DT_Accuracy = rbind(DT_Accuracy, c("Random Forest - k", max(rf$results$Accuracy),
                                   rf_validation_accuracy))

rm(rf_preds_validation, rf_table_validation, rf_validation_accuracy)

# RF with grid2 and ctrls, ROC
RNGversion("4.0.2")
set.seed(123)
comp_rf = train(Stage ~., data = train_set,
                method = "rf",
                metric = "ROC", trControl = ctrls, tuneGrid = grid2,
                na.action = "na.omit")

comp_rf_preds_validation = predict(comp_rf, validation_set)
comp_rf_table_validation = table(comp_rf_preds_validation, validation_set$Stage)
comp_rf_accuracy_validation = (comp_rf_table_validation[1] + comp_rf_table_validation[4])/
  sum(comp_rf_table_validation)

DT_Accuracy = rbind(DT_Accuracy, c("Random Forest - ROC", max(comp_rf$results$ROC),
                                   comp_rf_accuracy_validation))

rm(comp_rf_preds_validation, comp_rf_table_validation, comp_rf_accuracy_validation)

# C5.0 with grid_c50
RNGversion("4.0.2")
set.seed(123)
m_c50 = train(Stage ~., data = train_set,
              method="C5.0", metric = "ROC",
              trControl = ctrls, tuneGrid = grid_c50, 
              na.action = "na.omit")

m_c50_preds_validation = predict(m_c50, validation_set)
m_c50_table_validation = table(m_c50_preds_validation, validation_set$Stage)
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

DT_Accuracy = DT_Accuracy[2:nrow(DT_Accuracy),]
addWorksheet(dtwb, "DT_performance")
writeData(dtwb, "DT_performance", DT_Accuracy)
rm(C50_model, model_bag, boost, adaboost_cv, rf, comp_rf, m_c50, 
   DT_Accuracy, DT_models) 
saveWorkbook(dtwb, file = "ML/Tumor_stage/DT.xlsx",
             overwrite = TRUE); rm(dtwb)

# SVMs #####
# Preparing output
svmwb = createWorkbook()
cost_values = seq(from = 0.5, to = 20, by = 0.5)
weights = c(0.15, 0.2, 0.25, 0.5, 0.75, 0.8, 0.85, 0.9, 0.95)
names(weights) = rep("Late", 9)
linear_svm_ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 10,
                               selectionFunction = "best", allowParallel = TRUE)
linear_svm_tune = expand.grid(C = cost_values)
L2_linear_svm_ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 10,
                                  selectionFunction = "best", allowParallel = TRUE)
L2_linear_svm_tune = expand.grid(cost = cost_values,
                                 Loss = "L2")
L2_linear_svm_tune_weights = expand.grid(cost = cost_values,
                                 Loss = "L2", weight = weights)
rbf_svm_ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 10,
                            selectionFunction = "best", allowParallel = TRUE)
rbf_svm_tune = expand.grid(C = cost_values,
                           sigma = c(0, 0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9))
rbf_svm_tune_weights = expand.grid(C = cost_values,
                           sigma = c(0, 0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9),
                           Weight = weights)

SVM_Accuracy = data.frame(matrix(ncol=6,0))
colnames(SVM_Accuracy) = c("Kernel", "Cost_value_C", "Sigma", "Weight",
                           "Accuracy_train", "Accuracy_validation")

# Linear kernel - simple
RNGversion("4.0.2")
set.seed(123)
linear_svm_model = train(Stage ~., data = train_set, 
                         method = "svmLinear", 
                         preProcess = NULL,
                         trControl = linear_svm_ctrl,
                         tuneGrid = linear_svm_tune,
                         na.action = "na.omit")

train_pred = predict(linear_svm_model, train_set)
agree_train = ifelse(train_pred == train_set$Stage, 1, 0)
accuracy_train = sum(agree_train) / nrow(train_set)

validation_pred = predict(linear_svm_model, validation_set)
agree_validation = ifelse(validation_pred == validation_set$Stage, 1, 0)
accuracy_validation = sum(agree_validation) / nrow(validation_set)

SVM_Accuracy = rbind(SVM_Accuracy, c("Linear", linear_svm_model$bestTune[["C"]], "NA",
                                     "NA", accuracy_train, accuracy_validation))

rm(train_pred, agree_train, accuracy_train, validation_pred, agree_validation,
   accuracy_validation)

# Linear kernel - L2 regularised
RNGversion("4.0.2")
set.seed(123)
L2_linear_svm_model = train(Stage ~., data = train_set, 
                            method = "svmLinear3", 
                            preProcess = NULL,
                            trControl = L2_linear_svm_ctrl,
                            tuneGrid = L2_linear_svm_tune,
                            na.action = "na.omit")

train_pred = predict(L2_linear_svm_model, train_set)
agree_train = ifelse(train_pred == train_set$Stage, 1, 0)
accuracy_train = sum(agree_train) / nrow(train_set)

validation_pred = predict(L2_linear_svm_model, validation_set)
agree_validation = ifelse(validation_pred == validation_set$Stage, 1, 0)
accuracy_validation = sum(agree_validation) / nrow(validation_set)

SVM_Accuracy = rbind(SVM_Accuracy, c("L2_Linear", L2_linear_svm_model$bestTune[["cost"]], 
                                     "NA", "NA", accuracy_train, accuracy_validation))

rm(train_pred, agree_train, accuracy_train, validation_pred, agree_validation,
   accuracy_validation)

# Linear kernel - L2 regularised WITH class weights
RNGversion("4.0.2")
set.seed(123)
L2_weights_linear_svm_model = train(Stage ~., data = train_set, 
                                    method = "svmLinearWeights2", 
                                    preProcess = NULL,
                                    trControl = L2_linear_svm_ctrl,
                                    tuneGrid = L2_linear_svm_tune_weights,
                                    na.action = "na.omit")

train_pred = predict(L2_weights_linear_svm_model, train_set)
agree_train = ifelse(train_pred == train_set$Stage, 1, 0)
accuracy_train = sum(agree_train) / nrow(train_set)

validation_pred = predict(L2_weights_linear_svm_model, validation_set)
agree_validation = ifelse(validation_pred == validation_set$Stage, 1, 0)
accuracy_validation = sum(agree_validation) / nrow(validation_set)

SVM_Accuracy = rbind(SVM_Accuracy, c("L2_Linear_weights", 
                                     L2_weights_linear_svm_model$bestTune[["cost"]], 
                                     "NA",  L2_weights_linear_svm_model$bestTune[["weight"]],
                                     accuracy_train, accuracy_validation))

rm(train_pred, agree_train, accuracy_train, validation_pred, agree_validation,
   accuracy_validation)

# Radial Basis Function
RNGversion("4.0.2")
set.seed(123)
rbf_svm_model = train(Stage ~., data = train_set, 
                      method = "svmRadial", 
                      preProcess = NULL,
                      trControl = rbf_svm_ctrl,
                      tuneGrid = rbf_svm_tune,
                      na.action = "na.omit")

train_pred = predict(rbf_svm_model, train_set)
agree_train = ifelse(train_pred == train_set$Stage, 1, 0)
accuracy_train = sum(agree_train) / nrow(train_set)

validation_pred = predict(rbf_svm_model, validation_set)
agree_validation = ifelse(validation_pred == validation_set$Stage, 1, 0)
accuracy_validation = sum(agree_validation) / nrow(validation_set)

SVM_Accuracy = rbind(SVM_Accuracy, c("RBF", rbf_svm_model$bestTune[["C"]], 
                                     rbf_svm_model$bestTune[["sigma"]], "NA",
                                     accuracy_train, accuracy_validation))

rm(train_pred, agree_train, accuracy_train, validation_pred, agree_validation,
   accuracy_validation)

# Radial Basis Function WITH class weights
RNGversion("4.0.2")
set.seed(123)
rbf_weights_svm_model = train(Stage ~., data = train_set, 
                              method = "svmRadialWeights", 
                              preProcess = NULL,
                              trControl = rbf_svm_ctrl,
                              tuneGrid = rbf_svm_tune_weights,
                              na.action = "na.omit")

train_pred = predict(rbf_weights_svm_model, train_set)
agree_train = ifelse(train_pred == train_set$Stage, 1, 0)
accuracy_train = sum(agree_train) / nrow(train_set)

validation_pred = predict(rbf_weights_svm_model, validation_set)
agree_validation = ifelse(validation_pred == validation_set$Stage, 1, 0)
accuracy_validation = sum(agree_validation) / nrow(validation_set)

SVM_Accuracy = rbind(SVM_Accuracy, c("RBF_weights", rbf_weights_svm_model$bestTune[["C"]], 
                                     rbf_weights_svm_model$bestTune[["sigma"]],
                                     rbf_weights_svm_model$bestTune[["Weight"]],
                                     accuracy_train, accuracy_validation))

rm(train_pred, agree_train, accuracy_train, validation_pred, agree_validation,
   accuracy_validation)

# Saving
SVM_Accuracy = SVM_Accuracy[2:nrow(SVM_Accuracy),]
SVM_models = list(linear_svm_model, L2_linear_svm_model, L2_weights_linear_svm_model,
                  rbf_svm_model, rbf_weights_svm_model)
names(SVM_models) = c("Linear", "L2 Linear", "L2 Linear - weights", "RBF", "RBF - weights")
addWorksheet(svmwb, "SVM_performance")
writeData(svmwb, "SVM_performance", SVM_Accuracy)
saveWorkbook(svmwb, file = "ML/Tumor_stage/SVM.xlsx",
             overwrite = TRUE)
rm(svmwb, linear_svm_model, L2_linear_svm_model, L2_weights_linear_svm_model,
           rbf_svm_model, rbf_weights_svm_model); gc()
