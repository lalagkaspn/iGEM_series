## This script is used to download and pre-process series (GSE) from GEO for PDAC patients (normal and tumor tissues) with samples taken
## from blood (whole blood samples or circulating tumor samples)

#########################################################################################################################################
#########################################################################################################################################

## LIBRARIES
library(dplyr)
library(GEOquery)

#########################################################################################################################################
#########################################################################################################################################

## DOWNLOADING DATA SERIES FROM GEO DATABASE
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

#########################################################################################################################################
#########################################################################################################################################

## NA VALUES IN EACH EXPRESSION MATRIX

# Isolate expression data for each GSE
esets = list()
for (i in 1:length(GEOsets)) {
  esets[[i]] = exprs(GEOsets[[i]])
}; rm(i)
names(esets) = datasets; rm(i)

# Calculate NAs values
na_esets = c()
for(i in 1:length(esets)) {
  na_esets[i] = sum(is.na(esets[[i]]))
}; rm(i)
names(na_esets) = datasets
na_esets

#########################################################################################################################################
#########################################################################################################################################

## PHENOTYPIC DATA

pdata = list()
for(i in 1:length(GEOsets)) {
  pdata[[i]] = pData(GEOsets[[i]])
}
names(pdata) = datasets; rm(i)

## Isolate only needed pdata
filt_pdata = list()

###############
##           ##
## GSE125158 ##
##           ##
###############

## GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125158
## Paper: https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC6447845/
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

##############
##          ##
## GSE49641 ##
##          ##
##############

## GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49641
## Paper: https://pubmed.ncbi.nlm.nih.gov/25069573/
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

##############
##          ##
## GSE74629 ##
##          ##
##############

## GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74629
## Paper: https://pubmed.ncbi.nlm.nih.gov/29617451/
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

##############
##          ##
## GSE18670 ##
##          ##
##############

## GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18670
## Paper: https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC3599097/
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

## ESETS ##
# From GSE18670 keep only samples for tumor (CTC) and non_tumor (P)
esets[["GSE18670"]] = esets$GSE18670[, filt_pdata$GSE18670$GEO_accession]

