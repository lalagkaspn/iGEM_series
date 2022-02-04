## This script is used to download and pre-process series (GSE) from GEO for PDAC
## patients (normal and tumor tissues) who did not receive neoadjuvant chemotherapy. 
## Studies with samples annotated with tumor stage were kept only.

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
library(impute)

##### Downloading data #####
datasets = c("GSE21501", "GSE42952", "GSE18670", "GSE62452", "GSE62165", 
             "GSE102238", "GSE84219")

# Run this before getGEO
# There seems to be a bug with the newest version of readr and GEOquery
# https://github.com/seandavi/GEOquery/issues/114
readr::local_edition(1)

# Download data
GEOsets = list()
for (i in 1:length(datasets)){
  GEOsets[[i]] = getGEO(datasets[i])
}
GEOsets = unlist(GEOsets)
names(GEOsets) = datasets; rm(i)

###### pData #####
pdata = list()
for(i in 1:length(GEOsets)) {
  pdata[[i]] = pData(GEOsets[[i]])
}
names(pdata) = datasets; rm(i)

# TNM classification for PDAC
# https://www.cancerresearchuk.org/about-cancer/pancreatic-cancer/stages-types-grades/tnm-staging
# T: Describes the size of the cancer
#   Tis: carcinoma in situ --> very early stage pancreatic cancer
#   T1: the cancer is inside the pancreas and is 2cm or less in any direction
#    - Τ1a: <0.5cm
#    - Τ1b: more than 0.5cm but less than 1cm
#    - Τ1c: between 1cm and 2cm
#   T2: the cancer is between 2cm and 4cm in any direction
#   T3: the cancer is bigger than 4cm but is still within the pancreas
#   T4: the cancer has grown outside the pancreas, into the nearby large blood vessels
# N: Describes whether there are cancer cells in the lymph nodes
#   N0: there are no cancer cells in the nearby lymph nodes
#   N1: there are 1 to 3 lymph nodes that contain cancer cells
#   N2: there is cancer in 4 or more lymph nodes
# M: Describes whether the cancer has spread to a different part of the body
#   M0: cancer has not spread to other areas of the body
#   M1: cancer has spread to other areas of the body

# Association of AJCC classification (8th edition) and TNM classification for pancreatic cancer
# https://www.nature.com/articles/s41598-018-28193-4/tables/1

# Tumor grade for PDAC (this classification was used in all GSEs except for GSE62452)
# https://www.cancer.org/cancer/pancreatic-cancer/detection-diagnosis-staging/staging.html
#   Grade 1: cancer looks much like normal pancreas tissue
#     - Tend to grow and spread more slowly than high-grade (Grade 3) cancers
#   Grade 2: cancer falls somewhere in between Grade 1 and 3
#   Grade 3: cancer looks very abnormal
#     - Most of the times, tend to have poor diagnosis compared to Grade 1 or 2

# General tumor grades classification (this classification was used in GSE62452)
# https://www.cancer.gov/about-cancer/diagnosis-staging/prognosis/tumor-grade-fact-sheet
#   Grade X: Grade cannot be assessed (undetermined grade)
#   Grade 1: Well differentiated (low grade)
#   Grade 2: Moderately differentiated (intermediate grade)
#   Grade 3: Poorly differentiated (high grade)
#   Grade 4: Undifferentiated (high grade)

# Filter pdata for keeping only necessary information
filt_pdata = list()

# GSE21501
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21501
# Paper: https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC2903589/
# Comments:
#   - Clinical data from "University of Nebraska Medical Center Rapid Autopsy 
#     Pancreatic Program" (NEB) and "Yeh Lab-UNC-CH" (UNC) were 
#     not used in the original analysis. Therefore, authors did not include them
#     in the GSE dataset.
#   - GSM536033: it was used as reference and will be excluded from further analysis.
#   - 97 samples were primary PDAC and 15 were metastatic PDAC. For this reason,
#     we added an extra column for M_stage. 
#     However, for metastatic PDAC patients, they took samples from primary tumor.
#   - From "UNC" samples, 1 patient received neodjuvant chemotherapy and 4 
#     received adjuvant chemotherapy. NO GSM-specific information.
#     All UNC samples will be removed!
#   - From "NEB", 11 patients received chemotherapy <6 months prior to death. 
#     NO GSM-specific information.
#   - From "Northwestern Memorial Hospital" (NW) and "NorthShore University 
#     HealthSystem" (NSU) samples, 2 received neoadjuvant chemotherapy 
#     and 30 received adjuvant chemotherapy. NO GSM-specific information.
#     All NW and NSU patient samples will be removed!

# Select pdata
filt_pdata[["GSE21501"]] = pdata$GSE21501 %>%
  dplyr::select(GEO_accession = geo_accession, 
         Patient_ID = title,
         Platform = platform_id, 
         Biomaterial_provider = biomaterial_provider_ch2,
         Tissue_type = "tissue:ch2",
         T_stage = "t stage:ch2",
         N_stage = "n stage:ch2", 
         Survival_status = "os event:ch2", 
         Survival_months = "os time:ch2")

# Remove GSM536033 from the dataset (position of element: 6)
filt_pdata$GSE21501 = filt_pdata$GSE21501[-6,]

# Rename Biomaterial_provider to agree with the names given in the paper
filt_pdata$GSE21501$Biomaterial_provider = 
  ifelse(filt_pdata$GSE21501$Biomaterial_provider == 
           "University of Nebraska Medical Center Rapid Autopsy Pancreatic Program", "NEB", 
         ifelse(filt_pdata$GSE21501$Biomaterial_provider == 
                  "Yeh Lab-UNC-CH", "UNC",
                ifelse(filt_pdata$GSE21501$Biomaterial_provider == 
                         "Northwestern Memorial Hospital", "NW",
                       ifelse(filt_pdata$GSE21501$Biomaterial_provider == 
                                "NorthShore University HealthSystem" , "NSU", "JHMI"))))

# Annotate metastatic samples (M_stage = 1)
filt_pdata$GSE21501$M_stage = ifelse(filt_pdata$GSE21501$Biomaterial_provider ==
                                       "NEB", "1", "0") 

# Keep only JHMI samples (no neoadjuvant chemotherapy and clinical data availability)
filt_pdata$GSE21501 = filt_pdata$GSE21501 %>% 
  dplyr::filter(Biomaterial_provider == "JHMI")

# Add AJCC classification according to TNM classification
filt_pdata$GSE21501 = filt_pdata$GSE21501 %>%
  mutate(AJCC_classification = c("2b","2b","2b","2b","2b","2b","2b","2b","2b","2b",
                                 "2b","2b","2b","2b","2b","3","2b","2b","1b","2a",
                                 "2b","2b","2b","2b","2b","2b","2b","2b","2b","2b",
                                 "2b","2b","2b","2b"))

# Reorder columns
filt_pdata$GSE21501 = filt_pdata$GSE21501 %>%
  dplyr::select(GEO_accession, Patient_ID, Platform, Biomaterial_provider, 
                Tissue_type, AJCC_classification, T_stage, N_stage, M_stage,
                Survival_status, Survival_months)

# Transform to factors with consistent universal levels
for (i in 1:length(filt_pdata$GSE21501$Tissue_type)) {
  if(filt_pdata$GSE21501$Tissue_type[i] == "Patient Primary Pancreatic Tumor") {
    filt_pdata$GSE21501$Tissue_type[i] = "tumor"
  } else {
    filt_pdata$GSE21501$Tissue_type[i] = "non_tumor"
  }
}; rm(i)

filt_pdata$GSE21501$Tissue_type = factor(x = filt_pdata$GSE21501$Tissue_type, 
                                         levels = c("non_tumor", "tumor"), 
                                         labels = c("non_tumor", "tumor"))
filt_pdata$GSE21501$AJCC_classification = factor(x = filt_pdata$GSE21501$AJCC_classification,
                                                 levels = c("1a","1b","2a","2b","3","4"),
                                                 labels = c("1a","1b","2a","2b","3","4"))
filt_pdata$GSE21501$T_stage = factor(x = filt_pdata$GSE21501$T_stage,
                                     levels = c("1","2","3","4"), 
                                     labels = c("1","2","3","4"))
filt_pdata$GSE21501$N_stage = factor(x = filt_pdata$GSE21501$N_stage,
                                     levels = c("0","1","2"), 
                                     labels = c("0","1","2"))
filt_pdata$GSE21501$M_stage = factor(x = filt_pdata$GSE21501$M_stage, 
                                     levels = c("0", "1"), 
                                     labels = c("non_metastatic", "metastatic"))
filt_pdata$GSE21501$Survival_status = factor(x = filt_pdata$GSE21501$Survival_status,
                                             levels = c("1", "0"), 
                                             labels = c("alive", "dead"))
filt_pdata$GSE21501$Survival_months = as.numeric(filt_pdata$GSE21501$Survival_months)

# GSE42952
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42952
# Paper: https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC3511800/
# Comments:
#   - 34 samples
#       - 17 PDAC (classified as "good" and "bad" outcome)
#       - 6 surrounding non-tumoral pancreatic (control)
#           - 4 patients with tumor and non-tumor samples
#           - 2 patients with non-tumor samples
#       - 11 from metastatic PDAC (liver and peritoneal tissues)
#   - Metastatic samples (LM and PM) were contaminated with normal
#     liver and peritoneal tissue respectively. Due to this fact, it was argued that
#     significant genes identified by DGEA could be attributed to the different
#     tissues and therefore, only genes that were 
#     not differentially expressed between LM and PM samples, 
#     should be considered to be metastasis-specific genes, and these were subsequently used
#     for analysis between primary tumor and metastatic tissue in the original paper.

# Select pData
filt_pdata[["GSE42952"]] = pdata$GSE42952 %>%
  dplyr::select(GEO_accession = geo_accession,
         Patient_ID = title,
         Platform = platform_id,
         Age = "age:ch1",
         Gender = "gender:ch1",
         Prognosis = "prognosis:ch1",
         AJCC_classification = "Stage:ch1",
         Tissue_type = "tissue:ch1",
         TNM_classification = "tumor stage:ch1")

# Samples: "PDAC_P_117", "PDAC_P_104", "PDAC_T_105", "PDAC_P_105", "PDAC_T_91", 
# "PDAC_P_91", "PDAC_T_138", "PDAC_P_138", "PDAC_T_123", 
# are identical to samples from GSE18670. Thus, they will be removed from here 
# (position of elements: 1-10).
filt_pdata$GSE42952 = filt_pdata$GSE42952[-c(1:10),]

# Remove samples from non PDAC tissues (e.g. "liver metastasis from pancreatic cancer",
# "peritoneal metastasis from pancreatic cancer")
filt_pdata$GSE42952 = filt_pdata$GSE42952 %>% 
  dplyr::filter(`Tissue_type` != "liver metastasis from pancreatic cancer" & 
                  `Tissue_type` != "peritoneal metastasis from pancreatic cancer")

# Manual addition of extra data from Table 1 
# (https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC3511800/)
filt_pdata$GSE42952 = filt_pdata$GSE42952 %>% 
  mutate(Tumor_localization = c("head", "head", "head", "head", "head", "head",
                                "head", "head", "head", "head", "head", "head"),
         Tumor_grade = c(NA,1,3,3,2,3,3,3,2,2,3,3),
         T_stage = c(NA,3,3,3,3,3,3,3,3,2,3,3),
         N_stage = c(NA,0,0,1,1,0,0,1,1,0,1,1),
         M_stage = c(1,0,0,0,0,0,0,0,0,0,0,0),
         Perineural_invasion = c(NA,0,1,1,1,1,1,0,1,1,1,1),
         Vascular_invasion = c(NA,0,0,1,1,1,1,1,1,1,0,0)) %>%
  dplyr::select(-TNM_classification)

# Reorder columns
filt_pdata$GSE42952 = filt_pdata$GSE42952 %>%
  dplyr::select(GEO_accession, Patient_ID, Platform, Gender, Age, Tissue_type, 
         AJCC_classification, T_stage, N_stage, M_stage, Tumor_grade, 
         Perineural_invasion, Vascular_invasion, Prognosis, Tumor_localization)

# Transform to factors with consistent universal levels
for (i in 1:length(filt_pdata$GSE42952$Tissue_type)) {
  if(filt_pdata$GSE42952$Tissue_type[i] == "pancreatic ductal adenocarcinoma (PDAC)") {
    filt_pdata$GSE42952$Tissue_type[i] = "tumor"
  } else {
    filt_pdata$GSE42952$Tissue_type[i] = "non_tumor"
  }
}; rm(i)

filt_pdata$GSE42952$Tissue_type = factor(x = filt_pdata$GSE42952$Tissue_type,
                                         levels = c("non_tumor", "tumor"),
                                         labels = c("non_tumor", "tumor"))
filt_pdata$GSE42952$AJCC_classification = factor(x = filt_pdata$GSE42952$AJCC_classification,
                                                 levels = c("1a","1b","2a","2b","3","4"),
                                                 labels = c("1a","1b","2a","2b","3","4"))
filt_pdata$GSE42952$T_stage = factor(x = filt_pdata$GSE42952$T_stage,
                                     levels = c("1","2","3","4"),
                                     labels = c("1","2","3","4"))
filt_pdata$GSE42952$N_stage = factor(x = filt_pdata$GSE42952$N_stage,
                                     levels = c("0","1","2"),
                                     labels = c("0","1","2"))
filt_pdata$GSE42952$M_stage = factor(x = filt_pdata$GSE42952$M_stage,
                                     levels = c("0", "1"),
                                     labels = c("non_metastatic", "metastatic"))
filt_pdata$GSE42952$Tumor_grade = factor(x = filt_pdata$GSE42952$Tumor_grade,
                                         levels = c("1","2","3"),
                                         labels = c("1","2","3"))
filt_pdata$GSE42952$Perineural_invasion = factor(x = filt_pdata$GSE42952$Perineural_invasion,
                                                 levels = c(0, 1),
                                                 labels = c("no_PNI", "PNI"))
filt_pdata$GSE42952$Vascular_invasion = factor(x = filt_pdata$GSE42952$Vascular_invasion,
                                               levels = c(0,1),
                                               labels = c("no_VI", "VI"))
filt_pdata$GSE42952$Gender = factor(x = filt_pdata$GSE42952$Gender,
                                    levels = c("female", "male"),
                                    labels = c("female", "male"))
filt_pdata$GSE42952$Age = gsub("years", "", filt_pdata$GSE42952$Age)
filt_pdata$GSE42952$Age = as.numeric(filt_pdata$GSE42952$Age)
filt_pdata$GSE42952$Prognosis = factor(x = filt_pdata$GSE42952$Prognosis,
                                       levels = c("good", "bad"),
                                       labels = c("good", "bad"))
filt_pdata$GSE42952$Tumor_localization = factor(x = filt_pdata$GSE42952$Tumor_localization,
                                                levels = c("head", "body/tail"),
                                                labels = c("head", "body/tail"))

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

# Keep only tumor (T) and control (P) samples
patients_keep = grep("_T_",x = filt_pdata$GSE18670$Tissue_type)
patients_keep = c(patients_keep, grep("_P_",x = filt_pdata$GSE18670$Tissue_type))

filt_pdata$GSE18670 = filt_pdata$GSE18670[patients_keep, ] ; rm(patients_keep)

# Add M_stage column
filt_pdata$GSE18670$M_stage = "0"

# Replace Tissue_type with tumor and non_tumor values
tumor = which(grepl("_T_", filt_pdata$GSE18670$Tissue_type))
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

# GSE62452
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62452
# Paper: https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC4930741/
# Comments:
#   - 130 samples
#     - 69 pancreatic tumor
#     - 61 adjacent pancreatic non-tumor
#   - For tumor grading, the general classification was used (GX,G1,G2,G3,G4). 
#     In samples from other studies, the PDAC classification was used (G1,G2,G3,G4)

# Select pData
filt_pdata[["GSE62452"]] = pdata$GSE62452 %>%
  dplyr::select(GEO_accession = geo_accession,
         Patient_ID = title,
         Platform = platform_id,
         Tissue_type = "tissue:ch1",
         AJCC_classification = "Stage:ch1",
         Tumor_grade = "grading:ch1")

# We didn't keep information regarding survival status and survival time because 
# classification was not clear.

# Remove letter "G" from Tumor_grade
filt_pdata$GSE62452$Tumor_grade = substring(filt_pdata$GSE62452$Tumor_grade, 2)

# Replace Tissue_type with tumor and non_tumor values
for (i in 1:length(filt_pdata[["GSE62452"]][["Tissue_type"]])) {
  if(filt_pdata[["GSE62452"]][["Tissue_type"]][i] == "Pancreatic tumor") {
    filt_pdata[["GSE62452"]][["Tissue_type"]][i] = "tumor"
  } else {
    filt_pdata[["GSE62452"]][["Tissue_type"]][i] = "non_tumor"
  }
} ; rm(i)

# Transform to factors with consistent universal levels
filt_pdata$GSE62452$Tissue_type = factor(x = filt_pdata$GSE62452$Tissue_type,
                                         levels = c("non_tumor","tumor"),
                                         labels = c("non_tumor","tumor"))
filt_pdata$GSE62452$AJCC_classification = gsub(">IIB", "3", fixed = TRUE,
                                               filt_pdata$GSE62452$AJCC_classification)
filt_pdata$GSE62452$AJCC_classification = gsub("IIB", "2b", fixed = TRUE,
                                               filt_pdata$GSE62452$AJCC_classification)
filt_pdata$GSE62452$AJCC_classification = gsub("IB", "1b", fixed = TRUE,
                                               filt_pdata$GSE62452$AJCC_classification)
filt_pdata$GSE62452$AJCC_classification = gsub("IIA", "2a", fixed = TRUE,
                                               filt_pdata$GSE62452$AJCC_classification)
filt_pdata$GSE62452$AJCC_classification = gsub("III", "3", fixed = TRUE,
                                               filt_pdata$GSE62452$AJCC_classification)
filt_pdata$GSE62452$AJCC_classification = gsub("IVA", "4", fixed = TRUE,
                                               filt_pdata$GSE62452$AJCC_classification)
filt_pdata$GSE62452$AJCC_classification = gsub("IVB", "4", fixed = TRUE,
                                               filt_pdata$GSE62452$AJCC_classification)
filt_pdata$GSE62452$AJCC_classification = gsub("IV", "4", fixed = TRUE,
                                               filt_pdata$GSE62452$AJCC_classification)
filt_pdata$GSE62452$AJCC_classification = factor(x = filt_pdata$GSE62452$AJCC_classification,
                                                 levels = c("1a","1b","2a","2b","3","4"),
                                                 labels = c("1a","1b","2a","2b","3","4"))
filt_pdata$GSE62452$Tumor_grade = factor(x = filt_pdata$GSE62452$Tumor_grade,
                                         levels = c("x","1","2","3", "4"),
                                         labels = c("GX","G1","G2","G3","G4"))

# Add NAs to non_tumor samples' clinical data
for(i in 1:length(filt_pdata$GSE62452$Tissue_type)){
  if(filt_pdata$GSE62452$Tissue_type[i] == "non_tumor") {
    filt_pdata$GSE62452$AJCC_classification[i] = NA
    filt_pdata$GSE62452$Tumor_grade[i] = NA
  }
}; rm(i)

# GSE62165
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62165
# Paper: https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC4983037/
# Comments:
#   - 131 samples
#     - 118 PDAC samples
#     - 13 histologically normal pancreatic tissue samples

# Select pData
filt_pdata[["GSE62165"]] = pdata$GSE62165 %>%
  dplyr::select(GEO_accession = geo_accession,
         Patient_ID = title,
         Platform = platform_id,
         Tissue_type =  "tissue:ch1",
         AJCC_classification = "Stage:ch1")

# Replace Tissue_type with tumor and non_tumor values
for (i in 1:length(filt_pdata[["GSE62165"]][["Tissue_type"]])) {
  if(filt_pdata[["GSE62165"]][["Tissue_type"]][i] == "pancreatic tumor") {
    filt_pdata[["GSE62165"]][["Tissue_type"]][i] = "tumor"
  } else {
    filt_pdata[["GSE62165"]][["Tissue_type"]][i] = "non_tumor"
  }
}; rm(i)

# Transform to factors with consistent universal levels
filt_pdata$GSE62165$Tissue_type = factor(x = filt_pdata$GSE62165$Tissue_type,
                                         levels = c("non_tumor","tumor"),
                                         labels = c("non_tumor","tumor"))
filt_pdata$GSE62165$AJCC_classification = factor(x = filt_pdata$GSE62165$AJCC_classification,
                                                 levels = c("1a","1b","2a","2b","3","4"),
                                                 labels = c("1a","1b","2a","2b","3","4"))

# Add NAs to non_tumor samples clinical data
for(i in 1:length(filt_pdata$GSE62165$Tissue_type)){
  if(filt_pdata$GSE62165$Tissue_type[i] == "non_tumor") {
    filt_pdata$GSE62165$AJCC_classification[i] = NA
  }
}; rm(i)

# GSE102238
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102238
# Paper: https://cancerres.aacrjournals.org/content/80/10/1991.long
# Comments:
#   - 100 samples (50 patients)
#     - 50 tumor samples with/without perineural invasion
#     - 50 non-tumor samples with/without perineural invasion

# Select pData
filt_pdata[["GSE102238"]] = pdata$GSE102238 %>%
  dplyr::select(GEO_accession = geo_accession,
         Patient_ID = title,
         Platform = platform_id,
         Gender = "genderï¼ˆmale=0ï¼Œfemale=1ï¼‰:ch1",     # column 45               # 0: male           1: female
         Age_grouped = "ageï¼ˆâ‰¤65=0ï¼Œï¼ž65=1ï¼‰:ch1",   # column 43               # 0: <=65           1: >65
         Tissue_type =  "source_name_ch1",
         T_stage = "t stage:ch1",
         N_stage = "n stageï¼ˆâ€\u009dn0=0,n1=1):ch1",     # column 48 
         M_stage = "m stage(m0=0.m1=1):ch1",
         Vascular_invasion = "vessel invasion(absence=0,present=1):ch1",             # 0: absence        1: present
         Differentiation = "differentiation(well/moderate=0,poor=1):ch1",            # 0: well/moderate  1: poor
         Survival_status = "survival status(alive=0,dead=1):ch1",                    # 0: alive          1: dead
         Survival_days = "survival time(day):ch1",
         Tumor_localization = "localization of tumor(head=1,body/tail=2):ch1")       # 1: head           2: body/tail

# Separate information in Tissue_type
# Patient_ID
library(readr)
filt_pdata$GSE102238$Patient_ID = parse_number(filt_pdata[["GSE102238"]][["Tissue_type"]])

# PNI vs no_PNI
for (i in 1:length(filt_pdata$GSE102238$Tissue_type)) {
  if(grepl("without PNI", x = filt_pdata$GSE102238$Tissue_type[[i]]) == TRUE) {
    filt_pdata$GSE102238[i, "Perineural_invasion"] = "no_PNI"
  }
  if(grepl("with PNI", x = filt_pdata$GSE102238$Tissue_type[[i]] ) == TRUE) {
    filt_pdata$GSE102238[i, "Perineural_invasion"] = "PNI"    
  }
} ; rm(i)

# tumor vs non_tumor
for (i in 1:length(filt_pdata$GSE102238$Tissue_type)) {
  if(grepl("normal tissue", x = filt_pdata$GSE102238$Tissue_type[[i]]) == TRUE) {
    filt_pdata$GSE102238[i, "Tissue_type"] = "non_tumor"
  }
  if(grepl("tumor tissue", x = filt_pdata$GSE102238$Tissue_type[[i]] ) == TRUE) {
    filt_pdata$GSE102238[i, "Tissue_type"] = "tumor"    
  }
} ; rm(i)

# Transform to factors with consistent universal levels
filt_pdata$GSE102238$Gender = factor(x = filt_pdata$GSE102238$Gender, 
                                     levels = c("1","0"), 
                                     labels = c("female","male"))
filt_pdata$GSE102238$Age_grouped = factor(x = filt_pdata$GSE102238$Age_grouped,
                                          levels = c("0","1"), 
                                          labels = c("<=65", ">65"))
filt_pdata$GSE102238$Tissue_type = factor(x = filt_pdata$GSE102238$Tissue_type,
                                          levels = c("non_tumor","tumor"),
                                          labels = c("non_tumor","tumor"))
filt_pdata$GSE102238$T_stage = factor(x = filt_pdata$GSE102238$T_stage,
                                      levels = c("1","2","3","4"),
                                      labels = c("1","2","3","4"))
filt_pdata$GSE102238$N_stage = factor(x = filt_pdata$GSE102238$N_stage,
                                      levels = c("0","1","2"),
                                      labels = c("0","1","2"))
filt_pdata$GSE102238$M_stage = factor(x = filt_pdata$GSE102238$M_stage,
                                      levels = c("0", "1"),
                                      labels = c("non_metastatic", "metastatic"))
filt_pdata$GSE102238$Perineural_invasion = factor(x = filt_pdata$GSE102238$Perineural_invasion,
                                                  levels = c("no_PNI","PNI"),
                                                  labels = c("no_PNI","PNI"))
filt_pdata$GSE102238$Vascular_invasion = factor(x = filt_pdata$GSE102238$Vascular_invasion,
                                                levels = c("0","1"),
                                                labels = c("no_VI","VI"))
filt_pdata$GSE102238$Differentiation = factor(x = filt_pdata$GSE102238$Differentiation,
                                              levels = c("0","1"),
                                              labels = c("well/moderate","poor"))
filt_pdata$GSE102238$Survival_status = factor(x = filt_pdata$GSE102238$Survival_status,
                                              levels = c("0","1"),
                                              labels = c("alive","dead"))
filt_pdata$GSE102238$Tumor_localization = factor(x = filt_pdata$GSE102238$Tumor_localization,
                                                 levels = c("1","2"),
                                                 labels = c("head", "body/tail"))
# Transform survival_days to survival_months
filt_pdata$GSE102238$Survival_days = as.numeric(filt_pdata$GSE102238$Survival_days)
filt_pdata$GSE102238 = filt_pdata$GSE102238 %>% mutate(Survival_months = Survival_days/30)

# Add AJCC classification according to TNM classification
for (i in 1:length(filt_pdata$GSE102238$T_stage)) {
  if(filt_pdata$GSE102238$T_stage[i] == 1 &&
     filt_pdata$GSE102238$N_stage[i] == 0 &&
     filt_pdata$GSE102238$M_stage[i] == "non_metastatic") {
    filt_pdata$GSE102238$AJCC_classification[i] = "1a"
  }
  if(filt_pdata$GSE102238$T_stage[i] == 2 &&
     filt_pdata$GSE102238$N_stage[i] == 0 &&
     filt_pdata$GSE102238$M_stage[i] == "non_metastatic") {
    filt_pdata$GSE102238$AJCC_classification[i] = "1b"
  }
  if(filt_pdata$GSE102238$T_stage[i] == 3 &&
     filt_pdata$GSE102238$N_stage[i] == 0 &&
     filt_pdata$GSE102238$M_stage[i] == "non_metastatic") {
    filt_pdata$GSE102238$AJCC_classification[i] = "2a"
  }
  if(filt_pdata$GSE102238$T_stage[i] %in% 1:3 &&
     filt_pdata$GSE102238$N_stage[i] == 1 &&
     filt_pdata$GSE102238$M_stage[i] == "non_metastatic") {
    filt_pdata$GSE102238$AJCC_classification[i] = "2b"
  }
  if(filt_pdata$GSE102238$T_stage[i] == 4 &&
     filt_pdata$GSE102238$N_stage[i] %in% 0:2 &&
     filt_pdata$GSE102238$M_stage[i] == "non_metastatic") {
    filt_pdata$GSE102238$AJCC_classification[i] = "3"
  }
  if(filt_pdata$GSE102238$T_stage[i] %in% 1:4 &&
     filt_pdata$GSE102238$N_stage[i] == 2 &&
     filt_pdata$GSE102238$M_stage[i] == "non_metastatic") {
    filt_pdata$GSE102238$AJCC_classification[i] = "3"
  }
  if(filt_pdata$GSE102238$T_stage[i] %in% 1:4 &&
     filt_pdata$GSE102238$N_stage[i] %in% 0:2 &&
     filt_pdata$GSE102238$M_stage[i] == "metastatic") {
    filt_pdata$GSE102238$AJCC_classification[i] = "4"
  }
}; rm(i)

filt_pdata$GSE102238$AJCC_classification = factor(x = filt_pdata$GSE102238$AJCC_classification,
                                                  levels = c("1a","1b","2a","2b","3","4"),
                                                  labels = c("1a","1b","2a","2b","3","4"))

# Reorder columns
filt_pdata$GSE102238 = filt_pdata$GSE102238 %>%
  dplyr::select(GEO_accession:Tissue_type, AJCC_classification, T_stage:M_stage,
         Perineural_invasion, Vascular_invasion, Differentiation,
         Survival_status, Survival_days, Survival_months,Tumor_localization)

# Add NAs to non_tumor samples clinical data
for(i in 1:length(filt_pdata$GSE102238$Tissue_type)){
  if(filt_pdata$GSE102238$Tissue_type[i] == "non_tumor") {
    filt_pdata$GSE102238$AJCC_classification[i] = NA
    filt_pdata$GSE102238$T_stage[i] = NA
    filt_pdata$GSE102238$N_stage[i] = NA
    filt_pdata$GSE102238$M_stage[i] = NA
    filt_pdata$GSE102238$Perineural_invasion[i] = NA
    filt_pdata$GSE102238$Vascular_invasion[i] = NA
    filt_pdata$GSE102238$Differentiation[i] = NA
    filt_pdata$GSE102238$Survival_status[i] = NA
    filt_pdata$GSE102238$Survival_days[i] = NA
    filt_pdata$GSE102238$Survival_months[i] = NA
    filt_pdata$GSE102238$Tumor_localization[i] = NA
  }
}; rm(i)

# GSE84219
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84219
# Paper: https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC5739601/
# Comments:
#   - 30 samples (30 patients)
#     - 15 short survival (lower than 24 months)
#     - 15 long survival (higher than 24 months)

# Select pData
filt_pdata[["GSE84219"]] = pdata$GSE84219 %>%
  dplyr::select(GEO_accession = geo_accession,
         Patient_ID = description,
         Platform = platform_id,
         AJCC_classification = "tumor_stage:ch1",
         Survival_months = "months_survival:ch1")

# Add Tisse_type column
filt_pdata$GSE84219$Tissue_type = rep("tumor", length(filt_pdata$GSE84219$GEO_accession))

# Transform to factors with consistent universal levels
filt_pdata$GSE84219$Tissue_type = factor(x = filt_pdata$GSE84219$Tissue_type,
                                         levels = c("non_tumor","tumor"),
                                         labels = c("non_tumor","tumor"))
filt_pdata$GSE84219$AJCC_classification = gsub("IIB", "2b", fixed = TRUE,
                                               filt_pdata$GSE84219$AJCC_classification)
filt_pdata$GSE84219$AJCC_classification = gsub("IB", "1b", fixed = TRUE,
                                               filt_pdata$GSE84219$AJCC_classification)
filt_pdata$GSE84219$AJCC_classification = gsub("IIA", "2a", fixed = TRUE,
                                               filt_pdata$GSE84219$AJCC_classification)
filt_pdata$GSE84219$AJCC_classification = gsub("IA", "1a", fixed = TRUE,
                                               filt_pdata$GSE84219$AJCC_classification)
filt_pdata$GSE84219$AJCC_classification = factor(x = filt_pdata$GSE84219$AJCC_classification,
                                                 levels = c("1a","1b","2a","2b","3","4"),
                                                 labels = c("1a","1b","2a","2b","3","4"))

# Reorder columns
filt_pdata$GSE84219 = filt_pdata$GSE84219 %>%
  dplyr::select(GEO_accession, Patient_ID, Platform, Tissue_type, AJCC_classification,
         Survival_months)

# full_pdata
# Keep only information for Study, GEO_accession, Tissue_type and AJCC classification
# Useful for QC analysis
pdataGSE21501 = filt_pdata$GSE21501 %>%
  dplyr::select(GEO_accession, Tissue_type, AJCC_classification) %>%
  dplyr::mutate(Study = "GSE21501")
pdataGSE42952 = filt_pdata$GSE42952 %>%
  dplyr::select(GEO_accession, Tissue_type, AJCC_classification) %>%
  dplyr::mutate(Study = "GSE42952")
pdataGSE18670 = filt_pdata$GSE18670 %>%
  dplyr::select(GEO_accession, Tissue_type, AJCC_classification) %>%
  dplyr::mutate(Study = "GSE18670")
pdataGSE62452 = filt_pdata$GSE62452 %>%
  dplyr::select(GEO_accession, Tissue_type, AJCC_classification) %>%
  dplyr::mutate(Study = "GSE62452")
pdataGSE62165 = filt_pdata$GSE62165 %>%
  dplyr::select(GEO_accession, Tissue_type, AJCC_classification) %>%
  dplyr::mutate(Study = "GSE62165")
pdataGSE102238 = filt_pdata$GSE102238 %>%
  dplyr::select(GEO_accession, Tissue_type, AJCC_classification) %>%
  dplyr::mutate(Study = "GSE102238")
pdataGSE84219 = filt_pdata$GSE84219 %>%
  dplyr::select(GEO_accession, Tissue_type, AJCC_classification) %>%
  dplyr::mutate(Study = "GSE84219")

full_pdata = rbind(pdataGSE21501, pdataGSE42952, pdataGSE18670, pdataGSE62452,
                   pdataGSE62165, pdataGSE102238, pdataGSE84219)
rownames(full_pdata) = full_pdata$GEO_accession
rm(pdataGSE21501, pdataGSE42952, pdataGSE18670, pdataGSE62452, pdataGSE62165,
   pdataGSE102238, pdataGSE84219)

##### Expression data #####
# GSE42952 has a lot of missing values. Certain columns (~40% of the samples)
# have more than 80% missing values (on the same rows) and therefore imputation
# could yield quite biased results and is avoided in this case. The rows with
# more than 25% missing values are removed.
GEOsets[["GSE42952"]] = GEOsets[["GSE42952"]][rowSums(is.na(GEOsets[["GSE42952"]]@assayData[["exprs"]]))/
                              length(colnames(GEOsets[["GSE42952"]]@assayData[["exprs"]])) < 0.25, ]

# GSE102238 is the only other gene expression matrix with a few missing values (23)
# which due to the large sample (100) will be imputed using kNN imputation:
RNGversion("4.0.2")
eset102238 = GEOsets[["GSE102238"]]@assayData[["exprs"]]
eset102238 = impute.knn(eset102238, k = 10, maxp = nrow(eset102238),
                       rng.seed = 123)
eset102238 = eset102238[["data"]]

# Create the esets objects
esets = list()
for (i in 1:length(GEOsets)) {
  esets[[i]] = exprs(GEOsets[[i]])
  esets[[i]] = as.data.frame(esets[[i]])
}
names(esets) = datasets; rm(i)
esets[["GSE102238"]] = as.data.frame(eset102238); rm(eset102238)

# GSE21501
# Remove sample GSM506033 (6th in order) because it was used as reference 
# (excluded from pdata, too)
esets[["GSE21501"]] = esets[["GSE21501"]][, -6]
# Keep only expression data for JHMI samples (no neoadjuvant chemotherapy)
esets$GSE21501 = esets$GSE21501[, colnames(esets$GSE21501) %in%
                                  filt_pdata[["GSE21501"]][["GEO_accession"]]]

# GSE42952
# Remove first 10 samples because they are identical to the ones of GSE18670
# Remove non-PDAC tissue samples
esets[["GSE42952"]] = esets[["GSE42952"]][, -c(1:10, 13,14,17,18,21,22,25,26,28,30,33)]

# GSE18670 
# Keep only samples for tumor (T) and non_tumor (P)
esets[["GSE18670"]] = esets$GSE18670[, c("GSM463723", "GSM463724", "GSM463727",
                                         "GSM463728", "GSM463731", "GSM463732", 
                                         "GSM463735", "GSM463736", "GSM463739",
                                         "GSM463740", "GSM463743", "GSM463744")]

##### Annotation with Entrez ID's #####
# GSE21501
fdata21501 = fData(GEOsets$GSE21501) %>%
  dplyr::select(ID, ENTREZ_GENE_ID = GENE) %>%
  dplyr::filter(!grepl("///", ENTREZ_GENE_ID)) %>%
  dplyr::filter(nchar(ENTREZ_GENE_ID)>0)
esets[["GSE21501"]]$ID = as.numeric(rownames(esets[["GSE21501"]]))
esets[["GSE21501"]] = esets[["GSE21501"]] %>%
  inner_join(fdata21501) %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  summarise_all(mean, na.rm = TRUE)
rm(fdata21501)

# GSE42952
fdata42952 = fData(GEOsets$GSE42952) %>%
  dplyr::select(ID, ENTREZ_GENE_ID) %>%
  dplyr::filter(!grepl("///", ENTREZ_GENE_ID)) %>%
  dplyr::filter(nchar(ENTREZ_GENE_ID)>0)
esets[["GSE42952"]]$ID = rownames(esets[["GSE42952"]])
esets[["GSE42952"]] = esets[["GSE42952"]] %>%
  inner_join(fdata42952) %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  summarise_all(mean, na.rm = TRUE)
rm(fdata42952)

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

# GSE62452
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

fdata62452 = fData(GEOsets$GSE62452) %>%
  dplyr::select(ID, RefSeq = GB_LIST) %>%
  dplyr::filter(is.na(RefSeq)==F) %>%
  dplyr::filter(nchar(RefSeq)>0)

names = c()
for (i in 1:202){
  names = c(names, paste0("RS", i))
}
fdata62452$ID = as.character(fdata62452$ID)
fdata62452_sep = melt(separate(fdata62452, col = RefSeq, into = names, sep = ","),
                      id.vars = "ID") %>%
  dplyr::select(-variable) %>%
  dplyr::rename(RefSeq = value) %>%
  dplyr::filter(is.na(RefSeq)==F) %>%
  dplyr::filter(nchar(RefSeq)>0) %>%
  inner_join(ref_df, by = "RefSeq") %>%
  distinct() %>%
  group_by(ID)
fdata62452_sep$Filter = NA;

# Keep only the probes which match to a unique Entrez, even if these mapped to
# multiple RefSeq ID's
for (i in 1:nrow(fdata62452_sep)){
  if (length(unique(fdata62452_sep$ENTREZ_GENE_ID[fdata62452_sep$ID==fdata62452_sep$ID[i]]))==1){
    fdata62452_sep$Filter[i] = "keep"
  } else {
    fdata62452_sep$Filter[i] = "discard"
  }
}

final_fdata62452 = fdata62452_sep %>%
  dplyr::filter(Filter == "keep") %>%
  distinct() %>% 
  dplyr::select(-RefSeq, -Filter)
esets[["GSE62452"]]$ID = as.character(rownames(esets[["GSE62452"]]))
esets[["GSE62452"]] = esets[["GSE62452"]] %>%
  inner_join(final_fdata62452, by = "ID") %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  summarise_all(mean, na.rm = TRUE)
rm(fdata62452, fdata62452_sep, final_fdata62452, names)

# GSE62165
fdata62165 = fData(GEOsets$GSE62165) %>%
  dplyr::select(ID, ENTREZ_GENE_ID = 'Entrez Gene') %>%
  dplyr::filter(!grepl("///", ENTREZ_GENE_ID)) %>%
  dplyr::filter(nchar(ENTREZ_GENE_ID)>0) %>%
  dplyr::filter(!grepl("---", ENTREZ_GENE_ID))
esets[["GSE62165"]]$ID = rownames(esets[["GSE62165"]])
esets[["GSE62165"]] = esets[["GSE62165"]] %>%
  inner_join(fdata62165) %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  summarise_all(mean, na.rm = TRUE)
rm(fdata62165)

# GSE102238
# No information for corresponding genes provided in GEOset
# Two separate .fasta files were prepared based on the probe sequences
# available at the GPL file.
# BLASTN alignment (against Nucleotide sequences, "nt") was performed
# to annotate probes with RefSeq gene IDs
# https://blast.ncbi.nlm.nih.gov/Blast.cgi
# Read csv results from BLASTN

blastn_1 = read.csv(file = "Sequences_1.csv", header = FALSE)
blastn_2 = read.csv(file = "Sequences_2.csv", header = FALSE)
blastn_3 = read.csv(file = "Sequences_3.csv", header = FALSE)
blastn_4 = read.csv(file = "Sequences_4.csv", header = FALSE)
blastn = rbind(blastn_1, blastn_2, blastn_3, blastn_4)
blastn$V2 = gsub("\\..", "", blastn$V2) # RefSeq column: keep main nomenclature - ignore variants
blastn = blastn %>%
  dplyr::select(V1, V2, V3, V11) %>%
  dplyr::filter(V3 == 100.00) %>%
  distinct()
rm(blastn_1, blastn_2, blastn_3, blastn_4)
colnames(blastn) = c("probe", "RefSeq", "Alignment_perc", "E_value")
blastn$Alignment_perc = as.numeric(blastn$Alignment_perc)
blastn$E_value = as.numeric(blastn$E_value)
blastn = blastn[order(blastn$probe, 1/blastn$E_value),]

# Map through org.Hs.eg.db
mapped_blastn = inner_join(blastn, ref_df, by = "RefSeq") %>%
  dplyr::select(probe, RefSeq, ENTREZ_GENE_ID) %>%
  distinct() %>%
  dplyr::select(probe, ENTREZ_GENE_ID) # we do not use distinct() here for filtering purposes

# Probes which match to multiple ID's: 
# Check if >50% of Entrez ID's mapping to a probe are actually a unique Entrez ID
# If yes, map the probe to that Entrez ID. If not, discard the probe
# Process described here and suggested by Ensembl: 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3431719/

mapped_blastn = mapped_blastn[order(mapped_blastn$probe),]
mapped_blastn$unique_Entrez_perc = NA

# Add a column to filter as described previously
for (i in 1:nrow(mapped_blastn)){
  if (length(unique(mapped_blastn$ENTREZ_GENE_ID[mapped_blastn$probe == 
                                                 mapped_blastn$probe[i]])) == 1) {
    mapped_blastn$unique_Entrez_perc[i] = 100
  } else {
    mapped_blastn$unique_Entrez_perc[i] = 100*length(which(mapped_blastn$ENTREZ_GENE_ID[mapped_blastn$probe == 
                                                                                          mapped_blastn$probe[i]] == mapped_blastn$ENTREZ_GENE_ID[i]))/
      nrow(mapped_blastn[mapped_blastn$probe == mapped_blastn$probe[i],])
  }
}

# We now keep everything with unique Entrez mapping percentage over 50%
# That is quaranteed to keep one Entrez ID for each probe and discard probes 
# for which only lower mapping percentages exist

mapped_blastn_filt = mapped_blastn %>%
  dplyr::filter(unique_Entrez_perc > 50) %>%
  dplyr::select(probe, ENTREZ_GENE_ID) %>%
  distinct()

# > length(which(duplicated(mapped_blastn_filt$probe)))
# [1] 0
# > length(which(duplicated(mapped_blastn_filt$ENTREZ_GENE_ID)))
# [1] 4818

# That means each probe is mapped to a unique Entrez ID, but multiple probes
# may map to the same ID. These probes will be averaged after the gene expression
# matrix is annotated, as we did in previous cases.

esets[["GSE102238"]]$probe = rownames(esets[["GSE102238"]])
esets[["GSE102238"]] = esets[["GSE102238"]] %>%
  inner_join(mapped_blastn_filt, by = "probe") %>%
  dplyr::select(-probe) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  summarise_all(mean, na.rm = TRUE)
rm(blastn, mapped_blastn, mapped_blastn_filt)


# GSE84219
fdata84219 = fData(GEOsets$GSE84219) %>%
  dplyr::select(ID, ENTREZ_GENE_ID = Entrez_Gene_ID) %>%
  dplyr::filter(!grepl("///", ENTREZ_GENE_ID)) %>%
  dplyr::filter(nchar(ENTREZ_GENE_ID)>0)
esets[["GSE84219"]]$ID = rownames(esets[["GSE84219"]])
esets[["GSE84219"]] = esets[["GSE84219"]] %>%
  inner_join(fdata84219) %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  summarise_all(mean, na.rm = TRUE)
rm(fdata84219)

##### Calculate NAs values #####
na_esets = c()
for(i in 1:length(esets)) {
  na_esets[i] = sum(is.na(esets[[i]]))
}; rm(i)
names(na_esets) = datasets
na_esets

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
  z[[i]]$EntrezGene.ID = esets[[i]]$ENTREZ_GENE_ID
  rm(t, z_t, df)
}; rm(i)

##### Quality Control #####

# We keep the complete-case join for these plots (intersection)

# Depending on whether the poorly annotated experiment GSE102238 will be included
# we define our full pheno data differently. Here it is excluded.

# full_pdata_filt = full_pdata[full_pdata$Study != "GSE102238", ]
# In order to include GSE102238 change the previous line of code to:
full_pdata_filt = full_pdata

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
  inner_join(esets[[5]], by = "ENTREZ_GENE_ID") %>%
  inner_join(esets[[6]], by = "ENTREZ_GENE_ID") %>%
  inner_join(esets[[7]], by = "ENTREZ_GENE_ID") %>%
  dplyr::select(ENTREZ_GENE_ID, everything())

rows = original_exprs$ENTREZ_GENE_ID
original_exprs = as.matrix(original_exprs %>% dplyr::select(-ENTREZ_GENE_ID))
rownames(original_exprs) = rows; rm(rows) # 2601 x 449
# Making sure we do not have NAs in any row
original_exprs = original_exprs[rowSums(is.na(original_exprs)) != ncol(original_exprs), ]
original_exprs_nonas = na.omit(original_exprs) # 2601 x 449

# Joining in one expression matrix: z-score normalised version
for (i in 1:length(z)){
  z[[i]]$EntrezGene.ID = as.character(z[[i]]$EntrezGene.ID)
  z[[i]] = z[[i]][,c("EntrezGene.ID", 
                             intersect(colnames(z[[i]]),
                                       filt_pdata[[i]]$GEO_accession))]
}
z_exprs = z[[1]] %>% inner_join(z[[2]], by = "EntrezGene.ID") %>%
  inner_join(z[[3]], by = "EntrezGene.ID") %>%
  inner_join(z[[4]], by = "EntrezGene.ID") %>%
  inner_join(z[[5]], by = "EntrezGene.ID") %>%
  inner_join(z[[6]], by = "EntrezGene.ID") %>%
  inner_join(z[[7]], by = "EntrezGene.ID") %>%
  dplyr::select(EntrezGene.ID, everything())

rownames(z_exprs) = z_exprs$EntrezGene.ID
z_exprs = as.matrix(z_exprs %>% dplyr::select(-EntrezGene.ID)) # 2601 x 449
# Making sure we do not have NAs in any row
z_exprs = z_exprs[rowSums(is.na(z_exprs)) != ncol(z_exprs), ]
z_exprs_nonas = na.omit(z_exprs) # 2601 x 449

# Multidimensional scaling plot: original matrix #####
original_mds = plotMDS(original_exprs_nonas)
original_pca = data.frame(cbind(original_mds$x, original_mds$y, 
                                as.character(full_pdata_filt$Study), 
                                full_pdata_filt$GEO_accession, 
                                as.character(full_pdata_filt$Tissue_type)))
colnames(original_pca) = c("X1", "X2", "Study", "GSM", "Type")
original_pca$Study = factor(original_pca$Study)
original_pca$Type = factor(original_pca$Type)
original_pca$X1 = as.numeric(original_pca$X1)
original_pca$X2 = as.numeric(original_pca$X2)

png("Plots/QC/Tumor_stage/Original_MDS.png", width = 1024, height = 768)
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
                         as.character(full_pdata_filt$Study), full_pdata_filt$GEO_accession, 
                         as.character(full_pdata_filt$Tissue_type)))
colnames(z_pca) = c("X1", "X2", "Study", "GSM", "Type")
z_pca$Study = factor(z_pca$Study)
z_pca$Type = factor(z_pca$Type)
z_pca$X1 = as.numeric(z_pca$X1)
z_pca$X2 = as.numeric(z_pca$X2)

png("Plots/QC/Tumor_stage/KBZ_MDS.png", width = 1024, height = 768)
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
png("Plots/QC/Tumor_stage/Original_boxplot.png", width = 1920, height = 1080)
ggplot(melt(original_eset), aes(x=variable, y=value)) +
  geom_boxplot(outlier.size = 0.4, outlier.shape = 20,
               fill = c(rep("cyan", 34), rep("chartreuse", 12),
                        rep("orange", 12), rep("red", 130), rep("grey", 131),
                        rep("purple", 100), 
                        rep("pink", 30)), outlier.alpha = 0.1)+
  scale_y_continuous("Expression", limits = c(0,round(max(melt(original_eset)$value)+1)), 
                     breaks = seq(0,round(max(melt(original_eset)$value)+1), 1))+
  theme(plot.title = element_text(face = "bold", size = 27, hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.text.y = element_text(angle = 0, hjust = 1, margin = margin(t = 1, unit = "cm"),
                                   size = 5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5, 
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
png("Plots/QC/Tumor_stage/KBZ_boxplot.png", width = 1920, height = 1080)
ggplot(melt(z_eset), aes(x=variable, y=value)) +
  geom_boxplot(outlier.size = 0.4, outlier.shape = 20,
               fill = c(rep("cyan", 34), rep("chartreuse", 12),
                        rep("orange", 12), rep("red", 130), rep("grey", 131),
                        rep("purple", 100), 
                        rep("pink", 30)), outlier.alpha = 0.1)+
  scale_y_continuous("Expression", limits = c(0,round(max(melt(z_eset)$value)+1)), 
                     breaks = seq(0,round(max(melt(z_eset)$value)+1), 1))+
  theme(plot.title = element_text(face = "bold", size = 27, hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.text.y = element_text(angle = 0, hjust = 1, margin = margin(t = 1, unit = "cm"),
                                   size = 5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5, 
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

# Original
annotation_for_heatmap = full_pdata_filt[, c("Study", "Tissue_type")]
rownames(annotation_for_heatmap) = full_pdata_filt$GEO_accession

original_dists = as.matrix(dist(t(original_exprs_nonas), method = "manhattan"))

rownames(original_dists) = full_pdata_filt$GEO_accession
hmcol = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255))
colnames(original_dists) <- NULL
diag(original_dists) <- NA

ann_colors <- list(
  Tissue_type = c(tumor = "deeppink4", non_tumor = "dodgerblue4"),
  Study = c(GSE21501 = "darkseagreen", GSE42952 = "darkorange",
            GSE18670 = "darkcyan", GSE62452 = "darkred",
            GSE62165 = "grey", GSE102238 = "darkmagenta", 
            GSE84219 = "yellow")
)

original_heatmap = pheatmap(t(original_dists), col = hmcol,
                            annotation_col = annotation_for_heatmap,
                            annotation_colors = ann_colors,
                            legend = TRUE,
                            treeheight_col = 0,
                            legend_breaks = c(min(original_dists, na.rm = TRUE), 
                                              max(original_dists, na.rm = TRUE)), 
                            legend_labels = (c("small distance", "large distance")),
                            main = "Original heatmap")
save_pheatmap_png(original_heatmap, "Plots/QC/Tumor_stage/original_heatmap.png")

# Z-score version
annotation_for_heatmap = full_pdata_filt[, c("Study", "Tissue_type")]
rownames(annotation_for_heatmap) = full_pdata_filt$GEO_accession

z_dists = as.matrix(dist(t(z_exprs_nonas), method = "manhattan"))

rownames(z_dists) = full_pdata_filt$GEO_accession
hmcol = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255))
colnames(z_dists) <- NULL
diag(z_dists) <- NA

ann_colors <- list(
  Tissue_type = c(tumor = "deeppink4", non_tumor = "dodgerblue4"),
  Study = c(GSE21501 = "darkseagreen", GSE42952 = "darkorange",
            GSE18670 = "darkcyan", GSE62452 = "darkred",
            GSE62165 = "grey", GSE102238 = "darkmagenta", 
            GSE84219 = "yellow")
)

z_heatmap = pheatmap(t(z_dists), col = hmcol,
                     annotation_col = annotation_for_heatmap,
                     annotation_colors = ann_colors,
                     legend = TRUE,
                     treeheight_col = 0,
                     legend_breaks = c(min(z_dists, na.rm = TRUE), 
                                       max(z_dists, na.rm = TRUE)), 
                     legend_labels = (c("small distance", "large distance")),
                     main = "Z-score normalisation heatmap")
save_pheatmap_png(z_heatmap, "Plots/QC/Tumor_stage/KBZ_heatmap.png")

##### Differential Gene Expression (DGEA) #####
# 4 types of DGEA:l
#    - Stages 1/2 vs Stages 3/4
#    - Stage 1 vs Stage 4
#    - Stage 1 vs normal
#    - Stage 2 vs normal

full_pdata_filt$Study = as.factor(full_pdata_filt$Study)

# Annotation with official gene symbols
# official_df, Aliases, ID_Map
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

##### Union #####
# In this section we perform DGEA on the union of the gene expression matrices,
# no the intersection (as we did before), in order to keep all genes from all
# platforms. NA's will be introduced in this manner, but limma ignores them
# during model fitting.

# Generating our union z-score matrix
union_z_exprs = z[[1]] %>% full_join(z[[2]], by = "EntrezGene.ID") %>%
  full_join(z[[3]], by = "EntrezGene.ID") %>%
  full_join(z[[4]], by = "EntrezGene.ID") %>%
  full_join(z[[5]], by = "EntrezGene.ID") %>%
  full_join(z[[6]], by = "EntrezGene.ID") %>%
  full_join(z[[7]], by = "EntrezGene.ID") %>%
  dplyr::select(EntrezGene.ID, everything())

rownames(union_z_exprs) = union_z_exprs$EntrezGene.ID
union_z_exprs = as.matrix(union_z_exprs %>% dplyr::select(-EntrezGene.ID)) # 25699 x 449

# Stages 1/2 vs stages 3/4 #####
Stages_pdata = full_pdata_filt %>%
  dplyr::filter(Tissue_type == "tumor") %>%
  dplyr::filter(is.na(AJCC_classification) == FALSE) %>%
  mutate(Stage_group = NA)
Stages_pdata$Stage_group[Stages_pdata$AJCC_classification == 3 | 
                           Stages_pdata$AJCC_classification == 4] = "Advanced_stage"
Stages_pdata$Stage_group[Stages_pdata$AJCC_classification == "1a" |
                           Stages_pdata$AJCC_classification == "1b" |
                           Stages_pdata$AJCC_classification == "2a" |
                           Stages_pdata$AJCC_classification == "2b"] = "Early_stage"
Stages_pdata$Stage_group = as.factor(Stages_pdata$Stage_group)

# Stages_pdata can be used here for phenotypic data
union_stages_z_matrix = union_z_exprs[, Stages_pdata$GEO_accession] # 25699 x 318

# Design and contrast matrices
design = model.matrix(~0 + Stages_pdata$Stage_group + Stages_pdata$Study)
colnames(design) = c("Advanced_stage", "Early_stage", "GSE18670", "GSE21501", 
                      "GSE42952", "GSE62165", "GSE62452", "GSE84219")
rownames(design) = colnames(union_stages_z_matrix)
cont.matrix = makeContrasts(Advancedvsearly=Advanced_stage-Early_stage, levels=design)

# The design2 design matrix can be used here, along with cont.matrix2
union_Stages_z_fit = lmFit(union_stages_z_matrix, design)
union_Stages_z_fit2 = contrasts.fit(union_Stages_z_fit, cont.matrix)
union_Stages_z_fit2 = eBayes(union_Stages_z_fit2, robust = TRUE)
union_Stages_z_results = decideTests(union_Stages_z_fit2)
summary(union_Stages_z_results)
union_Stages_z_DE = as.data.frame(topTable(union_Stages_z_fit2, adjust.method="BH", 
                                           number = Inf))
union_Stages_z_DE$EntrezGene.ID = rownames(union_Stages_z_DE)

# Annotation with official gene symbols
ID_Map$EntrezGene.ID = as.character(ID_Map$EntrezGene.ID)
union_Stages_z_DE_mapped = union_Stages_z_DE %>% left_join(ID_Map, by = "EntrezGene.ID")
union_Stages_z_DE_mapped$Filter = NA
unmapped = which(is.na(union_Stages_z_DE_mapped$HGNC_Official))
union_Stages_z_DE_mapped$HGNC_Official[unmapped] = "unmapped"
for(i in 1:nrow(union_Stages_z_DE_mapped)){
  if(union_Stages_z_DE_mapped$HGNC_Official[i] == "Yes"){
    union_Stages_z_DE_mapped$Filter[i] = "Keep"
  } else if(length(unique(union_Stages_z_DE_mapped$HGNC_Official[union_Stages_z_DE_mapped$EntrezGene.ID ==
                         union_Stages_z_DE_mapped$EntrezGene.ID[i]])) > 1 &&
            union_Stages_z_DE_mapped$HGNC_Official[i] == "No"){
    union_Stages_z_DE_mapped$Filter[i] = "Discard"
  } else if(unique(union_Stages_z_DE_mapped$HGNC_Official[union_Stages_z_DE_mapped$EntrezGene.ID ==
            union_Stages_z_DE_mapped$EntrezGene.ID[i]]) == "No"){
    union_Stages_z_DE_mapped$Filter[i] = "Keep"
    union_Stages_z_DE_mapped$Gene.Symbol[i] = union_Stages_z_DE_mapped$EntrezGene.ID[i]
  } else if(union_Stages_z_DE_mapped$HGNC_Official[i] == "unmapped"){
    union_Stages_z_DE_mapped$Gene.Symbol[i] = union_Stages_z_DE_mapped$EntrezGene.ID[i]
    union_Stages_z_DE_mapped$Filter[i] = "Keep"
  }
}

union_Stages_z_DE_mapped = union_Stages_z_DE_mapped %>% 
  dplyr::filter(Filter == "Keep") %>%
  dplyr::select(EntrezGene.ID, Gene.Symbol, everything()) %>%
  dplyr::select(-Filter) %>%
  distinct()
union_Stages_z_DE_mapped = union_Stages_z_DE_mapped[order(union_Stages_z_DE_mapped$adj.P.Val),]
rownames(union_Stages_z_DE_mapped) = union_Stages_z_DE_mapped$EntrezGene.ID
write.xlsx(union_Stages_z_DE_mapped, "DGEA/Union/Late_stage_vs_Early_stage_z_DE_topTable.xlsx",
           overwrite = TRUE)

# Stage 1 vs stage 4 #####

# Filtered Stages_pdata can be used here for phenotypic data
one_four_pdata = Stages_pdata %>%
  dplyr::filter(AJCC_classification == "1a" |
                  AJCC_classification == "1b" |
                  AJCC_classification == "4") %>%
  dplyr::select(-Stage_group)
one_four_pdata$Stage = NA
one_four_pdata$Stage[one_four_pdata$AJCC_classification == "1a"] = "Stage_one"
one_four_pdata$Stage[one_four_pdata$AJCC_classification == "1b"] = "Stage_one"
one_four_pdata$Stage[one_four_pdata$AJCC_classification == "4"] = "Stage_four"

# Design and contrast matrices
union_one_four_stages_z_matrix = union_z_exprs[, one_four_pdata$GEO_accession] # 25699 x 57
design_1_4 = model.matrix(~0 + one_four_pdata$Stage + one_four_pdata$Study)
colnames(design_1_4) = c("Stage_four", "Stage_one", "GSE18670", "GSE21501", 
                      "GSE42952", "GSE62165", "GSE62452", "GSE84219")
rownames(design_1_4) = colnames(union_one_four_stages_z_matrix)
cont.matrix_1_4 = makeContrasts(Onevsfour = Stage_one - Stage_four, 
                                levels = design_1_4)

union_one_four_Stages_z_fit = lmFit(union_one_four_stages_z_matrix, design_1_4)
union_one_four_Stages_z_fit2 = contrasts.fit(union_one_four_Stages_z_fit, cont.matrix_1_4)
union_one_four_Stages_z_fit2 = eBayes(union_one_four_Stages_z_fit2, robust = TRUE)
union_one_four_Stages_z_results = decideTests(union_one_four_Stages_z_fit2)
summary(union_one_four_Stages_z_results)
union_one_four_Stages_z_DE = as.data.frame(topTable(union_one_four_Stages_z_fit2, 
                                                    adjust.method="BH", number = Inf))
union_one_four_Stages_z_DE$EntrezGene.ID = rownames(union_one_four_Stages_z_DE)

# Annotation with official gene symbols
union_one_four_Stages_z_DE_mapped = union_one_four_Stages_z_DE %>% left_join(ID_Map, by = "EntrezGene.ID")
union_one_four_Stages_z_DE_mapped$Filter = NA
unmapped = which(is.na(union_one_four_Stages_z_DE_mapped$HGNC_Official))
union_one_four_Stages_z_DE_mapped$HGNC_Official[unmapped] = "unmapped"
for(i in 1:nrow(union_one_four_Stages_z_DE_mapped)){
  if(union_one_four_Stages_z_DE_mapped$HGNC_Official[i] == "Yes"){
    union_one_four_Stages_z_DE_mapped$Filter[i] = "Keep"
  } else if(length(unique(union_one_four_Stages_z_DE_mapped$HGNC_Official[union_one_four_Stages_z_DE_mapped$EntrezGene.ID ==
                                                                 union_one_four_Stages_z_DE_mapped$EntrezGene.ID[i]])) > 1 &&
            union_one_four_Stages_z_DE_mapped$HGNC_Official[i] == "No"){
    union_one_four_Stages_z_DE_mapped$Filter[i] = "Discard"
  } else if(unique(union_one_four_Stages_z_DE_mapped$HGNC_Official[union_one_four_Stages_z_DE_mapped$EntrezGene.ID ==
                                                          union_one_four_Stages_z_DE_mapped$EntrezGene.ID[i]]) == "No"){
    union_one_four_Stages_z_DE_mapped$Filter[i] = "Keep"
    union_one_four_Stages_z_DE_mapped$Gene.Symbol[i] = union_one_four_Stages_z_DE_mapped$EntrezGene.ID[i]
  } else if(union_one_four_Stages_z_DE_mapped$HGNC_Official[i] == "unmapped"){
    union_one_four_Stages_z_DE_mapped$Gene.Symbol[i] = union_one_four_Stages_z_DE_mapped$EntrezGene.ID[i]
    union_one_four_Stages_z_DE_mapped$Filter[i] = "Keep"
  }
}

union_one_four_Stages_z_DE_mapped = union_one_four_Stages_z_DE_mapped %>% 
  dplyr::filter(Filter == "Keep") %>%
  dplyr::select(EntrezGene.ID, Gene.Symbol, everything()) %>%
  dplyr::select(-Filter) %>%
  distinct()
union_one_four_Stages_z_DE_mapped = union_one_four_Stages_z_DE_mapped[order(union_one_four_Stages_z_DE_mapped$adj.P.Val),]
rownames(union_one_four_Stages_z_DE_mapped) = union_one_four_Stages_z_DE_mapped$EntrezGene.ID
write.xlsx(union_one_four_Stages_z_DE_mapped, "DGEA/Union/Stage_1_vs_Stage_4_z_DE_topTable.xlsx",
           overwrite = TRUE)

# Stage 1 vs normal #####
one_normal_pdata = full_pdata_filt %>%
  dplyr::filter(AJCC_classification == "1a" |
                  AJCC_classification == "1b" |
                  is.na(AJCC_classification) == TRUE) %>%
  dplyr::filter(!(is.na(AJCC_classification) == TRUE & Tissue_type == "tumor"))
one_normal_z_matrix = z_exprs[, one_normal_pdata$GEO_accession] # 2601 x 165

union_one_normal_z_matrix = union_z_exprs[, one_normal_pdata$GEO_accession] # 25699 x 165

# Design and contrast matrices
design_1_normal = model.matrix(~0 + one_normal_pdata$Tissue_type + one_normal_pdata$Study)
colnames(design_1_normal) = c("non_tumor", "tumor", "GSE18670", "GSE21501", 
                              "GSE42952", "GSE62165", "GSE62452", "GSE84219")
rownames(design_1_normal) = colnames(union_one_normal_z_matrix)
cont.matrix_1_normal = makeContrasts(Onevsnormal = tumor - non_tumor, 
                                     levels = design_1_normal)

union_one_normal_z_fit = lmFit(union_one_normal_z_matrix, design_1_normal)
union_one_normal_z_fit2 = contrasts.fit(union_one_normal_z_fit, cont.matrix_1_normal)
union_one_normal_z_fit2 = eBayes(union_one_normal_z_fit2, robust = TRUE)
union_one_normal_z_results = decideTests(union_one_normal_z_fit2)
summary(union_one_normal_z_results) 

union_one_normal_z_DE = as.data.frame(topTable(union_one_normal_z_fit2, 
                                                      adjust.method="BH", number = Inf))
union_one_normal_z_DE$EntrezGene.ID = rownames(union_one_normal_z_DE)

# 10891 stat. sig. diff. expressed genes of which 1977 (out of 2051) are also found 
# as stat. sig. diff. expressed in the complete-case DGEA:
# length(intersect(union_one_normal_z_DE$EntrezGene.ID[union_one_normal_z_DE$adj.P.Val<0.05],
# one_normal_z_DE$EntrezGene.ID[one_normal_z_DE$adj.P.Val<0.05]))

# Annotation with official gene symbols
union_one_normal_z_DE_mapped = union_one_normal_z_DE %>% left_join(ID_Map, by = "EntrezGene.ID")
union_one_normal_z_DE_mapped$Filter = NA
unmapped = which(is.na(union_one_normal_z_DE_mapped$HGNC_Official))
union_one_normal_z_DE_mapped$HGNC_Official[unmapped] = "unmapped"
for(i in 1:nrow(union_one_normal_z_DE_mapped)){
  if(union_one_normal_z_DE_mapped$HGNC_Official[i] == "Yes"){
    union_one_normal_z_DE_mapped$Filter[i] = "Keep"
  } else if(length(unique(union_one_normal_z_DE_mapped$HGNC_Official[union_one_normal_z_DE_mapped$EntrezGene.ID ==
                                                                            union_one_normal_z_DE_mapped$EntrezGene.ID[i]])) > 1 &&
            union_one_normal_z_DE_mapped$HGNC_Official[i] == "No"){
    union_one_normal_z_DE_mapped$Filter[i] = "Discard"
  } else if(unique(union_one_normal_z_DE_mapped$HGNC_Official[union_one_normal_z_DE_mapped$EntrezGene.ID ==
                                                                     union_one_normal_z_DE_mapped$EntrezGene.ID[i]]) == "No"){
    union_one_normal_z_DE_mapped$Filter[i] = "Keep"
    union_one_normal_z_DE_mapped$Gene.Symbol[i] = union_one_normal_z_DE_mapped$EntrezGene.ID[i]
  } else if(union_one_normal_z_DE_mapped$HGNC_Official[i] == "unmapped"){
    union_one_normal_z_DE_mapped$Gene.Symbol[i] = union_one_normal_z_DE_mapped$EntrezGene.ID[i]
    union_one_normal_z_DE_mapped$Filter[i] = "Keep"
  }
}

union_one_normal_z_DE_mapped = union_one_normal_z_DE_mapped %>% 
  dplyr::filter(Filter == "Keep") %>%
  dplyr::select(EntrezGene.ID, Gene.Symbol, everything()) %>%
  dplyr::select(-Filter) %>%
  distinct()
union_one_normal_z_DE_mapped = union_one_normal_z_DE_mapped[order(union_one_normal_z_DE_mapped$adj.P.Val),]
rownames(union_one_normal_z_DE_mapped) = union_one_normal_z_DE_mapped$EntrezGene.ID
write.xlsx(union_one_normal_z_DE_mapped, "DGEA/Union/Stage_1_vs_Normal_z_DE_topTable.xlsx",
           overwrite = TRUE)

# Stage 2 vs normal #####
two_normal_pdata = full_pdata_filt %>%
  dplyr::filter(AJCC_classification == "2a" |
                  AJCC_classification == "2b" |
                  is.na(AJCC_classification) == TRUE) %>%
  dplyr::filter(!(is.na(AJCC_classification) == TRUE & Tissue_type == "tumor"))
union_two_normal_z_matrix = union_z_exprs[, two_normal_pdata$GEO_accession] # 25699 x 368

# Design and contrast matrices
design_2_normal = model.matrix(~0 + two_normal_pdata$Tissue_type + two_normal_pdata$Study)
colnames(design_2_normal) = c("non_tumor", "tumor", "GSE18670", "GSE21501", 
                              "GSE42952", "GSE62165", "GSE62452", "GSE84219")
rownames(design_2_normal) = colnames(union_two_normal_z_matrix)
cont.matrix_2_normal = makeContrasts(Twovsnormal = tumor - non_tumor, 
                                     levels = design_2_normal)

union_two_normal_z_fit = lmFit(union_two_normal_z_matrix, design_2_normal)
union_two_normal_z_fit2 = contrasts.fit(union_two_normal_z_fit, cont.matrix_2_normal)
union_two_normal_z_fit2 = eBayes(union_two_normal_z_fit2, robust = TRUE)
union_two_normal_z_results = decideTests(union_two_normal_z_fit2)
summary(union_two_normal_z_results) 

union_two_normal_z_DE = as.data.frame(topTable(union_two_normal_z_fit2, 
                                               adjust.method="BH", number = Inf))
union_two_normal_z_DE$EntrezGene.ID = rownames(union_two_normal_z_DE)

# Annotation with official gene symbols
union_two_normal_z_DE_mapped = union_two_normal_z_DE %>% left_join(ID_Map, by = "EntrezGene.ID")
union_two_normal_z_DE_mapped$Filter = NA
unmapped = which(is.na(union_two_normal_z_DE_mapped$HGNC_Official))
union_two_normal_z_DE_mapped$HGNC_Official[unmapped] = "unmapped"
for(i in 1:nrow(union_two_normal_z_DE_mapped)){
  if(union_two_normal_z_DE_mapped$HGNC_Official[i] == "Yes"){
    union_two_normal_z_DE_mapped$Filter[i] = "Keep"
  } else if(length(unique(union_two_normal_z_DE_mapped$HGNC_Official[union_two_normal_z_DE_mapped$EntrezGene.ID ==
                                                                     union_two_normal_z_DE_mapped$EntrezGene.ID[i]])) > 1 &&
            union_two_normal_z_DE_mapped$HGNC_Official[i] == "No"){
    union_two_normal_z_DE_mapped$Filter[i] = "Discard"
  } else if(unique(union_two_normal_z_DE_mapped$HGNC_Official[union_two_normal_z_DE_mapped$EntrezGene.ID ==
                                                              union_two_normal_z_DE_mapped$EntrezGene.ID[i]]) == "No"){
    union_two_normal_z_DE_mapped$Filter[i] = "Keep"
    union_two_normal_z_DE_mapped$Gene.Symbol[i] = union_two_normal_z_DE_mapped$EntrezGene.ID[i]
  } else if(union_two_normal_z_DE_mapped$HGNC_Official[i] == "unmapped"){
    union_two_normal_z_DE_mapped$Gene.Symbol[i] = union_two_normal_z_DE_mapped$EntrezGene.ID[i]
    union_two_normal_z_DE_mapped$Filter[i] = "Keep"
  }
}

union_two_normal_z_DE_mapped = union_two_normal_z_DE_mapped %>% 
  dplyr::filter(Filter == "Keep") %>%
  dplyr::select(EntrezGene.ID, Gene.Symbol, everything()) %>%
  dplyr::select(-Filter) %>%
  distinct()
union_two_normal_z_DE_mapped = union_two_normal_z_DE_mapped[order(union_two_normal_z_DE_mapped$adj.P.Val),]
rownames(union_two_normal_z_DE_mapped) = union_two_normal_z_DE_mapped$EntrezGene.ID
write.xlsx(union_two_normal_z_DE_mapped, "DGEA/Union/Stage_2_vs_Normal_z_DE_topTable.xlsx",
           overwrite = TRUE)

# 16191 stat. sig. diff. expressed genes of which 2285 (out of 2305) are also found as 
# stat. sig. diff. expressed in the complete-case DGEA:
# length(intersect(union_two_normal_z_DE$EntrezGene.ID[union_two_normal_z_DE$adj.P.Val<0.05],
# two_normal_z_DE$EntrezGene.ID[two_normal_z_DE$adj.P.Val<0.05]))

# 16191 stat. sig. diff. expressed genes of which 10621 (out of 10891) are also found as 
# stat. sig. diff. expressed in the Stage 1 vs Normal comparison:
# length(intersect(union_two_normal_z_DE$EntrezGene.ID[union_two_normal_z_DE$adj.P.Val<0.05],
# union_one_normal_z_DE$EntrezGene.ID[union_one_normal_z_DE$adj.P.Val<0.05]))

# Save differences between Stage1/Normal and Stage2/Normal DEGs
# in a variable called "discriminators". However, the real discriminators are both
# genes not in common in the two lists as well as genes which are found to be
# stat. sig. diff. expressed towards opposite directions

common_genes = intersect(union_two_normal_z_DE$EntrezGene.ID[union_two_normal_z_DE$adj.P.Val<0.05],
                         union_one_normal_z_DE$EntrezGene.ID[union_one_normal_z_DE$adj.P.Val<0.05])

diffDEGs = union_one_normal_z_DE_mapped[!union_one_normal_z_DE_mapped$EntrezGene.ID
                                        %in% common_genes, ] %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::select(EntrezGene.ID, Gene.Symbol)

# We need to establish which of the overlapping genes are differentially expressed
# towards the same direction (up-/down-regulated):

stage_1_subset = union_one_normal_z_DE_mapped[union_one_normal_z_DE_mapped$EntrezGene.ID
                                              %in% common_genes, ] %>%
  dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
  dplyr::rename(logFC_stage1 = logFC, adj_p_val_stage1 = adj.P.Val)
stage_2_subset = union_two_normal_z_DE_mapped[union_two_normal_z_DE_mapped$EntrezGene.ID
                                              %in% common_genes, ] %>%
  dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
  dplyr::rename(logFC_stage2 = logFC, adj_p_val_stage2 = adj.P.Val)
common_set = inner_join(stage_1_subset, stage_2_subset, by = c("EntrezGene.ID", "Gene.Symbol"))
common_set$concordance = ifelse(common_set$logFC_stage1*common_set$logFC_stage2 > 0, 1, 0)
concordant_set = common_set %>%
  dplyr::filter(concordance == 1)
discordant_set = common_set %>%
  dplyr::filter(concordance == 0)
concordant_set = concordant_set[order(concordant_set$adj_p_val_stage1, 
                                      concordant_set$adj_p_val_stage2), ]
discordant_set = discordant_set[order(discordant_set$adj_p_val_stage1, 
                                      discordant_set$adj_p_val_stage2), ]
cwb = createWorkbook()
addWorksheet(cwb, "Concordance")
writeData(cwb, "Concordance", concordant_set)
addWorksheet(cwb, "Discordance")
writeData(cwb, "Discordance", discordant_set)
saveWorkbook(cwb, file = "DGEA/Stage_1_Stage_2_union_comparison.xlsx",
             overwrite = TRUE); rm(cwb)

# 10618 genes are stat. sig. diff. expressed towards the same direction between
# the two comparisons. We are interested in the remaining genes which were 
# stat. sig. diff. expressed towards different directions. Just 2 genes are left:
# SNORD115-38, SNORD115-33
discordants = discordant_set %>%
  dplyr::select(EntrezGene.ID, Gene.Symbol)
discriminators = rbind(discordants, diffDEGs)
discriminators = discriminators %>%
  left_join(union_one_normal_z_DE_mapped, by = c("EntrezGene.ID", "Gene.Symbol")) %>%
  dplyr::select(EntrezGene.ID, Gene.Symbol, logFC, adj.P.Val) %>%
  dplyr::rename(logFC_Stage_1 = logFC, adj.P.Val_Stage_1 = adj.P.Val) %>%
  left_join(union_two_normal_z_DE_mapped, by = c("EntrezGene.ID", "Gene.Symbol")) %>%
  dplyr::select(EntrezGene.ID, Gene.Symbol, logFC_Stage_1, adj.P.Val_Stage_1, logFC, adj.P.Val) %>%
  dplyr::rename(logFC_Stage_2 = logFC, adj.P.Val_Stage_2 = adj.P.Val)

write.xlsx(discriminators, "DGEA/Union/discriminators.xlsx", overwrite = TRUE)

# Union volcanoes #####
union_stages_volcano = EnhancedVolcano(union_Stages_z_DE_mapped,
                                   lab = union_Stages_z_DE_mapped[, "Gene.Symbol"],
                                   x = 'logFC',
                                   y = 'adj.P.Val',
                                   title = "Stage 3/4 vs. Stage 1/2",
                                   pCutoff = 0.05,
                                   FCcutoff = 1,
                                   col=c('grey', 'pink', 'purple4', 'red4'),
                                   colAlpha = 0.7)
png("DGEA/Union/Union_Stages_Volcano.png", width = 1920, height = 1080)
union_stages_volcano
dev.off()

union_one_four_stages_volcano = EnhancedVolcano(union_one_four_Stages_z_DE_mapped,
                                       lab = union_one_four_Stages_z_DE_mapped[, "Gene.Symbol"],
                                       x = 'logFC',
                                       y = 'adj.P.Val',
                                       title = "Stage 1 vs. Stage 4",
                                       pCutoff = 0.05,
                                       FCcutoff = 1,
                                       col=c('grey', 'pink', 'purple4', 'red4'),
                                       colAlpha = 0.7)
png("DGEA/Union/Union_One_Four_Stages_Volcano.png", width = 1920, height = 1080)
union_one_four_stages_volcano
dev.off()

union_one_normal_volcano = EnhancedVolcano(union_one_normal_z_DE_mapped,
                                                lab = union_one_normal_z_DE_mapped[, "Gene.Symbol"],
                                                x = 'logFC',
                                                y = 'adj.P.Val',
                                                title = "Stage 1 vs. Normal",
                                                pCutoff = 0.05,
                                                FCcutoff = 1,
                                                col=c('grey', 'pink', 'purple4', 'red4'),
                                                colAlpha = 0.7)
png("DGEA/Union/Union_One_Normal_Volcano.png", width = 1920, height = 1080)
union_one_normal_volcano
dev.off()

union_two_normal_volcano = EnhancedVolcano(union_two_normal_z_DE_mapped,
                                           lab = union_two_normal_z_DE_mapped[, "Gene.Symbol"],
                                           x = 'logFC',
                                           y = 'adj.P.Val',
                                           title = "Stage 2 vs. Normal",
                                           pCutoff = 0.05,
                                           FCcutoff = 1,
                                           col=c('grey', 'pink', 'purple4', 'red4'),
                                           colAlpha = 0.7)
png("DGEA/Union/Union_Two_Normal_Volcano.png", width = 1920, height = 1080)
union_two_normal_volcano
dev.off()

