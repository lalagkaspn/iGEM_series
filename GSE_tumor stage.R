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
# could be yield quite biased results and is avoided in this case. The rows with
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
# BLASTN alignment (against RefSeq RNA sequences, "refseq_rna") was performed
# to annotate probes with RefSeq gene IDs
# https://blast.ncbi.nlm.nih.gov/Blast.cgi
# Read csv results from BLASTN

blastn_1 = read.csv(file = "GSE102238_BLASTN-1.csv", header = FALSE)
blastn_2 = read.csv(file = "GSE102238_BLASTN-2.csv", header = FALSE)
blastn = rbind(blastn_1, blastn_2)
blastn$V2 = gsub("\\..", "", blastn$V2) # RefSeq column: keep main nomenclature - ignore variants
blastn = blastn %>%
  dplyr::select(V1, V2, V3, V11) %>%
  distinct()
rm(blastn_1, blastn_2)
colnames(blastn) = c("probe", "RefSeq", "Alignment_perc", "E_value")
blastn$Alignment_perc = as.numeric(blastn$Alignment_perc)
blastn$E_value = as.numeric(blastn$E_value)
blastn = blastn[order(blastn$probe, 1/blastn$E_value),]
blastn$status = NA

# Add a column to filter for best matches for each probe, in the case of ties,
# keep them all - perhaps it's different RefSeq ID's but the same Entrez ID
for (i in 1:nrow(blastn)){
  if (blastn$Alignment_perc[i] == 
      max(blastn$Alignment_perc[blastn$probe == blastn$probe[i]])) {
    blastn$status[i] = "Keep"
  } else {
    blastn$status[i] = "Discard"
  }
}

blastn = blastn %>%
  dplyr::filter(status == "Keep")

# Map through org.Hs.eg.db
mapped_blastn = inner_join(blastn, ref_df, by = "RefSeq")
mapped_blastn_filt = mapped_blastn %>%
  dplyr::select(probe, ENTREZ_GENE_ID) %>%
  distinct()

# Probes which match to multiple ID's: remove them
mapped_dups = unique(mapped_blastn_filt$probe[which(duplicated(mapped_blastn_filt$probe))])
annot = mapped_blastn_filt[!mapped_blastn_filt$probe %in% mapped_dups,]
rm(blastn, mapped_blastn, mapped_blastn_filt, mapped_dups)

esets[["GSE102238"]]$probe = rownames(esets[["GSE102238"]])
esets[["GSE102238"]] = esets[["GSE102238"]] %>%
  inner_join(annot, by = "probe") %>%
  dplyr::select(-probe) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  summarise_all(mean, na.rm = TRUE)
rm(annot)


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
rownames(original_exprs) = rows; rm(rows) # 2628 x 449
# Making sure we do not have NAs in any row
original_exprs = original_exprs[rowSums(is.na(original_exprs)) != ncol(original_exprs), ]
original_exprs_nonas = na.omit(original_exprs) # 2628 x 449

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
z_exprs = as.matrix(z_exprs %>% dplyr::select(-EntrezGene.ID)) # 2628 x 449
# Making sure we do not have NAs in any row
z_exprs = z_exprs[rowSums(is.na(z_exprs)) != ncol(z_exprs), ]
z_exprs_nonas = na.omit(z_exprs) # 2628 x 449

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
                         as.character(full_pdata$Study), full_pdata$GEO_accession, 
                         as.character(full_pdata$Tissue_type)))
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
                        rep("purple", 100), rep("pink", 30)), outlier.alpha = 0.1)+
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
                        rep("purple", 100), rep("pink", 30)), outlier.alpha = 0.1)+
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
annotation_for_heatmap = full_pdata[, c("Study", "Tissue_type")]
rownames(annotation_for_heatmap) = full_pdata$GEO_accession

original_dists = as.matrix(dist(t(original_exprs_nonas), method = "manhattan"))

rownames(original_dists) = full_pdata$GEO_accession
hmcol = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255))
colnames(original_dists) <- NULL
diag(original_dists) <- NA

ann_colors <- list(
  Type = c(tumor = "deeppink4", non_tumor = "dodgerblue4"),
  Study = c(GSE21501 = "darkseagreen", GSE42952 = "darkorange",
            GSE18670 = "darkcyan", GSE62452 = "darkred",
            GSE62165 = "grey", GSE102238 = "darkmagenta", GSE84219 = "yellow")
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
annotation_for_heatmap = full_pdata[, c("Study", "Tissue_type")]
rownames(annotation_for_heatmap) = full_pdata$GEO_accession

z_dists = as.matrix(dist(t(z_exprs_nonas), method = "manhattan"))

rownames(z_dists) = full_pdata$GEO_accession
hmcol = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255))
colnames(z_dists) <- NULL
diag(z_dists) <- NA

ann_colors <- list(
  Type = c(tumor = "deeppink4", non_tumor = "dodgerblue4"),
  Study = c(GSE21501 = "darkseagreen", GSE42952 = "darkorange",
            GSE18670 = "darkcyan", GSE62452 = "darkred",
            GSE62165 = "grey", GSE102238 = "darkmagenta", GSE84219 = "yellow")
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
# 6 types of DGEA:
#    - Tumor vs Normal
#    - Stages 1/2 vs Stages 3/4
#    - Stage 1 vs Stages 2/3/4
#    - Stage 1 vs Stages 3/4
#    - Stage 2 vs Stages 3/4
#    - Stage 2 vs Stage 1

full_pdata$Study = as.factor(full_pdata$Study)

# Tumor vs Normal #####

# Original matrix
design1 = model.matrix(~0 + full_pdata$Tissue_type + full_pdata$Study)
colnames(design1) = c("non_tumor", "tumor", "GSE18670", "GSE21501", "GSE42952",
                      "GSE62165", "GSE62452", "GSE84219")
rownames(design1) = colnames(original_exprs_nonas)
cont.matrix1 = makeContrasts(tumorvsnontumor=tumor-non_tumor, levels=design1)

TN_fit = lmFit(original_exprs_nonas, design1)
TN_fit2 = contrasts.fit(TN_fit, cont.matrix1)
TN_fit2 = eBayes(TN_fit2, robust = TRUE)
TN_results = decideTests(TN_fit2)
summary(TN_results)
TN_DE = as.data.frame(topTable(TN_fit2, adjust.method ="BH", number = Inf))
TN_DE$EntrezGene.ID = rownames(TN_DE)

# z-score-normalised matrix
TN_z_fit = lmFit(z_exprs_nonas, design1)
TN_z_fit2 = contrasts.fit(TN_z_fit, cont.matrix1)
TN_z_fit2 = eBayes(TN_z_fit2, robust = TRUE)
TN_z_results = decideTests(TN_z_fit2)
summary(TN_z_results)
TN_z_DE = as.data.frame(topTable(TN_z_fit2, adjust.method="BH", number = Inf))
TN_z_DE$EntrezGene.ID = rownames(TN_z_DE)

# Concordance of results (%)
up_concordance = paste0(round(100*length(intersect(TN_z_DE$EntrezGene.ID[TN_z_DE$adj.P.Val < 0.05 &
                                                                           TN_z_DE$logFC > 0],
                                                   TN_DE$EntrezGene.ID[TN_DE$adj.P.Val < 0.05 & TN_DE$logFC > 0]))/
                                max(length(TN_z_DE$EntrezGene.ID[TN_z_DE$adj.P.Val < 0.05 & TN_z_DE$logFC > 0]), 
                                    length(TN_DE$EntrezGene.ID[TN_DE$adj.P.Val < 0.05 & TN_DE$logFC > 0])),2), "%")
down_concordance = paste0(round(100*length(intersect(TN_z_DE$EntrezGene.ID[TN_z_DE$adj.P.Val < 0.05 &
                                                                             TN_z_DE$logFC < 0],
                                                     TN_DE$EntrezGene.ID[TN_DE$adj.P.Val < 0.05 & TN_DE$logFC < 0]))/
                                  max(length(TN_z_DE$EntrezGene.ID[TN_z_DE$adj.P.Val < 0.05 & TN_z_DE$logFC < 0]), 
                                      length(TN_DE$EntrezGene.ID[TN_DE$adj.P.Val < 0.05 & TN_DE$logFC < 0])),2), "%")
ns_concordance = paste0(round(100*length(intersect(TN_z_DE$EntrezGene.ID[TN_z_DE$adj.P.Val > 0.05],
                                                   TN_DE$EntrezGene.ID[TN_DE$adj.P.Val > 0.05]))/
                                max(length(TN_z_DE$EntrezGene.ID[TN_z_DE$adj.P.Val > 0.05]), 
                                    length(TN_DE$EntrezGene.ID[TN_DE$adj.P.Val > 0.05])),2), "%")

# 3 (all) up-DEGs from original matrix are also found as up-DEGs in the z-matrix
# 1 (only) down-DEG from original matrix is also found as down-DEG in the z-matrix

# Annotation with official gene symbols
official = org.Hs.egSYMBOL
mapped_genes_official = mappedkeys(official)
official_df = as.data.frame(official[mapped_genes_official])
official_df = official_df %>% dplyr::rename(EntrezGene.ID = gene_id, Gene.Symbol = symbol)
official_df = official_df[-which(duplicated(official_df$Gene.Symbol)==T),]
official_df = distinct(official_df)

TN_DE_mapped = TN_DE %>% left_join(official_df, by = "EntrezGene.ID")
TN_DE_mapped = TN_DE_mapped %>% dplyr::select(EntrezGene.ID, Gene.Symbol, everything())
rownames(TN_DE_mapped) = TN_DE_mapped$EntrezGene.ID
write.xlsx(TN_DE_mapped, "DGEA/Tumor_stage_analysis/Tumor_vs_Normal/TN_DE_topTable.xlsx")

TN_z_DE_mapped = TN_z_DE %>% left_join(official_df, by = "EntrezGene.ID")
TN_z_DE_mapped = TN_z_DE_mapped %>% dplyr::select(EntrezGene.ID, Gene.Symbol, everything())
rownames(TN_z_DE_mapped) = TN_z_DE_mapped$EntrezGene.ID
write.xlsx(TN_z_DE_mapped, "DGEA/Tumor_stage_analysis/Tumor_vs_Normal/TN_z_DE_topTable.xlsx")

# Stages 1/2 vs stages 3/4 #####

# Filtering the matrices for tumor samples only (the ones with available
# AJCC classification info):
Stages_pdata = full_pdata %>%
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

stages_original_matrix = original_exprs_nonas[, Stages_pdata$GEO_accession] # 2628 x 318
stages_z_matrix = z_exprs_nonas[, Stages_pdata$GEO_accession] # 2628 x 318

# Original matrix
design2 = model.matrix(~0 + Stages_pdata$Stage_group + Stages_pdata$Study)
colnames(design2) = c("Advanced_stage", "Early_stage", "GSE18670", "GSE21501", 
                      "GSE42952", "GSE62165", "GSE62452", "GSE84219")
rownames(design2) = colnames(stages_original_matrix)
cont.matrix2 = makeContrasts(Advancedvsearly=Advanced_stage-Early_stage, levels=design2)

Stages_fit = lmFit(stages_original_matrix, design2)
Stages_fit2 = contrasts.fit(Stages_fit, cont.matrix2)
Stages_fit2 = eBayes(Stages_fit2, robust = TRUE)
Stages_results = decideTests(Stages_fit2)
summary(Stages_results)
Stages_DE = as.data.frame(topTable(Stages_fit2, adjust.method ="BH", number = Inf))
Stages_DE$EntrezGene.ID = rownames(Stages_DE)

# z-score-normalised matrix
Stages_z_fit = lmFit(stages_z_matrix, design2)
Stages_z_fit2 = contrasts.fit(Stages_z_fit, cont.matrix2)
Stages_z_fit2 = eBayes(Stages_z_fit2, robust = TRUE)
Stages_z_results = decideTests(Stages_z_fit2)
summary(Stages_z_results)
Stages_z_DE = as.data.frame(topTable(Stages_z_fit2, adjust.method="BH", number = Inf))
Stages_z_DE$EntrezGene.ID = rownames(Stages_z_DE)

# Concordance of results (%)
stages_up_concordance = paste0(round(100*length(intersect(Stages_z_DE$EntrezGene.ID[Stages_z_DE$adj.P.Val < 0.05 &
                                                                               Stages_z_DE$logFC > 0],
                                                   Stages_DE$EntrezGene.ID[Stages_DE$adj.P.Val < 0.05 & Stages_DE$logFC > 0]))/
                                max(length(Stages_z_DE$EntrezGene.ID[Stages_z_DE$adj.P.Val < 0.05 & Stages_z_DE$logFC > 0]), 
                                    length(Stages_DE$EntrezGene.ID[Stages_DE$adj.P.Val < 0.05 & Stages_DE$logFC > 0])),2), "%")
stages_down_concordance = paste0(round(100*length(intersect(Stages_z_DE$EntrezGene.ID[Stages_z_DE$adj.P.Val < 0.05 &
                                                                                 Stages_z_DE$logFC < 0],
                                                     Stages_DE$EntrezGene.ID[Stages_DE$adj.P.Val < 0.05 & Stages_DE$logFC < 0]))/
                                  max(length(Stages_z_DE$EntrezGene.ID[Stages_z_DE$adj.P.Val < 0.05 & Stages_z_DE$logFC < 0]), 
                                      length(Stages_DE$EntrezGene.ID[Stages_DE$adj.P.Val < 0.05 & Stages_DE$logFC < 0])),2), "%")
stages_ns_concordance = paste0(round(100*length(intersect(Stages_z_DE$EntrezGene.ID[Stages_z_DE$adj.P.Val > 0.05],
                                                   Stages_DE$EntrezGene.ID[Stages_DE$adj.P.Val > 0.05]))/
                                max(length(Stages_z_DE$EntrezGene.ID[Stages_z_DE$adj.P.Val > 0.05]), 
                                    length(Stages_DE$EntrezGene.ID[Stages_DE$adj.P.Val > 0.05])),2), "%")

# No significant results are output by limma on either of the two matrices

# Annotation with official gene symbols
Stages_DE_mapped = Stages_DE %>% left_join(official_df, by = "EntrezGene.ID")
Stages_DE_mapped = Stages_DE_mapped %>% dplyr::select(EntrezGene.ID, Gene.Symbol, everything())
rownames(Stages_DE_mapped) = Stages_DE_mapped$EntrezGene.ID
write.xlsx(Stages_DE_mapped, "DGEA/Tumor_stage_analysis/Late_stage_vs_Early_stage/Stages_DE_topTable.xlsx")

Stages_z_DE_mapped = Stages_z_DE %>% left_join(official_df, by = "EntrezGene.ID")
Stages_z_DE_mapped = Stages_z_DE_mapped %>% dplyr::select(EntrezGene.ID, Gene.Symbol, everything())
rownames(Stages_z_DE_mapped) = Stages_z_DE_mapped$EntrezGene.ID
write.xlsx(Stages_z_DE_mapped, "DGEA/Tumor_stage_analysis/Late_stage_vs_Early_stage/Stages_z_DE_topTable.xlsx")

# Stage 1 vs Stages 2/3/4 #####

# Adding a column for our comparison in Stages_pdata
Stages_pdata$One_vs_all[Stages_pdata$AJCC_classification == "2a" |
                          Stages_pdata$AJCC_classification == "2b" |
                          Stages_pdata$AJCC_classification == 3 | 
                          Stages_pdata$AJCC_classification == 4] = "Advanced_stage"
Stages_pdata$One_vs_all[Stages_pdata$AJCC_classification == "1a" |
                          Stages_pdata$AJCC_classification == "1b"] = "Early_stage"
Stages_pdata$One_vs_all = as.factor(Stages_pdata$One_vs_all)

# Original matrix
design3 = model.matrix(~0 + Stages_pdata$One_vs_all + Stages_pdata$Study)
colnames(design3) = c("Advanced_stage", "Early_stage", "GSE18670", "GSE21501", 
                      "GSE42952", "GSE62165", "GSE62452", "GSE84219")
rownames(design3) = colnames(stages_original_matrix)
cont.matrix3 = makeContrasts(Advancedvsearly=Advanced_stage-Early_stage, levels=design3)

Onevsall_fit = lmFit(stages_original_matrix, design3)
Onevsall_fit2 = contrasts.fit(Onevsall_fit, cont.matrix3)
Onevsall_fit2 = eBayes(Onevsall_fit2, robust = TRUE)
Onevsall_results = decideTests(Onevsall_fit2)
summary(Onevsall_results)
Onevsall_DE = as.data.frame(topTable(Onevsall_fit2, adjust.method ="BH", number = Inf))
Onevsall_DE$EntrezGene.ID = rownames(Onevsall_DE)

# z-score-normalised matrix
Onevsall_z_fit = lmFit(stages_z_matrix, design3)
Onevsall_z_fit2 = contrasts.fit(Onevsall_z_fit, cont.matrix3)
Onevsall_z_fit2 = eBayes(Onevsall_z_fit2, robust = TRUE)
Onevsall_z_results = decideTests(Onevsall_z_fit2)
summary(Onevsall_z_results)
Onevsall_z_DE = as.data.frame(topTable(Onevsall_z_fit2, adjust.method="BH", number = Inf))
Onevsall_z_DE$EntrezGene.ID = rownames(Onevsall_z_DE)

# Concordance of results (%)
Onevsall_up_concordance = paste0(round(100*length(intersect(Onevsall_z_DE$EntrezGene.ID[Onevsall_z_DE$adj.P.Val < 0.05 &
                                                                                          Onevsall_z_DE$logFC > 0],
                                                            Onevsall_DE$EntrezGene.ID[Onevsall_DE$adj.P.Val < 0.05 & Onevsall_DE$logFC > 0]))/
                                         max(length(Onevsall_z_DE$EntrezGene.ID[Onevsall_z_DE$adj.P.Val < 0.05 & Onevsall_z_DE$logFC > 0]), 
                                             length(Onevsall_DE$EntrezGene.ID[Onevsall_DE$adj.P.Val < 0.05 & Onevsall_DE$logFC > 0])),2), "%")
Onevsall_down_concordance = paste0(round(100*length(intersect(Onevsall_z_DE$EntrezGene.ID[Onevsall_z_DE$adj.P.Val < 0.05 &
                                                                                            Onevsall_z_DE$logFC < 0],
                                                              Onevsall_DE$EntrezGene.ID[Onevsall_DE$adj.P.Val < 0.05 & Onevsall_DE$logFC < 0]))/
                                           max(length(Onevsall_z_DE$EntrezGene.ID[Onevsall_z_DE$adj.P.Val < 0.05 & Onevsall_z_DE$logFC < 0]), 
                                               length(Onevsall_DE$EntrezGene.ID[Onevsall_DE$adj.P.Val < 0.05 & Onevsall_DE$logFC < 0])),2), "%")
Onevsall_ns_concordance = paste0(round(100*length(intersect(Onevsall_z_DE$EntrezGene.ID[Onevsall_z_DE$adj.P.Val > 0.05],
                                                            Onevsall_DE$EntrezGene.ID[Onevsall_DE$adj.P.Val > 0.05]))/
                                         max(length(Onevsall_z_DE$EntrezGene.ID[Onevsall_z_DE$adj.P.Val > 0.05]), 
                                             length(Onevsall_DE$EntrezGene.ID[Onevsall_DE$adj.P.Val > 0.05])),2), "%")

# No significant results are output by limma on the z-score-transformed matrix

# Annotation with official gene symbols
Onevsall_DE_mapped = Onevsall_DE %>% left_join(official_df, by = "EntrezGene.ID")
Onevsall_DE_mapped = Onevsall_DE_mapped %>% dplyr::select(EntrezGene.ID, Gene.Symbol, everything())
rownames(Onevsall_DE_mapped) = Onevsall_DE_mapped$EntrezGene.ID
write.xlsx(Onevsall_DE_mapped, "DGEA/Tumor_stage_analysis/Late_stage_vs_Early_stage/Onevsall_DE_topTable.xlsx")

Onevsall_z_DE_mapped = Onevsall_z_DE %>% left_join(official_df, by = "EntrezGene.ID")
Onevsall_z_DE_mapped = Onevsall_z_DE_mapped %>% dplyr::select(EntrezGene.ID, Gene.Symbol, everything())
rownames(Onevsall_z_DE_mapped) = Onevsall_z_DE_mapped$EntrezGene.ID
write.xlsx(Onevsall_z_DE_mapped, "DGEA/Tumor_stage_analysis/Late_stage_vs_Early_stage/Onevsall_z_DE_topTable.xlsx")
# Stage 1 vs Stages 3/4 #####

# Removing stage 2 samples
Stages_pdata_comp4 = Stages_pdata %>%
  dplyr::filter(AJCC_classification != "2a") %>%
  dplyr::filter(AJCC_classification != "2b")
# 80 x 6 pdata frame 

# Original matrix
design4 = model.matrix(~0 + Stages_pdata_comp4$One_vs_all + Stages_pdata_comp4$Study)
colnames(design4) = c("Advanced_stage", "Early_stage", "GSE18670", "GSE21501", 
                      "GSE42952", "GSE62165", "GSE62452", "GSE84219")

comp4_matrix = original_exprs_nonas[, Stages_pdata_comp4$GEO_accession]
z_comp4_matrix = z_exprs_nonas[, Stages_pdata_comp4$GEO_accession]

rownames(design4) = colnames(comp4_matrix)
cont.matrix4 = makeContrasts(Advancedvsearly=Advanced_stage-Early_stage, levels=design4)

Onevslate_fit = lmFit(comp4_matrix, design4)
Onevslate_fit2 = contrasts.fit(Onevslate_fit, cont.matrix4)
Onevslate_fit2 = eBayes(Onevslate_fit2, robust = TRUE)
Onevslate_results = decideTests(Onevslate_fit2)
summary(Onevslate_results)
Onevslate_DE = as.data.frame(topTable(Onevslate_fit2, adjust.method ="BH", number = Inf))
Onevslate_DE$EntrezGene.ID = rownames(Onevslate_DE)

# z-score-normalised matrix
Onevslate_z_fit = lmFit(z_comp4_matrix, design4)
Onevslate_z_fit2 = contrasts.fit(Onevslate_z_fit, cont.matrix4)
Onevslate_z_fit2 = eBayes(Onevslate_z_fit2, robust = TRUE)
Onevslate_z_results = decideTests(Onevslate_z_fit2)
summary(Onevslate_z_results)
Onevslate_z_DE = as.data.frame(topTable(Onevslate_z_fit2, adjust.method="BH", number = Inf))
Onevslate_z_DE$EntrezGene.ID = rownames(Onevslate_z_DE)

# Concordance of results (%)
Onevslate_up_concordance = paste0(round(100*length(intersect(Onevslate_z_DE$EntrezGene.ID[Onevslate_z_DE$adj.P.Val < 0.05 &
                                                                                            Onevslate_z_DE$logFC > 0],
                                                             Onevslate_DE$EntrezGene.ID[Onevslate_DE$adj.P.Val < 0.05 & Onevslate_DE$logFC > 0]))/
                                          max(length(Onevslate_z_DE$EntrezGene.ID[Onevslate_z_DE$adj.P.Val < 0.05 & Onevslate_z_DE$logFC > 0]), 
                                              length(Onevslate_DE$EntrezGene.ID[Onevslate_DE$adj.P.Val < 0.05 & Onevslate_DE$logFC > 0])),2), "%")
Onevslate_down_concordance = paste0(round(100*length(intersect(Onevslate_z_DE$EntrezGene.ID[Onevslate_z_DE$adj.P.Val < 0.05 &
                                                                                              Onevslate_z_DE$logFC < 0],
                                                               Onevslate_DE$EntrezGene.ID[Onevslate_DE$adj.P.Val < 0.05 & Onevslate_DE$logFC < 0]))/
                                            max(length(Onevslate_z_DE$EntrezGene.ID[Onevslate_z_DE$adj.P.Val < 0.05 & Onevslate_z_DE$logFC < 0]), 
                                                length(Onevslate_DE$EntrezGene.ID[Onevslate_DE$adj.P.Val < 0.05 & Onevslate_DE$logFC < 0])),2), "%")
Onevslate_ns_concordance = paste0(round(100*length(intersect(Onevslate_z_DE$EntrezGene.ID[Onevslate_z_DE$adj.P.Val > 0.05],
                                                             Onevslate_DE$EntrezGene.ID[Onevslate_DE$adj.P.Val > 0.05]))/
                                          max(length(Onevslate_z_DE$EntrezGene.ID[Onevslate_z_DE$adj.P.Val > 0.05]), 
                                              length(Onevslate_DE$EntrezGene.ID[Onevslate_DE$adj.P.Val > 0.05])),2), "%")

# No significant results are output by limma on either matrix

# Annotation with official gene symbols
Onevslate_DE_mapped = Onevslate_DE %>% left_join(official_df, by = "EntrezGene.ID")
Onevslate_DE_mapped = Onevslate_DE_mapped %>% dplyr::select(EntrezGene.ID, Gene.Symbol, everything())
rownames(Onevslate_DE_mapped) = Onevslate_DE_mapped$EntrezGene.ID
write.xlsx(Onevslate_DE_mapped, "DGEA/Tumor_stage_analysis/Late_stage_vs_Early_stage/Onevslate_DE_topTable.xlsx")

Onevslate_z_DE_mapped = Onevslate_z_DE %>% left_join(official_df, by = "EntrezGene.ID")
Onevslate_z_DE_mapped = Onevslate_z_DE_mapped %>% dplyr::select(EntrezGene.ID, Gene.Symbol, everything())
rownames(Onevslate_z_DE_mapped) = Onevslate_z_DE_mapped$EntrezGene.ID
write.xlsx(Onevslate_z_DE_mapped, "DGEA/Tumor_stage_analysis/Late_stage_vs_Early_stage/Onevslate_z_DE_topTable.xlsx")
# Stage 2 vs Stages 3/4 #####

# Removing stage 1 samples
Stages_pdata_comp5 = Stages_pdata %>%
  dplyr::filter(AJCC_classification != "1a") %>%
  dplyr::filter(AJCC_classification != "1b")
# 283 x 6 pdata frame 

# Original matrix
design5 = model.matrix(~0 + Stages_pdata_comp5$Stage_group + Stages_pdata_comp5$Study)
colnames(design5) = c("Advanced_stage", "Early_stage", "GSE18670", "GSE21501", 
                      "GSE42952", "GSE62165", "GSE62452", "GSE84219")

comp5_matrix = original_exprs_nonas[, Stages_pdata_comp5$GEO_accession]
z_comp5_matrix = z_exprs_nonas[, Stages_pdata_comp5$GEO_accession]

rownames(design5) = colnames(comp5_matrix)
cont.matrix5 = makeContrasts(Advancedvsearly=Advanced_stage-Early_stage, levels=design5)

Twovslate_fit = lmFit(comp5_matrix, design5)
Twovslate_fit2 = contrasts.fit(Twovslate_fit, cont.matrix5)
Twovslate_fit2 = eBayes(Twovslate_fit2, robust = TRUE)
Twovslate_results = decideTests(Twovslate_fit2)
summary(Twovslate_results)
Twovslate_DE = as.data.frame(topTable(Twovslate_fit2, adjust.method ="BH", number = Inf))
Twovslate_DE$EntrezGene.ID = rownames(Twovslate_DE)

# z-score-normalised matrix
Twovslate_z_fit = lmFit(z_comp5_matrix, design5)
Twovslate_z_fit2 = contrasts.fit(Twovslate_z_fit, cont.matrix5)
Twovslate_z_fit2 = eBayes(Twovslate_z_fit2, robust = TRUE)
Twovslate_z_results = decideTests(Twovslate_z_fit2)
summary(Twovslate_z_results)
Twovslate_z_DE = as.data.frame(topTable(Twovslate_z_fit2, adjust.method="BH", number = Inf))
Twovslate_z_DE$EntrezGene.ID = rownames(Twovslate_z_DE)

# Concordance of results (%)
Twovslate_up_concordance = paste0(round(100*length(intersect(Twovslate_z_DE$EntrezGene.ID[Twovslate_z_DE$adj.P.Val < 0.05 &
                                                                                            Twovslate_z_DE$logFC > 0],
                                                             Twovslate_DE$EntrezGene.ID[Twovslate_DE$adj.P.Val < 0.05 & Twovslate_DE$logFC > 0]))/
                                          max(length(Twovslate_z_DE$EntrezGene.ID[Twovslate_z_DE$adj.P.Val < 0.05 & Twovslate_z_DE$logFC > 0]), 
                                              length(Twovslate_DE$EntrezGene.ID[Twovslate_DE$adj.P.Val < 0.05 & Twovslate_DE$logFC > 0])),2), "%")
Twovslate_down_concordance = paste0(round(100*length(intersect(Twovslate_z_DE$EntrezGene.ID[Twovslate_z_DE$adj.P.Val < 0.05 &
                                                                                              Twovslate_z_DE$logFC < 0],
                                                               Twovslate_DE$EntrezGene.ID[Twovslate_DE$adj.P.Val < 0.05 & Twovslate_DE$logFC < 0]))/
                                            max(length(Twovslate_z_DE$EntrezGene.ID[Twovslate_z_DE$adj.P.Val < 0.05 & Twovslate_z_DE$logFC < 0]), 
                                                length(Twovslate_DE$EntrezGene.ID[Twovslate_DE$adj.P.Val < 0.05 & Twovslate_DE$logFC < 0])),2), "%")
Twovslate_ns_concordance = paste0(round(100*length(intersect(Twovslate_z_DE$EntrezGene.ID[Twovslate_z_DE$adj.P.Val > 0.05],
                                                             Twovslate_DE$EntrezGene.ID[Twovslate_DE$adj.P.Val > 0.05]))/
                                          max(length(Twovslate_z_DE$EntrezGene.ID[Twovslate_z_DE$adj.P.Val > 0.05]), 
                                              length(Twovslate_DE$EntrezGene.ID[Twovslate_DE$adj.P.Val > 0.05])),2), "%")

# No significant results are output by limma on either matrix

# Annotation with official gene symbols
Twovslate_DE_mapped = Twovslate_DE %>% left_join(official_df, by = "EntrezGene.ID")
Twovslate_DE_mapped = Twovslate_DE_mapped %>% dplyr::select(EntrezGene.ID, Gene.Symbol, everything())
rownames(Twovslate_DE_mapped) = Twovslate_DE_mapped$EntrezGene.ID
write.xlsx(Twovslate_DE_mapped, "DGEA/Tumor_stage_analysis/Late_stage_vs_Early_stage/Twovslate_DE_topTable.xlsx")

Twovslate_z_DE_mapped = Twovslate_z_DE %>% left_join(official_df, by = "EntrezGene.ID")
Twovslate_z_DE_mapped = Twovslate_z_DE_mapped %>% dplyr::select(EntrezGene.ID, Gene.Symbol, everything())
rownames(Twovslate_z_DE_mapped) = Twovslate_z_DE_mapped$EntrezGene.ID
write.xlsx(Twovslate_z_DE_mapped, "DGEA/Tumor_stage_analysis/Late_stage_vs_Early_stage/Twovslate_z_DE_topTable.xlsx")

# Stage 2 vs Stage 1 #####

# Removing stage 3/4 samples
Stages_pdata_comp6 = Stages_pdata %>%
  dplyr::filter(AJCC_classification != 3) %>%
  dplyr::filter(AJCC_classification != 4)
Stages_pdata_comp6$Stage = NA
Stages_pdata_comp6$Stage[Stages_pdata_comp6$AJCC_classification == "1a" |
                           Stages_pdata_comp6$AJCC_classification == "1b"] = "Stage_1"
Stages_pdata_comp6$Stage[Stages_pdata_comp6$AJCC_classification == "2a" |
                           Stages_pdata_comp6$AJCC_classification == "2b"] = "Stage_2"
Stages_pdata_comp6$Stage = as.factor(Stages_pdata_comp6$Stage)
# 273 x 7 pdata frame 

# Original matrix
design6 = model.matrix(~0 + Stages_pdata_comp6$Stage + Stages_pdata_comp6$Study)
colnames(design6) = c("Stage_1", "Stage_2", "GSE18670", "GSE21501", 
                      "GSE42952", "GSE62165", "GSE62452", "GSE84219")

comp6_matrix = original_exprs_nonas[, Stages_pdata_comp6$GEO_accession]
z_comp6_matrix = z_exprs_nonas[, Stages_pdata_comp6$GEO_accession]

rownames(design6) = colnames(comp6_matrix)
cont.matrix6 = makeContrasts(TwovsOne=Stage_2-Stage_1, levels=design6)

Twovsone_fit = lmFit(comp6_matrix, design6)
Twovsone_fit2 = contrasts.fit(Twovsone_fit, cont.matrix6)
Twovsone_fit2 = eBayes(Twovsone_fit2, robust = TRUE)
Twovsone_results = decideTests(Twovsone_fit2)
summary(Twovsone_results)
Twovsone_DE = as.data.frame(topTable(Twovsone_fit2, adjust.method ="BH", number = Inf))
Twovsone_DE$EntrezGene.ID = rownames(Twovsone_DE)

# z-score-normalised matrix
Twovsone_z_fit = lmFit(z_comp6_matrix, design6)
Twovsone_z_fit2 = contrasts.fit(Twovsone_z_fit, cont.matrix6)
Twovsone_z_fit2 = eBayes(Twovsone_z_fit2, robust = TRUE)
Twovsone_z_results = decideTests(Twovsone_z_fit2)
summary(Twovsone_z_results)
Twovsone_z_DE = as.data.frame(topTable(Twovsone_z_fit2, adjust.method="BH", number = Inf))
Twovsone_z_DE$EntrezGene.ID = rownames(Twovsone_z_DE)

# Concordance of results (%)
Twovsone_up_concordance = paste0(round(100*length(intersect(Twovsone_z_DE$EntrezGene.ID[Twovsone_z_DE$adj.P.Val < 0.05 &
                                                                                          Twovsone_z_DE$logFC > 0],
                                                            Twovsone_DE$EntrezGene.ID[Twovsone_DE$adj.P.Val < 0.05 & Twovsone_DE$logFC > 0]))/
                                         max(length(Twovsone_z_DE$EntrezGene.ID[Twovsone_z_DE$adj.P.Val < 0.05 & Twovsone_z_DE$logFC > 0]), 
                                             length(Twovsone_DE$EntrezGene.ID[Twovsone_DE$adj.P.Val < 0.05 & Twovsone_DE$logFC > 0])),2), "%")
Twovsone_down_concordance = paste0(round(100*length(intersect(Twovsone_z_DE$EntrezGene.ID[Twovsone_z_DE$adj.P.Val < 0.05 &
                                                                                            Twovsone_z_DE$logFC < 0],
                                                              Twovsone_DE$EntrezGene.ID[Twovsone_DE$adj.P.Val < 0.05 & Twovsone_DE$logFC < 0]))/
                                           max(length(Twovsone_z_DE$EntrezGene.ID[Twovsone_z_DE$adj.P.Val < 0.05 & Twovsone_z_DE$logFC < 0]), 
                                               length(Twovsone_DE$EntrezGene.ID[Twovsone_DE$adj.P.Val < 0.05 & Twovsone_DE$logFC < 0])),2), "%")
Twovsone_ns_concordance = paste0(round(100*length(intersect(Twovsone_z_DE$EntrezGene.ID[Twovsone_z_DE$adj.P.Val > 0.05],
                                                            Twovsone_DE$EntrezGene.ID[Twovsone_DE$adj.P.Val > 0.05]))/
                                         max(length(Twovsone_z_DE$EntrezGene.ID[Twovsone_z_DE$adj.P.Val > 0.05]), 
                                             length(Twovsone_DE$EntrezGene.ID[Twovsone_DE$adj.P.Val > 0.05])),2), "%")

# No significant results are output by limma on the z-score-normalised matrix

# Annotation with official gene symbols
Twovsone_DE_mapped = Twovsone_DE %>% left_join(official_df, by = "EntrezGene.ID")
Twovsone_DE_mapped = Twovsone_DE_mapped %>% dplyr::select(EntrezGene.ID, Gene.Symbol, everything())
rownames(Twovsone_DE_mapped) = Twovsone_DE_mapped$EntrezGene.ID
write.xlsx(Twovsone_DE_mapped, "DGEA/Tumor_stage_analysis/Late_stage_vs_Early_stage/Twovsone_DE_topTable.xlsx")

Twovsone_z_DE_mapped = Twovsone_z_DE %>% left_join(official_df, by = "EntrezGene.ID")
Twovsone_z_DE_mapped = Twovsone_z_DE_mapped %>% dplyr::select(EntrezGene.ID, Gene.Symbol, everything())
rownames(Twovsone_z_DE_mapped) = Twovsone_z_DE_mapped$EntrezGene.ID
write.xlsx(Twovsone_z_DE_mapped, "DGEA/Tumor_stage_analysis/Late_stage_vs_Early_stage/Twovsone_z_DE_topTable.xlsx")

##### Volcano plots #####
# Tumor vs Normal
TN_volcano = EnhancedVolcano(TN_DE_mapped,
                             lab = TN_DE_mapped[, "Gene.Symbol"],
                             x = 'logFC',
                             y = 'adj.P.Val',
                             title = "Tumor vs. Non-tumor",
                             pCutoff = 0.05,
                             FCcutoff = 2,
                             col=c('grey', 'pink', 'purple4', 'red4'),
                             colAlpha = 0.7)
png("DGEA/Tumor_stage_analysis/Tumor_vs_Normal/TN_Volcano.png", width = 1920, height = 1080)
TN_volcano
dev.off()

TN_z_volcano = EnhancedVolcano(TN_z_DE_mapped,
                               lab = TN_z_DE_mapped[, "Gene.Symbol"],
                               x = 'logFC',
                               y = 'adj.P.Val',
                               title = "Tumor vs. Non-tumor (z-normalised)",
                               pCutoff = 0.05,
                               FCcutoff = 0.5,
                               col=c('grey', 'pink', 'purple4', 'red4'),
                               colAlpha = 0.7)
png("DGEA/Tumor_stage_analysis/Tumor_vs_Normal/TN_z_Volcano.png", width = 1920, height = 1080)
TN_z_volcano
dev.off()

# Stages 1/2 vs Stages 3/4
Stages_volcano = EnhancedVolcano(Stages_DE_mapped,
                                 lab = Stages_DE_mapped[, "Gene.Symbol"],
                                 x = 'logFC',
                                 y = 'adj.P.Val',
                                 title = "Stages 3/4 vs. Stages 1/2",
                                 pCutoff = 0.05,
                                 FCcutoff = 2,
                                 col=c('grey', 'pink', 'purple4', 'red4'),
                                 colAlpha = 0.7)
png("DGEA/Tumor_stage_analysis/Late_stage_vs_Early_stage/Stages_Volcano.png", width = 1920, height = 1080)
Stages_volcano
dev.off()

Stages_z_volcano = EnhancedVolcano(Stages_z_DE_mapped,
                                   lab = Stages_z_DE_mapped[, "Gene.Symbol"],
                                   x = 'logFC',
                                   y = 'adj.P.Val',
                                   title = "Stages 3/4 vs. Stages 1/2 (z-normalised)",
                                   pCutoff = 0.05,
                                   FCcutoff = 0.5,
                                   col=c('grey', 'pink', 'purple4', 'red4'),
                                   colAlpha = 0.7)
png("DGEA/Tumor_stage_analysis/Late_stage_vs_Early_stage/Stages_z_Volcano.png", width = 1920, height = 1080)
Stages_z_volcano
dev.off()

# Stage 1 vs Stages 2/3/4
Onevsall_volcano = EnhancedVolcano(Onevsall_DE_mapped,
                                   lab = Onevsall_DE_mapped[, "Gene.Symbol"],
                                   x = 'logFC',
                                   y = 'adj.P.Val',
                                   title = "Stages 2/3/4 vs. Stage 1",
                                   pCutoff = 0.05,
                                   FCcutoff = 2,
                                   col=c('grey', 'pink', 'purple4', 'red4'),
                                   colAlpha = 0.7)
png("DGEA/Tumor_stage_analysis/Late_stage_vs_Early_stage/Onevsall_Volcano.png", width = 1920, height = 1080)
Onevsall_volcano
dev.off()

Onevsall_z_volcano = EnhancedVolcano(Onevsall_z_DE_mapped,
                                     lab = Onevsall_z_DE_mapped[, "Gene.Symbol"],
                                     x = 'logFC',
                                     y = 'adj.P.Val',
                                     title = "Stages 2/3/4 vs. Stage 1 (z-normalised)",
                                     pCutoff = 0.05,
                                     FCcutoff = 0.5,
                                     col=c('grey', 'pink', 'purple4', 'red4'),
                                     colAlpha = 0.7)
png("DGEA/Tumor_stage_analysis/Late_stage_vs_Early_stage/Onevsall_z_Volcano.png", width = 1920, height = 1080)
Onevsall_z_volcano
dev.off()

# Stage 1 vs Stages 3/4
Onevslate_volcano = EnhancedVolcano(Onevslate_DE_mapped,
                                    lab = Onevslate_DE_mapped[, "Gene.Symbol"],
                                    x = 'logFC',
                                    y = 'adj.P.Val',
                                    title = "Stages 3/4 vs. Stage 1",
                                    pCutoff = 0.05,
                                    FCcutoff = 2,
                                    col=c('grey', 'pink', 'purple4', 'red4'),
                                    colAlpha = 0.7)
png("DGEA/Tumor_stage_analysis/Late_stage_vs_Early_stage/Onevslate_Volcano.png", width = 1920, height = 1080)
Onevslate_volcano
dev.off()

Onevslate_z_volcano = EnhancedVolcano(Onevslate_z_DE_mapped,
                                      lab = Onevslate_z_DE_mapped[, "Gene.Symbol"],
                                      x = 'logFC',
                                      y = 'adj.P.Val',
                                      title = "Stages 3/4 vs. Stage 1 (z-normalised)",
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      col=c('grey', 'pink', 'purple4', 'red4'),
                                      colAlpha = 0.7)
png("DGEA/Tumor_stage_analysis/Late_stage_vs_Early_stage/Onevslate_z_Volcano.png", width = 1920, height = 1080)
Onevslate_z_volcano
dev.off()

# Stage 2 vs Stages 3/4
Twovslate_volcano = EnhancedVolcano(Twovslate_DE_mapped,
                                    lab = Twovslate_DE_mapped[, "Gene.Symbol"],
                                    x = 'logFC',
                                    y = 'adj.P.Val',
                                    title = "Stages 3/4 vs. Stage 2",
                                    pCutoff = 0.05,
                                    FCcutoff = 2,
                                    col=c('grey', 'pink', 'purple4', 'red4'),
                                    colAlpha = 0.7)
png("DGEA/Tumor_stage_analysis/Late_stage_vs_Early_stage/Twovslate_Volcano.png", width = 1920, height = 1080)
Twovslate_volcano
dev.off()

Twovslate_z_volcano = EnhancedVolcano(Twovslate_z_DE_mapped,
                                      lab = Twovslate_z_DE_mapped[, "Gene.Symbol"],
                                      x = 'logFC',
                                      y = 'adj.P.Val',
                                      title = "Stages 3/4 vs. Stage 2 (z-normalised)",
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      col=c('grey', 'pink', 'purple4', 'red4'),
                                      colAlpha = 0.7)
png("DGEA/Tumor_stage_analysis/Late_stage_vs_Early_stage/Twovslate_z_Volcano.png", width = 1920, height = 1080)
Twovslate_z_volcano
dev.off()

# Stage 2 vs Stage 1
Twovsone_volcano = EnhancedVolcano(Twovsone_DE_mapped,
                                   lab = Twovsone_DE_mapped[, "Gene.Symbol"],
                                   x = 'logFC',
                                   y = 'adj.P.Val',
                                   title = "Stage 2 vs. Stage 1",
                                   pCutoff = 0.05,
                                   FCcutoff = 2,
                                   col=c('grey', 'pink', 'purple4', 'red4'),
                                   colAlpha = 0.7)
png("DGEA/Tumor_stage_analysis/Late_stage_vs_Early_stage/Twovsone_Volcano.png", width = 1920, height = 1080)
Twovsone_volcano
dev.off()

Twovsone_z_volcano = EnhancedVolcano(Twovsone_z_DE_mapped,
                                     lab = Twovsone_z_DE_mapped[, "Gene.Symbol"],
                                     x = 'logFC',
                                     y = 'adj.P.Val',
                                     title = "Stage 2 vs. Stage 1 (z-normalised)",
                                     pCutoff = 0.05,
                                     FCcutoff = 0.5,
                                     col=c('grey', 'pink', 'purple4', 'red4'),
                                     colAlpha = 0.7)
png("DGEA/Tumor_stage_analysis/Late_stage_vs_Early_stage/Twovsone_z_Volcano.png", width = 1920, height = 1080)
Twovsone_z_volcano
dev.off()

# Writing out the z-score-normalised expression matrix for machine learning purposes
tumor_matrix_pre = as.data.frame(z_exprs_nonas)
tumor_matrix_pre$EntrezGene.ID = rownames(z_exprs_nonas)
tumor_matrix_pre = tumor_matrix_pre %>%
  inner_join(TN_z_DE_mapped, by = "EntrezGene.ID") %>%
  dplyr::select(-logFC, -AveExpr, -t, -P.Value, -adj.P.Val, -B, -EntrezGene.ID)
rownames(tumor_matrix_pre) = tumor_matrix_pre$Gene.Symbol
tumor_matrix = t(tumor_matrix_pre)
tumor_frame = as.data.frame(tumor_matrix) %>%
  mutate(GEO_accession = rownames(tumor_matrix)) %>%
  inner_join(full_pdata, by = "GEO_accession")
write.xlsx(tumor_frame, "Tumor_samples_z_expression_matrix.xlsx")
rm(tumor_matrix, tumor_matrix_pre)
