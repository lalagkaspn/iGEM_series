## This script is used to download and pre-process series (GSE) from GEO for PDAC patients (normal and tumor tissues) who did not take
## neoadjuvant chemotherapy and mention tumor stage information.

library(dplyr)
library(readr)
library(GEOquery)
library(org.Hs.eg.db)

##### Downloading data #####
datasets = c("GSE21501", "GSE42952", "GSE18670", "GSE62452", "GSE62165", "GSE102238", "GSE84219")

# Run this before getGEO
# Probably, there is a bug with the newest version of readr and GEOquery
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
#   N0: there are no cancer cells in the nearby lymph
#   N1: there are 1 to 3 lymph nodes that contain cancer cells
#   N2: there is cancer in 4 or more lymph nodes
# M: Describes whether the cancer has spread to a different part of the body
#   M0: cancer has not spread to other areas of the body
#   M1: cancer has spread to other areas of the body

# Association of AJCC classification (8th edition) and TNM classification for pancreatic cancer
# https://www.nature.com/articles/s41598-018-28193-4/tables/1

# Tumor grade for PDAC (this classification was used in all GSEs except GSE62452)
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

# Filter pdata for keeping only needed information
filt_pdata = list()

# GSE21501
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21501
# Paper: https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC2903589/
# Comments:
#   - Clinical data from "University of Nebraska Medical Center Rapid Autopsy Pancreatic Program" (NEB) and "Yeh Lab-UNC-CH" (UNC) were 
#     not used in the original analysis. Therefore, authors did not include them in the GSE dataset.
#   - GSM536033: it was used as reference and will be excluded from further analysis.
#   - 97 samples were primary PDAC and 15 were metastatic PDAC. For this reason, I will add an extra row for M_stage. 
#     However, for metastatic PDAC patients, they took samples from primary tumor.
#   - From "UNC" samples, 1 patient received neodjuvant chemotherapy and 4 received adjuvant chemotherapy. NO GSM-specific information.
#     All UNC samples will be removed!
#   - From "NEB", 11 patients received chemotherapy <6 months prior to death. NO GSM-specific information.
#   - From "Northwestern Memorial Hospital" (NW) and "NorthShore University HealthSystem" (NSU) samples, 2 received neoadjuvant chemotherapy 
#     and 30 received adjuvant chemotherapy. NO GSM-specific information.
#     All NW and NSU patient samples will be removed!

# Select pdata
filt_pdata[["GSE21501"]] = pdata$GSE21501 %>%
  select(GEO_accession = geo_accession, 
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
filt_pdata$GSE21501$Biomaterial_provider = ifelse(filt_pdata$GSE21501$Biomaterial_provider == "University of Nebraska Medical Center Rapid Autopsy Pancreatic Program", "NEB", 
                                                  ifelse(filt_pdata$GSE21501$Biomaterial_provider == "Yeh Lab-UNC-CH", "UNC",
                                                         ifelse(filt_pdata$GSE21501$Biomaterial_provider == "Northwestern Memorial Hospital", "NW",
                                                                ifelse(filt_pdata$GSE21501$Biomaterial_provider == "NorthShore University HealthSystem" , "NSU", "JHMI"))))

# Annotate metastatic samples (M_stage = 1)
filt_pdata$GSE21501$M_stage = ifelse(filt_pdata$GSE21501$Biomaterial_provider == "NEB", "1", "0") 

# Keep only JHMI samples (no neoadjuvant chemotherapy and clinical data availability)
filt_pdata$GSE21501 = filt_pdata$GSE21501 %>% 
  filter(Biomaterial_provider == "JHMI")

# Add AJCC classification according to TNM classification
filt_pdata$GSE21501 = filt_pdata$GSE21501 %>%
  mutate(AJCC_classification = c("2b","2b","2b","2b","2b","2b","2b","2b","2b","2b","2b","2b","2b","2b","2b","3","2b","2b","1b","2a","2b","2b","2b","2b","2b","2b","2b","2b","2b","2b","2b","2b","2b","2b"))

# Reorder columns
filt_pdata$GSE21501 = filt_pdata$GSE21501 %>%
  select(GEO_accession, Patient_ID, Platform, Biomaterial_provider, Tissue_type, AJCC_classification, T_stage, N_stage, M_stage, Survival_status, Survival_months)

# Transform pData to universal levels
for (i in 1:length(filt_pdata$GSE21501$Tissue_type)) {
  if(filt_pdata$GSE21501$Tissue_type[i] == "Patient Primary Pancreatic Tumor") {
    filt_pdata$GSE21501$Tissue_type[i] = "tumor"
  } else {
    filt_pdata$GSE21501$Tissue_type[i] = "non_tumor"
  }
}; rm(i)

filt_pdata$GSE21501$Tissue_type = factor(x = filt_pdata$GSE21501$Tissue_type, levels = c("non_tumor", "tumor"), labels = c("non_tumor", "tumor"))
filt_pdata$GSE21501$AJCC_classification = factor(x = filt_pdata$GSE21501$AJCC_classification, levels = c("1a","1b","2a","2b","3","4"), labels = c("1a","1b","2a","2b","3","4"))
filt_pdata$GSE21501$T_stage = factor(x = filt_pdata$GSE21501$T_stage, levels = c("1","2","3","4"), labels = c("1","2","3","4"))
filt_pdata$GSE21501$N_stage = factor(x = filt_pdata$GSE21501$N_stage, levels = c("0","1","2"), labels = c("0","1","2"))
filt_pdata$GSE21501$M_stage = factor(x = filt_pdata$GSE21501$M_stage, levels = c("0", "1"), labels = c("non_metastatic", "metastatic"))
filt_pdata$GSE21501$Survival_status = factor(x = filt_pdata$GSE21501$Survival_status, levels = c("1", "0"), labels = c("alive", "dead"))
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
#   - Metastatic samples (LM and PM) were contaminated with respectively normal liver and peritoneal tissue, reflecting in upregulation 
#     of liver and peritoneal specific genes. Therefore only genes that were not differentially expressed between LM and PM samples, 
#     considered as metastatic specific genes, were used for analysis between primary tumor and metastatic tissue.

# Select pData
filt_pdata[["GSE42952"]] = pdata$GSE42952 %>%
  select(GEO_accession = geo_accession,
         Patient_ID = title,
         Platform = platform_id,
         Age = "age:ch1",
         Gender = "gender:ch1",
         Prognosis = "prognosis:ch1",
         AJCC_classification = "Stage:ch1",
         Tissue_type = "tissue:ch1",
         TNM_classification = "tumor stage:ch1")

# Samples: "PDAC_P_117", "PDAC_P_104", "PDAC_T_105", "PDAC_P_105", "PDAC_T_91", "PDAC_P_91", "PDAC_T_138", "PDAC_P_138", "PDAC_T_123", 
# "PDAC_P_123" are identical to samples from GSE18670. Thus, they will be removed from here (position of elements: 1-10).
filt_pdata$GSE42952 = filt_pdata$GSE42952[-c(1:10),]

# Remove samples from non PDAC tissues (e.g. "liver metastasis from pancreatic cancer", "peritoneal metastasis from pancreatic cancer")
filt_pdata$GSE42952 = filt_pdata$GSE42952 %>% filter(`Tissue_type` != "liver metastasis from pancreatic cancer" &  `Tissue_type` != "peritoneal metastasis from pancreatic cancer")

# Manually addition of extra data from Table 1 (https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC3511800/)
filt_pdata$GSE42952 = filt_pdata$GSE42952 %>% 
  mutate(Tumor_localization = c("head", "head", "head", "head", "head", "head", "head", "head", "head", "head", "head", "head"),
         Tumor_grade = c(NA,1,3,3,2,3,3,3,2,2,3,3),
         T_stage = c(NA,3,3,3,3,3,3,3,3,2,3,3),
         N_stage = c(NA,0,0,1,1,0,0,1,1,0,1,1),
         M_stage = c(1,0,0,0,0,0,0,0,0,0,0,0),
         Perineural_invasion = c(NA,0,1,1,1,1,1,0,1,1,1,1),
         Vascular_invasion = c(NA,0,0,1,1,1,1,1,1,1,0,0)) %>%
  select(-TNM_classification)

# Reorder columns
filt_pdata$GSE42952 = filt_pdata$GSE42952 %>%
  select(GEO_accession, Patient_ID, Platform, Gender, Age, Tissue_type, AJCC_classification, T_stage, N_stage, M_stage, Tumor_grade, Perineural_invasion, Vascular_invasion, Prognosis, Tumor_localization)

# Trasnform pData to universal levels
for (i in 1:length(filt_pdata$GSE42952$Tissue_type)) {
  if(filt_pdata$GSE42952$Tissue_type[i] == "pancreatic ductal adenocarcinoma (PDAC)") {
    filt_pdata$GSE42952$Tissue_type[i] = "tumor"
  } else {
    filt_pdata$GSE42952$Tissue_type[i] = "non_tumor"
  }
}; rm(i)

filt_pdata$GSE42952$Tissue_type = factor(x = filt_pdata$GSE42952$Tissue_type, levels = c("non_tumor", "tumor"), labels = c("non_tumor", "tumor"))
filt_pdata$GSE42952$AJCC_classification = factor(x = filt_pdata$GSE42952$AJCC_classification, levels = c("1a","1b","2a","2b","3","4"), labels = c("1a","1b","2a","2b","3","4"))
filt_pdata$GSE42952$T_stage = factor(x = filt_pdata$GSE42952$T_stage, levels = c("1","2","3","4"), labels = c("1","2","3","4"))
filt_pdata$GSE42952$N_stage = factor(x = filt_pdata$GSE42952$N_stage, levels = c("0","1","2"), labels = c("0","1","2"))
filt_pdata$GSE42952$M_stage = factor(x = filt_pdata$GSE42952$M_stage, levels = c("0", "1"), labels = c("non_metastatic", "metastatic"))
filt_pdata$GSE42952$Tumor_grade = factor(x = filt_pdata$GSE42952$Tumor_grade, levels = c("1","2","3"), labels = c("1","2","3"))
filt_pdata$GSE42952$Perineural_invasion = factor(x = filt_pdata$GSE42952$Perineural_invasion, levels = c(0, 1), labels = c("no_PNI", "PNI"))
filt_pdata$GSE42952$Vascular_invasion = factor(x = filt_pdata$GSE42952$Vascular_invasion, levels = c(0,1), labels = c("no_VI", "VI"))
filt_pdata$GSE42952$Gender = factor(x = filt_pdata$GSE42952$Gender, levels = c("female", "male"), labels = c("female", "male"))
filt_pdata$GSE42952$Age = gsub("years", "", filt_pdata$GSE42952$Age)
filt_pdata$GSE42952$Age = as.numeric(filt_pdata$GSE42952$Age)
filt_pdata$GSE42952$Prognosis = factor(x = filt_pdata$GSE42952$Prognosis, levels = c("good", "bad"), labels = c("good", "bad"))
filt_pdata$GSE42952$Tumor_localization = factor(x = filt_pdata$GSE42952$Tumor_localization, levels = c("head", "body/tail"), labels = c("head", "body/tail"))

# GSE18670
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18670
# Paper: https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC3599097/
# Comments:
#   - 24 samples (6 patients)
#       - 6 circulating tumor (CTC)
#       - 6 haematological (G)
#       - 6 original tumor (T)
#       - 6 non-tumor pancreatic control (P)

# Select pData
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
  select(GEO_accession, Patient_ID, Platform, Gender, Age, Tissue_type, AJCC_classification, T_stage, N_stage, M_stage, Tumor_grade)

# Transform pData to universal levels
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

# GSE62452
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62452
# Paper: https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC4930741/
# Comments:
#   - 130 samples
#     - 69 pancreatic tumor
#     - 61 adjacent pancreatic non-tumor
#   - For tumor grading, the general classification was used (GX,G1,G2,G3,G4). In all the other samples, the PDAC classification
#     was used (G1,G2,G3,G4)

# Select pData
filt_pdata[["GSE62452"]] = pdata$GSE62452 %>%
  select(GEO_accession = geo_accession,
         Patient_ID = title,
         Platform = platform_id,
         Tissue_type = "tissue:ch1",
         AJCC_classification = "Stage:ch1",
         Tumor_grade = "grading:ch1")
# Didn't keep information regarding survival status and survival time because classification was not clear.

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

# Transform pData to universall levels
filt_pdata$GSE62452$Tissue_type = factor(x = filt_pdata$GSE62452$Tissue_type, levels = c("non_tumor","tumor"), labels = c("non_tumor","tumor"))
filt_pdata$GSE62452$AJCC_classification = gsub(">IIB", "3", fixed = TRUE, filt_pdata$GSE62452$AJCC_classification)
filt_pdata$GSE62452$AJCC_classification = gsub("IIB", "2b", fixed = TRUE, filt_pdata$GSE62452$AJCC_classification)
filt_pdata$GSE62452$AJCC_classification = gsub("IB", "1b", fixed = TRUE, filt_pdata$GSE62452$AJCC_classification)
filt_pdata$GSE62452$AJCC_classification = gsub("IIA", "2a", fixed = TRUE, filt_pdata$GSE62452$AJCC_classification)
filt_pdata$GSE62452$AJCC_classification = gsub("III", "3", fixed = TRUE, filt_pdata$GSE62452$AJCC_classification)
filt_pdata$GSE62452$AJCC_classification = gsub("IVA", "4", fixed = TRUE, filt_pdata$GSE62452$AJCC_classification)
filt_pdata$GSE62452$AJCC_classification = gsub("IVB", "4", fixed = TRUE, filt_pdata$GSE62452$AJCC_classification)
filt_pdata$GSE62452$AJCC_classification = gsub("IV", "4", fixed = TRUE, filt_pdata$GSE62452$AJCC_classification)
filt_pdata$GSE62452$AJCC_classification = factor(x = filt_pdata$GSE62452$AJCC_classification, levels = c("1a","1b","2a","2b","3","4"), labels = c("1a","1b","2a","2b","3","4"))
filt_pdata$GSE62452$Tumor_grade = factor(x = filt_pdata$GSE62452$Tumor_grade, levels = c("x","1","2","3", "4"), labels = c("GX","G1","G2","G3","G4"))

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
  select(GEO_accession = geo_accession,
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

# Transform pdata to universal levels
filt_pdata$GSE62165$Tissue_type = factor(x = filt_pdata$GSE62165$Tissue_type, levels = c("non_tumor","tumor"), labels = c("non_tumor","tumor"))
filt_pdata$GSE62165$AJCC_classification = factor(x = filt_pdata$GSE62165$AJCC_classification, levels = c("1a","1b","2a","2b","3","4"), labels = c("1a","1b","2a","2b","3","4"))

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
  select(GEO_accession = geo_accession,
         Patient_ID = title,
         Platform = platform_id,
         Gender = "genderοΌmale=0οΌfemale=1οΌ‰:ch1",                               # 0: male           1: female
         Age_grouped = "ageοΌβ‰¤65=0οΌοΌ65=1οΌ‰:ch1",                             # 0: <=65           1: >65
         Tissue_type =  "source_name_ch1",
         T_stage = "t stage:ch1",
         N_stage = "n stageοΌβ€n0=0,n1=1):ch1",
         M_stage = "m stage(m0=0.m1=1):ch1",
         Vascular_invasion = "vessel invasion(absence=0,present=1):ch1",             # 0: absence        1: present
         Differentiation = "differentiation(well/moderate=0,poor=1):ch1",            # 0: well/moderate  1: poor
         Survival_status = "survival status(alive=0,dead=1):ch1",                    # 0: alive          1: dead
         Survival_days = "survival time(day):ch1",
         Tumor_localization = "localization of tumor(head=1,body/tail=2):ch1")       # 1: head           2: body/tail

# Separate information in Tissue_type
# Patient_ID
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

# Transform pdata to universal levels
filt_pdata$GSE102238$Gender = factor(x = filt_pdata$GSE102238$Gender, levels = c("1","0"), labels = c("female","male"))
filt_pdata$GSE102238$Age_grouped = factor(x = filt_pdata$GSE102238$Age_grouped, levels = c("0","1"), labels = c("<=65", ">65"))
filt_pdata$GSE102238$Tissue_type = factor(x = filt_pdata$GSE102238$Tissue_type, levels = c("non_tumor","tumor"), labels = c("non_tumor","tumor"))
filt_pdata$GSE102238$T_stage = factor(x = filt_pdata$GSE102238$T_stage, levels = c("1","2","3","4"), labels = c("1","2","3","4"))
filt_pdata$GSE102238$N_stage = factor(x = filt_pdata$GSE102238$N_stage, levels = c("0","1","2"), labels = c("0","1","2"))
filt_pdata$GSE102238$M_stage = factor(x = filt_pdata$GSE102238$M_stage, levels = c("0", "1"), labels = c("non_metastatic", "metastatic"))
filt_pdata$GSE102238$Perineural_invasion = factor(x = filt_pdata$GSE102238$Perineural_invasion, levels = c("no_PNI","PNI"), labels = c("no_PNI","PNI"))
filt_pdata$GSE102238$Vascular_invasion = factor(x = filt_pdata$GSE102238$Vascular_invasion, levels = c("0","1"), labels = c("no_VI","VI"))
filt_pdata$GSE102238$Differentiation = factor(x = filt_pdata$GSE102238$Differentiation, levels = c("0","1"), labels = c("well/moderate","poor"))
filt_pdata$GSE102238$Survival_status = factor(x = filt_pdata$GSE102238$Survival_status, levels = c("0","1"), labels = c("alive","dead"))
filt_pdata$GSE102238$Tumor_localization = factor(x = filt_pdata$GSE102238$Tumor_localization, levels = c("1","2"), labels = c("head", "body/tail"))
# Transform survival_days to survival_months
filt_pdata$GSE102238$Survival_days = as.numeric(filt_pdata$GSE102238$Survival_days)
filt_pdata$GSE102238 = filt_pdata$GSE102238 %>% mutate(Survival_months = Survival_days/30)

# Add AJCC classification according to TNM classification
for (i in 1:length(filt_pdata$GSE102238$T_stage)) {
  if(filt_pdata$GSE102238$T_stage[i] == 1 && filt_pdata$GSE102238$N_stage[i] == 0 && filt_pdata$GSE102238$M_stage[i] == "non_metastatic") {
    filt_pdata$GSE102238$AJCC_classification[i] = "1a"
  }
  if(filt_pdata$GSE102238$T_stage[i] == 2 && filt_pdata$GSE102238$N_stage[i] == 0 && filt_pdata$GSE102238$M_stage[i] == "non_metastatic") {
    filt_pdata$GSE102238$AJCC_classification[i] = "1b"
  }
  if(filt_pdata$GSE102238$T_stage[i] == 3 && filt_pdata$GSE102238$N_stage[i] == 0 && filt_pdata$GSE102238$M_stage[i] == "non_metastatic") {
    filt_pdata$GSE102238$AJCC_classification[i] = "2a"
  }
  if(filt_pdata$GSE102238$T_stage[i] %in% 1:3 && filt_pdata$GSE102238$N_stage[i] == 1 && filt_pdata$GSE102238$M_stage[i] == "non_metastatic") {
    filt_pdata$GSE102238$AJCC_classification[i] = "2b"
  }
  if(filt_pdata$GSE102238$T_stage[i] == 4 && filt_pdata$GSE102238$N_stage[i] %in% 0:2 && filt_pdata$GSE102238$M_stage[i] == "non_metastatic") {
    filt_pdata$GSE102238$AJCC_classification[i] = "3"
  }
  if(filt_pdata$GSE102238$T_stage[i] %in% 1:4 && filt_pdata$GSE102238$N_stage[i] == 2 && filt_pdata$GSE102238$M_stage[i] == "non_metastatic") {
    filt_pdata$GSE102238$AJCC_classification[i] = "3"
  }
  if(filt_pdata$GSE102238$T_stage[i] %in% 1:4 && filt_pdata$GSE102238$N_stage[i] %in% 0:2 && filt_pdata$GSE102238$M_stage[i] == "metastatic") {
    filt_pdata$GSE102238$AJCC_classification[i] = "4"
  }
}; rm(i)

filt_pdata$GSE102238$AJCC_classification = factor(x = filt_pdata$GSE102238$AJCC_classification, levels = c("1a","1b","2a","2b","3","4"), labels = c("1a","1b","2a","2b","3","4"))

# Reorder columns
filt_pdata$GSE102238 = filt_pdata$GSE102238 %>%
  select(GEO_accession:Tissue_type, AJCC_classification, T_stage:M_stage, Perineural_invasion, Vascular_invasion, Differentiation, Survival_status, Survival_days, Survival_months,Tumor_localization)

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
  select(GEO_accession = geo_accession,
         Patient_ID = description,
         Platform = platform_id,
         AJCC_classification = "tumor_stage:ch1",
         Survival_months = "months_survival:ch1")

# Add Tisse_type column
filt_pdata$GSE84219$Tissue_type = rep("tumor", length(filt_pdata$GSE84219$GEO_accession))

# Transform pdata to universall levels
filt_pdata$GSE84219$Tissue_type = factor(x = filt_pdata$GSE84219$Tissue_type, levels = c("non_tumor","tumor"), labels = c("non_tumor","tumor"))
filt_pdata$GSE84219$AJCC_classification = gsub("IIB", "2b", fixed = TRUE, filt_pdata$GSE84219$AJCC_classification)
filt_pdata$GSE84219$AJCC_classification = gsub("IB", "1b", fixed = TRUE, filt_pdata$GSE84219$AJCC_classification)
filt_pdata$GSE84219$AJCC_classification = gsub("IIA", "2a", fixed = TRUE, filt_pdata$GSE84219$AJCC_classification)
filt_pdata$GSE84219$AJCC_classification = gsub("IA", "1a", fixed = TRUE, filt_pdata$GSE84219$AJCC_classification)
filt_pdata$GSE84219$AJCC_classification = factor(x = filt_pdata$GSE84219$AJCC_classification, levels = c("1a","1b","2a","2b","3","4"), labels = c("1a","1b","2a","2b","3","4"))

# Reorder columns
filt_pdata$GSE84219 = filt_pdata$GSE84219 %>%
  select(GEO_accession, Patient_ID, Platform, Tissue_type, AJCC_classification, Survival_months)

##### Expression data #####
# Create the esets objects
esets = list()
for (i in 1:length(GEOsets)) {
  esets[[i]] = exprs(GEOsets[[i]])
  esets[[i]] = as.data.frame(esets[[i]])
}
names(esets) = datasets; rm(i)

## Remove unneeded/inappropriate samples
# GSE21501
# Remove sample GSM506033 because it was used as reference (excluded from pdata, too)
esets[["GSE21501"]] = esets[["GSE21501"]][, -6]
# Keep only expression data for JHMI samples (no neoadjuvant chemotherapy)
esets$GSE21501 = esets$GSE21501[, colnames(esets$GSE21501) %in% filt_pdata[["GSE21501"]][["GEO_accession"]]]

# GSE42952
# Remove first 10 samples because they are identical to the ones of GSE18670
# Remove non PDAC tissue samples
esets[["GSE42952"]] = esets[["GSE42952"]][, -c(1:10, 13,14,17,18,21,22,25,26,28,30,33)]

# GSE18670 
# Keep only samples for tumor (T) and non_tumor (P)
esets[["GSE18670"]] = esets$GSE18670[, c("GSM463723", "GSM463724", "GSM463727", "GSM463728", "GSM463731", "GSM463732", "GSM463735", "GSM463736", "GSM463739", "GSM463740", "GSM463743", "GSM463744")]

##### Annotation esets with Entrez ID's #####
# GSE21501
fdata21501 = fData(GEOsets$GSE21501) %>%
  dplyr::select(ID, ENTREZ_GENE_ID = GENE) %>%
  dplyr::filter(!grepl("///", ENTREZ_GENE_ID)) %>%
  dplyr::filter(nchar(ENTREZ_GENE_ID)>0)
esets[[1]]$ID = as.numeric(rownames(esets[[1]]))
esets[[1]] = esets[[1]] %>%
  inner_join(fdata21501) %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  summarise_all(mean, na.rm = TRUE)
esets[[1]]$ENTREZ_GENE_ID = as.numeric(esets[[1]]$ENTREZ_GENE_ID)
rm(fdata21501)

# GSE42952
fdata42952 = fData(GEOsets$GSE42952) %>%
  dplyr::select(ID, ENTREZ_GENE_ID) %>%
  dplyr::filter(!grepl("///", ENTREZ_GENE_ID)) %>%
  dplyr::filter(nchar(ENTREZ_GENE_ID)>0)
esets[[2]]$ID = rownames(esets[[2]])
esets[[2]] = esets[[2]] %>%
  inner_join(fdata42952) %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  summarise_all(mean, na.rm = TRUE)
esets[[2]]$ENTREZ_GENE_ID = as.numeric(esets[[2]]$ENTREZ_GENE_ID)
rm(fdata42952)

# GSE18670
fdata18670 = fData(GEOsets$GSE18670) %>%
  dplyr::select(ID, ENTREZ_GENE_ID) %>%
  dplyr::filter(!grepl("///", ENTREZ_GENE_ID)) %>%
  dplyr::filter(nchar(ENTREZ_GENE_ID)>0)
esets[[3]]$ID = rownames(esets[[3]])
esets[[3]] = esets[[3]] %>%
  inner_join(fdata18670) %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  summarise_all(mean, na.rm = TRUE)
esets[[3]]$ENTREZ_GENE_ID = as.numeric(esets[[3]]$ENTREZ_GENE_ID)
rm(fdata18670)

# GSE62452
# Convert GB_LIST (RefSeq) to ENTREZ IDs
ref = org.Hs.egREFSEQ2EG
mapped_genes_official = mappedkeys(ref)
ref_df = as.data.frame(ref[mapped_genes_official])
ref_df = ref_df %>% dplyr::rename(EntrezGene.ID = gene_id, RefSeq = accession)
length(unique(ref_df$RefSeq)) == length(ref_df$RefSeq) # TRUE --> No duplicates in RefSeq

# From GB_LIST remove:
#   - Values that have comma (,) (multiple genes to a probe)
#   - Everything after the bullet (.)
fdata62452 = fData(GEOsets$GSE62452) %>%
  dplyr::select(ID, RefSeq = GB_LIST) %>%
  dplyr::filter(!grepl(",", RefSeq)) %>%
  dplyr::filter(nchar(RefSeq)>0) %>%
  mutate(RefSeq = gsub("\\..*","", RefSeq)) %>%
  inner_join(ref_df, by = "RefSeq") %>%
  dplyr::select(ID, ENTREZ_GENE_ID = EntrezGene.ID)
esets[[4]]$ID = as.numeric(rownames(esets[[4]]))
esets[[4]] = esets[[4]] %>%
  inner_join(fdata62452) %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  summarise_all(mean, na.rm = TRUE)
esets[[4]]$ENTREZ_GENE_ID = as.numeric(esets[[4]]$ENTREZ_GENE_ID)
rm(fdata62452)

# GSE62165
fdata62165 = fData(GEOsets$GSE62165) %>%
  dplyr::select(ID, ENTREZ_GENE_ID = 'Entrez Gene') %>%
  dplyr::filter(!grepl("///", ENTREZ_GENE_ID)) %>%
  dplyr::filter(nchar(ENTREZ_GENE_ID)>0)
esets[[5]]$ID = rownames(esets[[5]])
esets[[5]] = esets[[5]] %>%
  inner_join(fdata62165) %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  dplyr::filter(ENTREZ_GENE_ID != "---")
esets[[5]]$ENTREZ_GENE_ID = as.numeric(esets[[5]]$ENTREZ_GENE_ID)
rm(fdata62165)

# GSE102238
# No information for corresponding genes provided in GEOset
# Manually perform BLASTN alignment to annotate probes with RefSeq gene IDs
# https://blast.ncbi.nlm.nih.gov/Blast.cgi
# Read csv results from BLASTN
blastn_1 = read.csv(file = "GSE102238_BLASTN-1.csv", header = FALSE)
blastn_2 = read.csv(file = "GSE102238_BLASTN-2.csv", header = FALSE)
blastn = rbind(blastn_1, blastn_2); rm(blastn_1, blastn_2)

# Remove duplicate probes (V1 columns)
dups = blastn$V1[which(duplicated(blastn$V1))]
final_blastn = blastn[!blastn$V1 %in% dups, ] %>%
  select(ID = V1, RefSeq = V2)

# Convert RefSeq IDs to ENTREZ IDs
ref = org.Hs.egREFSEQ2EG
mapped_genes_official = mappedkeys(ref)
ref_df = as.data.frame(ref[mapped_genes_official])
ref_df = ref_df %>% dplyr::rename(EntrezGene.ID = gene_id, RefSeq = accession)
length(unique(ref_df$RefSeq)) == length(ref_df$RefSeq) # TRUE --> No duplicates in RefSeq

annot = inner_join(final_blastn, ref_df, by = "RefSeq")

# GSE84219
fdata84219 = fData(GEOsets$GSE84219) %>%
  dplyr::select(ID, ENTREZ_GENE_ID = Entrez_Gene_ID) %>%
  dplyr::filter(!grepl("///", ENTREZ_GENE_ID)) %>%
  dplyr::filter(nchar(ENTREZ_GENE_ID)>0)
esets[[7]]$ID = rownames(esets[[7]])
esets[[7]] = esets[[7]] %>%
  inner_join(fdata84219) %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  group_by(ENTREZ_GENE_ID) %>%
  summarise_all(mean, na.rm = TRUE)
esets[[7]]$ENTREZ_GENE_ID = as.numeric(esets[[7]]$ENTREZ_GENE_ID)
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






