# Circos preparations
# Matching our genes with DrugBank and preparing Circos input
# Code for handling DrugBank data and matching it to our genes

# Libraries #####
library(dplyr)
library(readr)
library(reshape2)
library(org.Hs.eg.db)
library(openxlsx)
library(stringr)
library(tidyr)

# Importing .RData file produced after running the "Pathway_ontology_driver_genes.R" script
# load("C:/your/path/to/the/environment.RData") 

# Keep only the clustered results plus the miRNA analysis results:
rm(list=setdiff(ls(), c("clustered_results_stage_1", "clustered_results_stage_2",
                        "clustered_results_stage_3", "clustered_results_stage_4",
                        "clustered_results_blood", "clean_miRNA_venn", "miRNA_venn", "sig_DGEA")))

##### Importing ATC drug category data #####
ATC = read.xlsx("ATC_Drugs.xlsx", sheet = 3) # full set
ATC$ATC_Name = gsub(" ", "_", ATC$ATC_Name)
ATC = ATC %>% dplyr::rename(API = ATC_Name)
antineoplastics = ATC %>%
  dplyr::filter(Name1 == "ANTINEOPLASTIC AGENTS")
non_antineoplastics = anti_join(ATC, antineoplastics)

##### Importing DrugBank data #####
# We downloaded the approved protein identifier data from the website
# along with the open data section directory. Here we import the directory as
# "directory_pre" because it needs processing and the interaction data is
# similarly imported as "interaction_pre". Files where drug identifiers were in
# one column but separated by semicolons, were subjected to "text-to-column"
# transformation in MS Excel firstly.

directory_pre = read_csv("DrugBank_Directory_pre.csv")
interaction_pre_pharmacological = read_csv("DrugBank_Interactions_pharmacological_pre.csv")
interaction_pre = read_csv("DrugBank_Interactions_all_pre.csv")

interaction_pre_pharmacological = interaction_pre_pharmacological %>%
  dplyr::filter(Species=="Humans") %>%
  dplyr::select(-ID, -Species, -`GeneCard ID`, -`PDB ID`, -`Uniprot Title`,
                -`GenBank Gene ID`, -`GenBank Protein ID`, -`HGNC ID`, -`GenAtlas ID`)
interaction_pre = interaction_pre %>%
  dplyr::filter(Species=="Humans") %>%
  dplyr::select(-ID, -Species, -`GeneCard ID`, -`PDB ID`, -`Uniprot Title`,
                -`GenBank Gene ID`, -`GenBank Protein ID`, -`HGNC ID`, -`GenAtlas ID`)

# Melting #####
molten_pharmacological = melt(interaction_pre_pharmacological, 
                              id.vars = colnames(interaction_pre_pharmacological)[1:3])
molten_pharmacological = molten_pharmacological %>% 
  dplyr::select(-variable) %>%
  dplyr::rename(DrugBank_ID = value) %>%
  dplyr::filter(is.na(DrugBank_ID)==F)


molten_all = melt(interaction_pre, 
                  id.vars = colnames(interaction_pre)[1:3])
molten_all = molten_all %>% 
  dplyr::select(-variable) %>%
  dplyr::rename(DrugBank_ID = value) %>%
  dplyr::filter(is.na(DrugBank_ID)==F)

##### Gene annotation #####
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

ID_Map$EntrezGene.ID = as.character(ID_Map$EntrezGene.ID)
ID_Map = ID_Map %>% dplyr::rename(Gene.Symbol = probe)
rm(alias, alias_df, aliases_for_join, official,
   mapped_genes_alias, mapped_genes_official)

# Account for non-official gene names in the DrugBank sets
dups = setdiff(unique(molten_all$`Gene Name`),
               official_df$Gene.Symbol)
dups = dups[!dups %in% c("VMAT2", "PDE4", "ASS", "NOV", "GIG18", "GAD65", "TLH6",
                         "EEF1A1L14", "HIST1H2BC", "EFTUD1", "DKFZp686P18130")]

# ATP5L has been discontinued in January 2022

dups = na.omit(dups)
dups_actual = c("EPRS1", "ND3", "KARS1", "VARS1", "WARS1", "YARS1", "TARS1", "ND4L",
                "KYAT1", "COX1", "ND6", "ATP5F1D", "ND1", "ND4", "ND5", "NARS1",
                "ND2", "GFUS", "CARS1", "DARS1", "MMUT", "ALPG", "GARS1", "ACP3",
                "PLPBP", "KYAT3", "ADSS2", "COX2", "ATP5F1A", "COX3", "H1-4", "SELENOP",
                "MAP3K20", "OGA", "ATP6", "ATP5MG", "ATP5PO", "CILK1", "ADSS1")
replacements = as.data.frame(dups)
replacements$actual = dups_actual

molten_all = molten_all %>% dplyr::filter(is.na(`Gene Name`)==F)
molten_pharmacological = molten_pharmacological %>% dplyr::filter(is.na(`Gene Name`)==F)
molten_all$newname = molten_all$`Gene Name`
molten_pharmacological$newname = molten_pharmacological$`Gene Name`

for (i in 1:nrow(replacements)){
  for (j in 1:nrow(molten_all)){
    if (molten_all$`Gene Name`[j] == replacements$dups[i]){
      molten_all$newname[j] = replacements$actual[i]
    }
  }
}

for (i in 1:nrow(replacements)){
  for (j in 1:nrow(molten_pharmacological)){
    if (molten_pharmacological$`Gene Name`[j] == replacements$dups[i]){
      molten_pharmacological$newname[j] = replacements$actual[i]
    }
  }
}

molten_all_annot = molten_all %>% dplyr::rename(Gene.Symbol = newname) %>%
  dplyr::select(-`Gene Name`) %>% 
  inner_join(official_df, by = "Gene.Symbol")
molten_pharmacological_annot = molten_pharmacological %>% dplyr::rename(Gene.Symbol = newname) %>%
  dplyr::select(-`Gene Name`) %>% 
  inner_join(official_df, by = "Gene.Symbol")

molten_all_annot = molten_all_annot %>% dplyr::select(-`UniProt ID`)
molten_pharmacological_annot = molten_pharmacological_annot %>% 
  dplyr::select(-`UniProt ID`)

##### Refining our drug directory #####
directory = directory_pre %>% dplyr::rename(DrugBank_ID = `DrugBank ID`,
                                            API = `Common name`) %>%
  dplyr::select(DrugBank_ID, API) %>% distinct()

final_all = molten_all_annot %>% inner_join(directory, by = "DrugBank_ID") %>%
  dplyr::select(-DrugBank_ID)
final_pharmacological = molten_pharmacological_annot %>% 
  inner_join(directory, by = "DrugBank_ID") %>%
  dplyr::select(-DrugBank_ID)

##### Matching DrugBank with ATC prior to joining with DE genes for drug ontology #####
final_all$API = tolower(final_all$API)
final_pharmacological$API = tolower(final_pharmacological$API)
# Replace spaces with underscores here too:
final_all$API = gsub(" ", "_", final_all$API)
final_pharmacological$API = gsub(" ", "_", final_pharmacological$API)

ultimate_pharmacological_pre = final_pharmacological %>%
  left_join(ATC, by = "API")
ultimate_all_pre = final_all %>%
  left_join(ATC, by = "API")
unmatched_pharmacological = unique(ultimate_pharmacological_pre$API[
  which(is.na(ultimate_pharmacological_pre$Sub4))])
unmatched_all = unique(ultimate_all_pre$API[
  which(is.na(ultimate_all_pre$Sub4))])

ultimate_pharmacological_pre$NewName = NA
for (i in 1:nrow(ultimate_pharmacological_pre)){
  if(ultimate_pharmacological_pre$API[i] %in%
     unmatched_pharmacological){
    ultimate_pharmacological_pre$NewName[i] = NA
  } else {
    ultimate_pharmacological_pre$NewName[i] = ultimate_pharmacological_pre$API[i]
  }
}

ultimate_all_pre$NewName = NA
for (i in 1:nrow(ultimate_all_pre)){
  if(ultimate_all_pre$API[i] %in%
     unmatched_all){
    ultimate_all_pre$NewName[i] = NA
  } else {
    ultimate_all_pre$NewName[i] = ultimate_all_pre$API[i]
  }
}

# Here, manual matching must be performed for the unmatched drug names
# Hide the section if it takes too much space

##### Manual matching, just for drugs (not amino-acids, chemicals, radio etc.) #####
ultimate_all_pre$NewName[ultimate_all_pre$API=="nadh"] = "nicotinamide"
ultimate_all_pre$NewName[ultimate_all_pre$API=="drotrecogin_alfa"] = "drotrecogin_alfa_(activated)"
ultimate_all_pre$NewName[ultimate_all_pre$API=="vitamin_a"] = "retinol_(vit_A)"
ultimate_all_pre$NewName[ultimate_all_pre$API=="insulin_human"] = "insulin_(human)"
ultimate_all_pre$NewName[ultimate_all_pre$API=="copper"] = "copper_sulfate" # most appropriate
ultimate_all_pre$NewName[ultimate_all_pre$API=="pyridoxine"] = "pyridoxine_(vit_B6)"
ultimate_all_pre$NewName[ultimate_all_pre$API=="somatotropin"] = "somatropin"
ultimate_all_pre$NewName[ultimate_all_pre$API=="vitamin_e"] = "tocopherol_(vit_E)"
ultimate_all_pre$NewName[ultimate_all_pre$API=="telotristat_ethyl"] = "telotristat"
ultimate_all_pre$NewName[ultimate_all_pre$API=="amphetamine"] = "amfetamine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="glyburide"] = "glibenclamide"
ultimate_all_pre$NewName[ultimate_all_pre$API=="thiamine"] = "thiamine_(vit_B1)"
ultimate_all_pre$NewName[ultimate_all_pre$API=="leuprolide"] = "leuprorelin"
ultimate_all_pre$NewName[ultimate_all_pre$API=="zinc"] = "zinc_sulfate" # same as copper
ultimate_all_pre$NewName[ultimate_all_pre$API=="adenosine_phosphate"] = "adenosine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="muromonab"] = "muromonab-CD3"
ultimate_all_pre$NewName[ultimate_all_pre$API=="ardeparin"] = "enoxaparin"
ultimate_all_pre$NewName[ultimate_all_pre$API=="liothyronine"] = "liothyronine_sodium"
ultimate_all_pre$NewName[ultimate_all_pre$API=="secretin_human"] = "secretin"
ultimate_all_pre$NewName[ultimate_all_pre$API=="pentosan_polysulfate"] = "pentosan_polysulfate_sodium"
ultimate_all_pre$NewName[ultimate_all_pre$API=="coagulation_factor_viia_recombinant_human"] = "coagulation_factor_VIIa"
ultimate_all_pre$NewName[ultimate_all_pre$API=="cyclosporine"] = "ciclosporin"
ultimate_all_pre$NewName[ultimate_all_pre$API=="botulinum_toxin_type_a"] = "botulinum_toxin"
ultimate_all_pre$NewName[ultimate_all_pre$API=="salmon_calcitonin"] = "calcitonin_(salmon_synthetic)"
ultimate_all_pre$NewName[ultimate_all_pre$API=="interferon_gamma-1b"] = "interferon_gamma"
ultimate_all_pre$NewName[ultimate_all_pre$API=="niacin"] = "nicotinic_acid"
ultimate_all_pre$NewName[ultimate_all_pre$API=="capromab_pendetide"] = "indium_(111In)_capromab_pendetide"
ultimate_all_pre$NewName[ultimate_all_pre$API=="fluticasone_propionate"] = "fluticasone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="insulin_pork"] = "insulin_(pork)"
ultimate_all_pre$NewName[ultimate_all_pre$API=="lithium_citrate"] = "lithium" # only two lithium salts in the DrugBank files
ultimate_all_pre$NewName[ultimate_all_pre$API=="coagulation_factor_ix_(recombinant)"] = "coagulation_factor_IX"
ultimate_all_pre$NewName[ultimate_all_pre$API=="medroxyprogesterone_acetate"] = "medroxyprogesterone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="sipuleucel-t"] = "sipuleucel-T"
ultimate_all_pre$NewName[ultimate_all_pre$API=="ferrous_sulfate_anhydrous"] = "ferrous_sulfate"
ultimate_all_pre$NewName[ultimate_all_pre$API=="vitamin_d"] = "colecalciferol"
ultimate_all_pre$NewName[ultimate_all_pre$API=="indomethacin"] = "indometacin"
ultimate_all_pre$NewName[ultimate_all_pre$API=="botulinum_toxin_type_b"] = "botulinum_toxin"
ultimate_all_pre$NewName[ultimate_all_pre$API=="omega-3-carboxylic_acids"] = "omega-3-triglycerides_incl._other_esters_and_acids"
ultimate_all_pre$NewName[ultimate_all_pre$API=="dyphylline"] = "diprophylline"
ultimate_all_pre$NewName[ultimate_all_pre$API=="ezogabine"] = "retigabine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="glatiramer"] = "glatiramer_acetate"
ultimate_all_pre$NewName[ultimate_all_pre$API=="cysteamine"] = "mercaptamine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="acetaminophen"] = "paracetamol"
ultimate_all_pre$NewName[ultimate_all_pre$API=="pentamidine"] = "pentamidine_isethionate"
ultimate_all_pre$NewName[ultimate_all_pre$API=="rutin"] = "rutoside"
ultimate_all_pre$NewName[ultimate_all_pre$API=="anthralin"] = "dithranol"
ultimate_all_pre$NewName[ultimate_all_pre$API=="dalfampridine"] = "fampridine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="mycophenolate_mofetil"] = "mycophenolic_acid"
ultimate_all_pre$NewName[ultimate_all_pre$API=="fluciclovine_(18f)"] = "fluciclovine_(18F)"
ultimate_all_pre$NewName[ultimate_all_pre$API=="ouabain"] = "g-strophanthin"
ultimate_all_pre$NewName[ultimate_all_pre$API=="eslicarbazepine_acetate"] = "eslicarbazepine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="human_calcitonin"] = "calcitonin_(human_synthetic)"
ultimate_all_pre$NewName[ultimate_all_pre$API=="ergoloid_mesylate"] = "ergoloid_mesylates"
ultimate_all_pre$NewName[ultimate_all_pre$API=="butarbital"] = "butobarbital"
ultimate_all_pre$NewName[ultimate_all_pre$API=="estradiol_acetate"] = "estradiol"
ultimate_all_pre$NewName[ultimate_all_pre$API=="xanthinol"] = "xantinol_nicotinate"
ultimate_all_pre$NewName[ultimate_all_pre$API=="cholecystokinin"] = "pancreozymin_(cholecystokinin)"
ultimate_all_pre$NewName[ultimate_all_pre$API=="olmesartan"] = "olmesartan_medoxomil"
ultimate_all_pre$NewName[ultimate_all_pre$API=="carboprost_tromethamine"] = "carboprost"
ultimate_all_pre$NewName[ultimate_all_pre$API=="ibritumomab_tiuxetan"] = "ibritumomab_tiuxetan_(90Y)"
ultimate_all_pre$NewName[ultimate_all_pre$API=="dexamethasone_acetate"] = "dexamethasone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="sodium_ferric_gluconate_complex"] = "ferrous_gluconate"
ultimate_all_pre$NewName[ultimate_all_pre$API=="vasopressin"] = "vasopressin_(argipressin)"
ultimate_all_pre$NewName[ultimate_all_pre$API=="thimerosal"] = "thiomersal"
ultimate_all_pre$NewName[ultimate_all_pre$API=="chlorthalidone"] = "chlortalidone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="levothyroxine"] = "levothyroxine_sodium"
ultimate_all_pre$NewName[ultimate_all_pre$API=="meperidine"] = "pethidine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="isoflurophate"] = "fluostigmine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="nitroglycerin"] = "glyceryl_trinitrate"
ultimate_all_pre$NewName[ultimate_all_pre$API=="isoetharine"] = "isoetarine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="phylloquinone"] = "isoetarine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="lithium_carbonate"] = "lithium" # only two lithium salts in the DrugBank files
ultimate_all_pre$NewName[ultimate_all_pre$API=="coagulation_factor_vii"] = "coagulation_factor_VIIa"
ultimate_all_pre$NewName[ultimate_all_pre$API=="podofilox"] = "podophyllotoxin"
ultimate_all_pre$NewName[ultimate_all_pre$API=="human_thrombin"] = "thrombin"
ultimate_all_pre$NewName[ultimate_all_pre$API=="phylloquinone"] = "phytomenadione"
ultimate_all_pre$NewName[ultimate_all_pre$API=="gamma-aminobutyric_acid"] = "aminobutyric_acid"
ultimate_all_pre$NewName[ultimate_all_pre$API=="difluocortolone"] = "diflucortolone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="meclizine"] = "meclozine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="omega-3-_acid_ethyl_esters"] = "omega-3-triglycerides_incl._other_esters_and_acids"
ultimate_all_pre$NewName[ultimate_all_pre$API=="cyproterone_acetate"] = "cyproterone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="fenofibric_acid"] = "fenofibrate"
ultimate_all_pre$NewName[ultimate_all_pre$API=="estradiol_benzoate"] = "estradiol"
ultimate_all_pre$NewName[ultimate_all_pre$API=="hydrocortisone_acetate"] = "hydrocortisone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="tositumomab"] = "tositumomab/iodine_(131I)_tositumomab"
ultimate_all_pre$NewName[ultimate_all_pre$API=="orphenadrine"] = "orphenadrine_(citrate)"
ultimate_all_pre$NewName[ultimate_all_pre$API=="dinoprost_tromethamine"] = "dinoprost"
ultimate_all_pre$NewName[ultimate_all_pre$API=="megestrol_acetate"] = "megestrol"
ultimate_all_pre$NewName[ultimate_all_pre$API=="pentaerythritol_tetranitrate"] = "pentaerithrityl_tetranitrate"
ultimate_all_pre$NewName[ultimate_all_pre$API=="histamine"] = "histamine_phosphate"
ultimate_all_pre$NewName[ultimate_all_pre$API=="hyaluronidase_(human_recombinant)"] = "hyaluronidase"
ultimate_all_pre$NewName[ultimate_all_pre$API=="dextroamphetamine"] = "dexamfetamine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="methyldopa"] = "methyldopa_(racemic)"
ultimate_all_pre$NewName[ultimate_all_pre$API=="dipivefrin"] = "dipivefrine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="ibandronate"] = "ibandronic_acid"
ultimate_all_pre$NewName[ultimate_all_pre$API=="ethanolamine_oleate"] = "monoethanolamine_oleate"
ultimate_all_pre$NewName[ultimate_all_pre$API=="omega-3_fatty_acids"] = "omega-3-triglycerides_incl._other_esters_and_acids"
ultimate_all_pre$NewName[ultimate_all_pre$API=="estradiol_cypionate"] = "estradiol"
ultimate_all_pre$NewName[ultimate_all_pre$API=="protein_c"] = "protein_C"
ultimate_all_pre$NewName[ultimate_all_pre$API=="von_willebrand_factor_human"] = "von_Willebrand_factor"
ultimate_all_pre$NewName[ultimate_all_pre$API=="candesartan_cilexetil"] = "candesartan"
ultimate_all_pre$NewName[ultimate_all_pre$API=="liotrix"] = "combinations_of_levothyroxine_and_liothyronine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="hydroxyurea"] = "hydroxycarbamide"
ultimate_all_pre$NewName[ultimate_all_pre$API=="methimazole"] = "thiamazole"
ultimate_all_pre$NewName[ultimate_all_pre$API=="ipratropium"] = "ipratropium_bromide"
ultimate_all_pre$NewName[ultimate_all_pre$API=="metyrosine"] = "metirosine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="cholecalciferol"] = "colecalciferol"
ultimate_all_pre$NewName[ultimate_all_pre$API=="tegafur-uracil"] = "tegafur"
ultimate_all_pre$NewName[ultimate_all_pre$API=="coagulation_factor_ix_human"] = "coagulation_factor_IX"
ultimate_all_pre$NewName[ultimate_all_pre$API=="tromethamine"] = "trometamol"
ultimate_all_pre$NewName[ultimate_all_pre$API=="estradiol_dienanthate"] = "estradiol"
ultimate_all_pre$NewName[ultimate_all_pre$API=="hydrocortisone_cypionate"] = "hydrocortisone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="interferon_alfa-2a,_recombinant"] = "interferon_alfa-2a"
ultimate_all_pre$NewName[ultimate_all_pre$API=="clobetasol_propionate"] = "clobetasol"
ultimate_all_pre$NewName[ultimate_all_pre$API=="gallamine_triethiodide"] = "gallamine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="norelgestromin"] = "norelgestromin_and_ethinylestradiol"
ultimate_all_pre$NewName[ultimate_all_pre$API=="estradiol_valerate"] = "estradiol"
ultimate_all_pre$NewName[ultimate_all_pre$API=="hydrocortisone_phosphate"] = "hydrocortisone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="guaiacol"] = "guaiacolsulfonate"
ultimate_all_pre$NewName[ultimate_all_pre$API=="trastuzumab_deruxtecan"] = "trastuzumab"
ultimate_all_pre$NewName[ultimate_all_pre$API=="gabapentin_enacarbil"] = "gabapentin"
ultimate_all_pre$NewName[ultimate_all_pre$API=="nylidrin"] = "buphenine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="hydrocortisone_valerate"] = "hydrocortisone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="proparacaine"] = "proxymetacaine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="hydrocortisone_probutate"] = "hydrocortisone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="chorionic_gonadotropin_(human)"] = "choriogonadotropin_alfa"
ultimate_all_pre$NewName[ultimate_all_pre$API=="methsuximide"] = "mesuximide"
ultimate_all_pre$NewName[ultimate_all_pre$API=="lumacaftor"] = "ivacaftor_and_lumacaftor"
ultimate_all_pre$NewName[ultimate_all_pre$API=="florbetaben_(18f)"] = "florbetaben_(18F)"
ultimate_all_pre$NewName[ultimate_all_pre$API=="aripiprazol_lauroxil"] = "aripiprazole"
ultimate_all_pre$NewName[ultimate_all_pre$API=="florbetapir_(18f)"] = "florbetapir_(18F)"
ultimate_all_pre$NewName[ultimate_all_pre$API=="clorazepic_acid"] = "potassium_clorazepate"
ultimate_all_pre$NewName[ultimate_all_pre$API=="angiotensin_ii"] = "angiotensin_II"
ultimate_all_pre$NewName[ultimate_all_pre$API=="moricizine"] = "moracizine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="bretylium"] = "bretylium_tosilate"
ultimate_all_pre$NewName[ultimate_all_pre$API=="beclometasone_dipropionate"] = "beclometasone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="flutemetamol_(18f)"] = "flutemetamol_(18F)"
ultimate_all_pre$NewName[ultimate_all_pre$API=="human_interferon_beta"] = "interferon_beta_natural"
ultimate_all_pre$NewName[ultimate_all_pre$API=="rocuronium"] = "rocuronium_bromide"
ultimate_all_pre$NewName[ultimate_all_pre$API=="protamine_sulfate"] = "protamine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="desoxycorticosterone"] = "desoxycortone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="atracurium_besylate"] = "atracurium"
ultimate_all_pre$NewName[ultimate_all_pre$API=="antithrombin_iii_human"] = "antithrombin_III"
ultimate_all_pre$NewName[ultimate_all_pre$API=="ethynodiol_diacetate"] = "etynodiol"
ultimate_all_pre$NewName[ultimate_all_pre$API=="doxacurium"] = "doxacurium_chloride"
ultimate_all_pre$NewName[ultimate_all_pre$API=="estrone_sulfate"] = "estrone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="norgestimate"] = "norgestimate_and_estrogen"
ultimate_all_pre$NewName[ultimate_all_pre$API=="prednisolone_sulfate"] = "prednisolone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="testosterone_cypionate"] = "testosterone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="norethynodrel"] = "noretynodrel_and_estrogen"
ultimate_all_pre$NewName[ultimate_all_pre$API=="echothiophate"] = "ecothiopate"
ultimate_all_pre$NewName[ultimate_all_pre$API=="zimelidine"] = "zimeldine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="nandrolone_phenpropionate"] = "nandrolone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="methotrimeprazine"] = "levomepromazine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="testosterone_enanthate"] = "testosterone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="mivacurium"] = "mivacurium_chloride"
ultimate_all_pre$NewName[ultimate_all_pre$API=="dihydro-alpha-ergocryptine"] = "dihydroergocryptine_mesylate"
ultimate_all_pre$NewName[ultimate_all_pre$API=="betamethasone_phosphate"] = "betamethasone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="testosterone_undecanoate"] = "testosterone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="methscopolamine_bromide"] = "methylscopolamine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="racepinephrine"] = "epinephrine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="pipecuronium"] = "pipecuronium_bromide"
ultimate_all_pre$NewName[ultimate_all_pre$API=="flumethasone"] = "flumetasone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="hydroxyprogesterone_caproate"] = "hydroxyprogesterone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="testosterone_propionate"] = "testosterone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="sulthiame"] = "sultiame"
ultimate_all_pre$NewName[ultimate_all_pre$API=="chlorpheniramine"] = "dexchlorpheniramine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="tiotropium"] = "tiotropium_bromide"
ultimate_all_pre$NewName[ultimate_all_pre$API=="fluoroestradiol_f-18"] = "fluoroestradiol_(18F)"
ultimate_all_pre$NewName[ultimate_all_pre$API=="aclidinium"] = "aclidinium_bromide"
ultimate_all_pre$NewName[ultimate_all_pre$API=="umeclidinium"] = "umeclidinium_bromide"
ultimate_all_pre$NewName[ultimate_all_pre$API=="dicyclomine"] = "dicycloverine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="mometasone_furoate"] = "mometasone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="thiothixene"] = "tiotixene"
ultimate_all_pre$NewName[ultimate_all_pre$API=="nandrolone_decanoate"] = "nandrolone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="methylnaltrexone"] = "methylnaltrexone_bromide"
ultimate_all_pre$NewName[ultimate_all_pre$API=="clidinium"] = "clidinium_and_psycholeptics"
ultimate_all_pre$NewName[ultimate_all_pre$API=="antipyrine"] = "phenazone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="cortisone_acetate"] = "cortisone"
ultimate_all_pre$NewName[ultimate_all_pre$API=="metamizole"] = "metamizole_sodium"
ultimate_all_pre$NewName[ultimate_all_pre$API=="synthetic_conjugated_estrogens,_a"] = "conjugated_estrogens"
ultimate_all_pre$NewName[ultimate_all_pre$API=="synthetic_conjugated_estrogens,_b"] = "conjugated_estrogens"
ultimate_all_pre$NewName[ultimate_all_pre$API=="trolamine_salicylate"] = "trolamine"
ultimate_all_pre$NewName[ultimate_all_pre$API=="vilanterol"] = "vilanterol,_umeclidinium_bromide_and_fluticasone_furoate"
ultimate_all_pre$NewName[ultimate_all_pre$API=="loteprednol_etabonate"] = "loteprednol"
ultimate_all_pre$NewName[ultimate_all_pre$API=="prednisolone_acetate"] = "prednisolone"

# Pharmacological
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="nadh"] = "nicotinamide"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="drotrecogin_alfa"] = "drotrecogin_alfa_(activated)"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="vitamin_a"] = "retinol_(vit_A)"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="insulin_human"] = "insulin_(human)"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="copper"] = "copper_sulfate" # most appropriate
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="pyridoxine"] = "pyridoxine_(vit_B6)"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="somatotropin"] = "somatropin"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="vitamin_e"] = "tocopherol_(vit_E)"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="telotristat_ethyl"] = "telotristat"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="amphetamine"] = "amfetamine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="glyburide"] = "glibenclamide"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="thiamine"] = "thiamine_(vit_B1)"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="leuprolide"] = "leuprorelin"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="zinc"] = "zinc_sulfate" # same as copper
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="adenosine_phosphate"] = "adenosine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="muromonab"] = "muromonab-CD3"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="ardeparin"] = "enoxaparin"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="liothyronine"] = "liothyronine_sodium"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="secretin_human"] = "secretin"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="pentosan_polysulfate"] = "pentosan_polysulfate_sodium"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="coagulation_factor_viia_recombinant_human"] = "coagulation_factor_VIIa"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="cyclosporine"] = "ciclosporin"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="botulinum_toxin_type_a"] = "botulinum_toxin"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="salmon_calcitonin"] = "calcitonin_(salmon_synthetic)"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="interferon_gamma-1b"] = "interferon_gamma"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="niacin"] = "nicotinic_acid"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="capromab_pendetide"] = "indium_(111In)_capromab_pendetide"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="fluticasone_propionate"] = "fluticasone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="insulin_pork"] = "insulin_(pork)"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="lithium_citrate"] = "lithium" # only two lithium salts in the DrugBank files
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="coagulation_factor_ix_(recombinant)"] = "coagulation_factor_IX"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="medroxyprogesterone_acetate"] = "medroxyprogesterone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="sipuleucel-t"] = "sipuleucel-T"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="ferrous_sulfate_anhydrous"] = "ferrous_sulfate"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="vitamin_d"] = "colecalciferol"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="indomethacin"] = "indometacin"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="botulinum_toxin_type_b"] = "botulinum_toxin"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="omega-3-carboxylic_acids"] = "omega-3-triglycerides_incl._other_esters_and_acids"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="dyphylline"] = "diprophylline"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="ezogabine"] = "retigabine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="glatiramer"] = "glatiramer_acetate"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="cysteamine"] = "mercaptamine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="acetaminophen"] = "paracetamol"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="pentamidine"] = "pentamidine_isethionate"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="rutin"] = "rutoside"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="anthralin"] = "dithranol"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="dalfampridine"] = "fampridine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="mycophenolate_mofetil"] = "mycophenolic_acid"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="fluciclovine_(18f)"] = "fluciclovine_(18F)"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="ouabain"] = "g-strophanthin"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="eslicarbazepine_acetate"] = "eslicarbazepine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="human_calcitonin"] = "calcitonin_(human_synthetic)"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="ergoloid_mesylate"] = "ergoloid_mesylates"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="butarbital"] = "butobarbital"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="estradiol_acetate"] = "estradiol"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="xanthinol"] = "xantinol_nicotinate"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="cholecystokinin"] = "pancreozymin_(cholecystokinin)"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="olmesartan"] = "olmesartan_medoxomil"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="carboprost_tromethamine"] = "carboprost"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="ibritumomab_tiuxetan"] = "ibritumomab_tiuxetan_(90Y)"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="dexamethasone_acetate"] = "dexamethasone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="sodium_ferric_gluconate_complex"] = "ferrous_gluconate"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="vasopressin"] = "vasopressin_(argipressin)"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="thimerosal"] = "thiomersal"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="chlorthalidone"] = "chlortalidone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="levothyroxine"] = "levothyroxine_sodium"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="meperidine"] = "pethidine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="isoflurophate"] = "fluostigmine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="nitroglycerin"] = "glyceryl_trinitrate"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="isoetharine"] = "isoetarine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="phylloquinone"] = "isoetarine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="lithium_carbonate"] = "lithium" # only two lithium salts in the DrugBank files
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="coagulation_factor_vii"] = "coagulation_factor_VIIa"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="podofilox"] = "podophyllotoxin"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="human_thrombin"] = "thrombin"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="phylloquinone"] = "phytomenadione"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="gamma-aminobutyric_acid"] = "aminobutyric_acid"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="difluocortolone"] = "diflucortolone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="meclizine"] = "meclozine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="omega-3-_acid_ethyl_esters"] = "omega-3-triglycerides_incl._other_esters_and_acids"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="cyproterone_acetate"] = "cyproterone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="fenofibric_acid"] = "fenofibrate"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="estradiol_benzoate"] = "estradiol"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="hydrocortisone_acetate"] = "hydrocortisone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="tositumomab"] = "tositumomab/iodine_(131I)_tositumomab"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="orphenadrine"] = "orphenadrine_(citrate)"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="dinoprost_tromethamine"] = "dinoprost"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="megestrol_acetate"] = "megestrol"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="pentaerythritol_tetranitrate"] = "pentaerithrityl_tetranitrate"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="histamine"] = "histamine_phosphate"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="hyaluronidase_(human_recombinant)"] = "hyaluronidase"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="dextroamphetamine"] = "dexamfetamine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="methyldopa"] = "methyldopa_(racemic)"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="dipivefrin"] = "dipivefrine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="ibandronate"] = "ibandronic_acid"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="ethanolamine_oleate"] = "monoethanolamine_oleate"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="omega-3_fatty_acids"] = "omega-3-triglycerides_incl._other_esters_and_acids"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="estradiol_cypionate"] = "estradiol"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="protein_c"] = "protein_C"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="von_willebrand_factor_human"] = "von_Willebrand_factor"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="candesartan_cilexetil"] = "candesartan"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="liotrix"] = "combinations_of_levothyroxine_and_liothyronine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="hydroxyurea"] = "hydroxycarbamide"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="methimazole"] = "thiamazole"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="ipratropium"] = "ipratropium_bromide"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="metyrosine"] = "metirosine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="cholecalciferol"] = "colecalciferol"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="tegafur-uracil"] = "tegafur"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="coagulation_factor_ix_human"] = "coagulation_factor_IX"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="tromethamine"] = "trometamol"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="estradiol_dienanthate"] = "estradiol"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="hydrocortisone_cypionate"] = "hydrocortisone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="interferon_alfa-2a,_recombinant"] = "interferon_alfa-2a"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="clobetasol_propionate"] = "clobetasol"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="gallamine_triethiodide"] = "gallamine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="norelgestromin"] = "norelgestromin_and_ethinylestradiol"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="estradiol_valerate"] = "estradiol"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="hydrocortisone_phosphate"] = "hydrocortisone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="guaiacol"] = "guaiacolsulfonate"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="trastuzumab_deruxtecan"] = "trastuzumab"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="gabapentin_enacarbil"] = "gabapentin"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="nylidrin"] = "buphenine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="hydrocortisone_valerate"] = "hydrocortisone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="proparacaine"] = "proxymetacaine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="hydrocortisone_probutate"] = "hydrocortisone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="chorionic_gonadotropin_(human)"] = "choriogonadotropin_alfa"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="methsuximide"] = "mesuximide"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="lumacaftor"] = "ivacaftor_and_lumacaftor"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="florbetaben_(18f)"] = "florbetaben_(18F)"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="aripiprazol_lauroxil"] = "aripiprazole"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="florbetapir_(18f)"] = "florbetapir_(18F)"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="clorazepic_acid"] = "potassium_clorazepate"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="angiotensin_ii"] = "angiotensin_II"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="moricizine"] = "moracizine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="bretylium"] = "bretylium_tosilate"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="beclometasone_dipropionate"] = "beclometasone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="flutemetamol_(18f)"] = "flutemetamol_(18F)"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="human_interferon_beta"] = "interferon_beta_natural"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="rocuronium"] = "rocuronium_bromide"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="protamine_sulfate"] = "protamine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="desoxycorticosterone"] = "desoxycortone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="atracurium_besylate"] = "atracurium"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="antithrombin_iii_human"] = "antithrombin_III"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="ethynodiol_diacetate"] = "etynodiol"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="doxacurium"] = "doxacurium_chloride"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="estrone_sulfate"] = "estrone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="norgestimate"] = "norgestimate_and_estrogen"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="prednisolone_sulfate"] = "prednisolone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="testosterone_cypionate"] = "testosterone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="norethynodrel"] = "noretynodrel_and_estrogen"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="echothiophate"] = "ecothiopate"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="zimelidine"] = "zimeldine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="nandrolone_phenpropionate"] = "nandrolone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="methotrimeprazine"] = "levomepromazine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="testosterone_enanthate"] = "testosterone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="mivacurium"] = "mivacurium_chloride"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="dihydro-alpha-ergocryptine"] = "dihydroergocryptine_mesylate"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="betamethasone_phosphate"] = "betamethasone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="testosterone_undecanoate"] = "testosterone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="methscopolamine_bromide"] = "methylscopolamine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="racepinephrine"] = "epinephrine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="pipecuronium"] = "pipecuronium_bromide"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="flumethasone"] = "flumetasone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="hydroxyprogesterone_caproate"] = "hydroxyprogesterone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="testosterone_propionate"] = "testosterone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="sulthiame"] = "sultiame"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="chlorpheniramine"] = "dexchlorpheniramine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="tiotropium"] = "tiotropium_bromide"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="fluoroestradiol_f-18"] = "fluoroestradiol_(18F)"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="aclidinium"] = "aclidinium_bromide"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="umeclidinium"] = "umeclidinium_bromide"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="dicyclomine"] = "dicycloverine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="mometasone_furoate"] = "mometasone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="thiothixene"] = "tiotixene"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="nandrolone_decanoate"] = "nandrolone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="methylnaltrexone"] = "methylnaltrexone_bromide"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="clidinium"] = "clidinium_and_psycholeptics"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="antipyrine"] = "phenazone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="cortisone_acetate"] = "cortisone"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="metamizole"] = "metamizole_sodium"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="synthetic_conjugated_estrogens,_a"] = "conjugated_estrogens"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="synthetic_conjugated_estrogens,_b"] = "conjugated_estrogens"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="trolamine_salicylate"] = "trolamine"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="vilanterol"] = "vilanterol,_umeclidinium_bromide_and_fluticasone_furoate"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="loteprednol_etabonate"] = "loteprednol"
ultimate_pharmacological_pre$NewName[ultimate_pharmacological_pre$API=="prednisolone_acetate"] = "prednisolone"

# Joining again
ultimate_pharmacological = ultimate_pharmacological_pre %>%
  dplyr::select(Name, Gene.Symbol, EntrezGene.ID, API, NewName) %>%
  distinct() %>%
  left_join(ATC, by = c("NewName" = "API"))
ultimate_all = ultimate_all_pre %>%
  dplyr::select(Name, Gene.Symbol, EntrezGene.ID, API, NewName) %>%
  distinct() %>%
  left_join(ATC, by = c("NewName" = "API"))

# Allocate previously unmatched TKIs (with the risk of misclassifying some immunosuppressant ones:
for (i in 1:nrow(ultimate_all)){
  if (str_sub(ultimate_all$API[i], (nchar(ultimate_all$API[i])-3),
              nchar(ultimate_all$API[i])) %in% c("enib", "anib", "inib") && 
      is.na(ultimate_all$NewName[i])==TRUE){
    ultimate_all$Name0[i] = "ANTINEOPLASTIC AND IMMUNOMODULATING AGENTS"
    ultimate_all$Name1[i] = "ANTINEOPLASTIC AGENTS"
    ultimate_all$Name2[i] = "OTHER ANTINEOPLASTIC AGENTS"
    ultimate_all$Name3[i] = "Protein kinase inhibitors"
  }
}

for (i in 1:nrow(ultimate_pharmacological)){
  if (str_sub(ultimate_pharmacological$API[i], (nchar(ultimate_pharmacological$API[i])-3),
              nchar(ultimate_pharmacological$API[i])) %in% c("enib", "anib", "inib") && 
      is.na(ultimate_pharmacological$NewName[i])==TRUE){
    ultimate_pharmacological$Name0[i] = "ANTINEOPLASTIC AND IMMUNOMODULATING AGENTS"
    ultimate_pharmacological$Name1[i] = "ANTINEOPLASTIC AGENTS"
    ultimate_pharmacological$Name2[i] = "OTHER ANTINEOPLASTIC AGENTS"
    ultimate_pharmacological$Name3[i] = "Protein kinase inhibitors"
  }
}

# Finally label everything else as unmatched
for (i in 1:nrow(ultimate_all)){
  if (is.na(ultimate_all$NewName[i])==TRUE){
    ultimate_all$Name0[i] = "UNMATCHED"
    ultimate_all$Name1[i] = "UNMATCHED"
    ultimate_all$Name2[i] = "UNMATCHED"
    ultimate_all$Name3[i] = "UNMATCHED"
  }
}

for (i in 1:nrow(ultimate_pharmacological)){
  if (is.na(ultimate_pharmacological$NewName[i])==TRUE){
    ultimate_pharmacological$Name0[i] = "UNMATCHED"
    ultimate_pharmacological$Name1[i] = "UNMATCHED"
    ultimate_pharmacological$Name2[i] = "UNMATCHED"
    ultimate_pharmacological$Name3[i] = "UNMATCHED"
  }
}
rm(ultimate_pharmacological_pre, ultimate_all_pre, final_pharmacological, dups_actual,
   final_all, molten_all, molten_all_annot, molten_pharmacological, dups,
   molten_pharmacological_annot, unmatched_all, unmatched_pharmacological, i, j)

# Manually renaming the Name1 column levels of ultimate_pharmacological:
ultimate_pharmacological$Name0 = str_to_sentence(ultimate_pharmacological$Name0)
ultimate_pharmacological$Name0 = gsub(" ", "_", ultimate_pharmacological$Name0)
ultimate_pharmacological$Name0 = gsub("-", "_", ultimate_pharmacological$Name0)
ultimate_pharmacological$Name0 = gsub("___", "_", ultimate_pharmacological$Name0)
ultimate_pharmacological$Name0 = gsub("__", "_", ultimate_pharmacological$Name0)
ultimate_pharmacological$Name0 = gsub(",", "_", ultimate_pharmacological$Name0)
ultimate_pharmacological$Name0 = gsub("'", "_", ultimate_pharmacological$Name0)
ultimate_pharmacological$Name0 = gsub("Cardiovascular_System", "Cardiovascular", ultimate_pharmacological$Name0)
ultimate_pharmacological$Name0 = gsub("Antineoplastic_And_Immunomodulating_Agents", "Antineoplastics_and_Immunomodulating", ultimate_pharmacological$Name0)
ultimate_pharmacological$Name0 = gsub("Respiratory_System", "Respiratory", ultimate_pharmacological$Name0)
ultimate_pharmacological$Name0 = gsub("Unmatched", "Various", ultimate_pharmacological$Name0)
ultimate_pharmacological$Name0 = gsub("Genito_Urinary_System_And_Sex_Hormones", "Genito_urinary", ultimate_pharmacological$Name0)
ultimate_pharmacological$Name0 = gsub("Musculo_Skeletal_System", "Musculoskeletal", ultimate_pharmacological$Name0)
ultimate_pharmacological$Name0 = gsub("Alimentary_Tract_And_Metabolism", "GI_and_metabolism", ultimate_pharmacological$Name0)
ultimate_pharmacological$Name0 = gsub("Systemic_Hormonal_Preparations__Excl._Sex_Hormones_And_Insulins", "Hormonal_etc", ultimate_pharmacological$Name0)
ultimate_pharmacological$Name0 = gsub("Antiparasitic_Products__Insecticides_And_Repellents", "Antiparasitic_etc", ultimate_pharmacological$Name0)
ultimate_pharmacological$Name0 = gsub("Antiinfectives_For_Systemic_Use", "Antiinfectives", ultimate_pharmacological$Name0)

ultimate_pharmacological$API = gsub(",", "_", ultimate_pharmacological$API)
ultimate_pharmacological$API = gsub("-", "_", ultimate_pharmacological$API)
ultimate_pharmacological$API = str_to_sentence(ultimate_pharmacological$API)

##### Organize genes (drug targets) by networks #####
# We will focus on BioCarta and KEGG clustered terms (top 10 in each case). The plan is
# to make lists of up-/down-regulated genes from the clustering results. Each list
# will belong to a unique representative BioCarta/KEGG term. We will then create drug-gene
# links for each term. We will do this for each stage. We will also incorporate a miRNA
# ideogram on each plot to show interactions of these genes with miRNAs.

# A list of clustered results
clustered_results = list(clustered_results_stage_1, 
                         clustered_results_stage_2,
                         clustered_results_stage_3, 
                         clustered_results_stage_4,
                         clustered_results_blood)
names(clustered_results) = c("Stage_1", "Stage_2", "Stage_3", "Stage_4", "Blood")
rm(clustered_results_stage_1, clustered_results_stage_2,
   clustered_results_stage_3, clustered_results_stage_4,
   clustered_results_blood); gc()

Stage_1_BioCarta = list(); Stage_1_KEGG = list()
Stage_2_BioCarta = list(); Stage_2_KEGG = list()
Stage_3_BioCarta = list(); Stage_3_KEGG = list()
Stage_4_BioCarta = list(); Stage_4_KEGG = list()
Blood_BioCarta = list(); Blood_KEGG = list()

Stage_1_lists = list(Stage_1_BioCarta, Stage_1_KEGG); rm(Stage_1_BioCarta, Stage_1_KEGG)
Stage_2_lists = list(Stage_2_BioCarta, Stage_2_KEGG); rm(Stage_2_BioCarta, Stage_2_KEGG)
Stage_3_lists = list(Stage_3_BioCarta, Stage_3_KEGG); rm(Stage_3_BioCarta, Stage_3_KEGG)
Stage_4_lists = list(Stage_4_BioCarta, Stage_4_KEGG); rm(Stage_4_BioCarta, Stage_4_KEGG)
Blood_lists = list(Blood_BioCarta, Blood_KEGG); rm(Blood_BioCarta, Blood_KEGG)
names(Stage_1_lists) = names(Stage_2_lists) = names(Stage_3_lists) = names(Stage_4_lists) =
  names(Blood_lists) = c("BioCarta", "KEGG")

circos_path = list(Stage_1_lists, Stage_2_lists, Stage_3_lists, Stage_4_lists, Blood_lists)
names(circos_path) = c("Stage_1", "Stage_2", "Stage_3", "Stage_4", "Blood")
rm(Stage_1_lists, Stage_2_lists, Stage_3_lists, Stage_4_lists, Blood_lists); gc()

for (i in 1:5){ # 4 stages plus blood
  for (k in c("BioCarta", "KEGG")){ # databases of interest
    for (j in 1:10){ # top 10 representative terms
      up = scan(text = clustered_results[[i]][[k]]$Up_regulated[clustered_results[[i]][[k]]$Status == "Representative"][j], what = "")
      down = scan(text = clustered_results[[i]][[k]]$Down_regulated[clustered_results[[i]][[k]]$Status == "Representative"][j], what = "")
      
      if (length(up) > 1){
        str_sub(up[1:(length(up) - 1)], -1, -1) = ""
        checker_up = "yes"
      } else if (length(up) == 1) {
        checker_up = "yes"
      } else if (length(up) == 0) {
        checker_up = "not"
      }
      
      if (length(down) > 1){
        str_sub(down[1:(length(down) - 1)], -1, -1) = ""
        checker_down = "yes"
      } else if (length(down) == 1){
        checker_down = "yes"
      } else if (length(down) == 0){
        checker_down = "not"
      }
      
      if (checker_up == "yes" && checker_down == "yes"){
        circos_path[[i]][[k]][[j]] = sort(c(up, down))
      } else if (checker_up == "yes" && checker_down == "not") {
        circos_path[[i]][[k]][[j]] = sort(up)
      } else if (checker_up == "not" && checker_down == "yes") {
        circos_path[[i]][[k]][[j]] = sort(down)
      } else if (checker_up == "not" && checker_down == "not") {
        circos_path[[i]][[k]][[j]] = NA
      }
    }
    names(circos_path[[i]][[k]]) = clustered_results[[i]][[k]][["Term_Description"]][1:10]
  }
}

##### Circos #####
# Prepare the ideograms; which in this case will be networks of genes. We will prepare
# two karyotypes: one with BioCarta ideograms and one with KEGG ideograms. Drug ideograms
# separated by classes will be included in both karyotypes. Same for miRNAs.

# The above will be generated five times: 4 stages + blood

Circos_universe = list()
circos_stage_1 = list(); circos_stage_2 = list(); circos_stage_3 = list()
circos_stage_4 = list(); circos_blood = list()
circos_stage_1_BioCarta = list(); circos_stage_1_KEGG = list()
circos_stage_2_BioCarta = list(); circos_stage_2_KEGG = list()
circos_stage_3_BioCarta = list(); circos_stage_3_KEGG = list()
circos_stage_4_BioCarta = list(); circos_stage_4_KEGG = list()
circos_blood_BioCarta = list(); circos_blood_KEGG = list()
circos_stage_1 = list(circos_stage_1_BioCarta, circos_stage_1_KEGG)
circos_stage_2 = list(circos_stage_2_BioCarta, circos_stage_2_KEGG)
circos_stage_3 = list(circos_stage_3_BioCarta, circos_stage_3_KEGG)
circos_stage_4 = list(circos_stage_4_BioCarta, circos_stage_4_KEGG)
circos_blood = list(circos_blood_BioCarta, circos_blood_KEGG)
names(circos_stage_1) = names(circos_stage_2) = names(circos_stage_3) = 
  names(circos_stage_4) = names(circos_blood) = c("BioCarta", "KEGG")
Circos_universe = list(circos_stage_1, circos_stage_2, circos_stage_3,
                       circos_stage_4, circos_blood)

rm(circos_stage_1, circos_stage_2, circos_stage_3, circos_stage_4,
   circos_stage_1_BioCarta, circos_stage_1_KEGG, circos_stage_2_BioCarta,
   circos_stage_2_KEGG, circos_stage_3_BioCarta, circos_stage_3_KEGG,
   circos_stage_4_BioCarta, circos_stage_4_KEGG, circos_blood,
   circos_blood_BioCarta, circos_blood_KEGG)

names(Circos_universe) = c("Stage_1", "Stage_2", "Stage_3", "Stage_4", "Blood")

# miRNA predictions from Mienturnet #####
# One variable that we would like to enrich our genes with is the "degree", i.e.
# the number of miRNAs each gene is interacting with. We will do the same, but 
# in the opposite way for the miRNAs. We start with the genes (one needs to enable
# editing on all .xlsx files loaded here):

gene_results_mienturnet = list()
mienturnet_dirs = c("Stage_1", "Stage_2", "Stage_3", "Stage_4", "Blood")

for (i in 1:5){ # 4 stages + blood
  mienturnet = read.xlsx(paste0("Mienturnet/", mienturnet_dirs[i], 
                                "/", mienturnet_dirs[i], 
                                "_genes_as_columns_miRTarBase.xlsx"),
                         startRow = 2, colNames = TRUE)
  genes_mienturnet = data.frame(matrix(ncol = 2))
  colnames(genes_mienturnet) = c("Gene.Symbol", "Degree")
  for(j in 1:ncol(mienturnet)){
    degree = length(which(is.na(mienturnet[,j])==F))
    gene = colnames(mienturnet)[j]
    genes_mienturnet = rbind(genes_mienturnet, c(gene, degree))
  }
  rm(gene, degree, mienturnet)
  gene_results_mienturnet[[i]] = genes_mienturnet[2:nrow(genes_mienturnet),]
  rm(genes_mienturnet)
}

names(gene_results_mienturnet) = names(Circos_universe)

miRNA_results_mienturnet = list()

for (i in 1:5){ # 4 stages
  mienturnet = read.xlsx(paste0("Mienturnet/", mienturnet_dirs[i], 
                                "/", mienturnet_dirs[i], 
                                "_Mienturnet_Enrichment_results_miRTarBase.xlsx"),
                         startRow = 2, colNames = TRUE) %>%
    dplyr::select(microRNA, FDR, Number.of.interactions) %>%
    dplyr::rename(degree = Number.of.interactions, adjpval = FDR) %>%
    dplyr::filter(as.numeric(adjpval) < 0.05)
  miRNA_results_mienturnet[[i]] = mienturnet
  rm(mienturnet)
}

names(miRNA_results_mienturnet) = names(Circos_universe)

# COSMIC data #####
COSMIC_pre = read.csv("COSMIC_PDC_22022022.csv")
COSMIC = COSMIC_pre %>%
  dplyr::mutate(Ratio = Mutated.samples/Samples.tested) %>%
  dplyr::filter(Ratio > 0.01) %>% # filter for mutations with frequency of over 1%
  dplyr::filter(!grepl("\\.", Gene.name)) %>%
  dplyr::filter(!grepl("_", Gene.name))
COSMIC = COSMIC[order(COSMIC$Ratio, decreasing = TRUE), ] %>%
  distinct(Gene.name, .keep_all = TRUE); rm(COSMIC_pre)

# Common Gene symbols between COSMIC and the official gene symbols from the 
# org.Hs.eg.db database:

concordance = length(which(COSMIC$Gene.name %in% official_df$Gene.Symbol == TRUE))/nrow(COSMIC)
mismatch = COSMIC$Gene.name[which(COSMIC$Gene.name %in% official_df$Gene.Symbol == FALSE)]
a_concordance = length(which(mismatch %in% ID_Map$probe == TRUE))/length(mismatch)
discordance = mismatch[-which(mismatch %in% ID_Map$probe == TRUE)]
discordance

# discordance contains 19 genes some of which can be matched to Entrez ID's
# (https://www.genenames.org/tools/search/#!/?query=)

COSMIC$probe = COSMIC$Gene.name
COSMIC$probe[COSMIC$Gene.name == "36951"] = "RPL37P10"  # pseudogene
COSMIC$probe[COSMIC$Gene.name == "40603"] = "UPK1A-AS1" # long, non-coding RNA
COSMIC$probe[COSMIC$Gene.name == "37226"] = "PHGR1"     # coding
COSMIC$probe[COSMIC$Gene.name == "38047"] = "UQCRHP4"   # pseudogene
COSMIC$probe[COSMIC$Gene.name == "37681"] = "unclear" 
COSMIC$probe[COSMIC$Gene.name == "41833"] = "MIR4683"   # miRNA
COSMIC$probe[COSMIC$Gene.name == "40238"] = "unclear"
COSMIC$probe[COSMIC$Gene.name == "39326"] = "SNRPGP7"   # pseudogene
COSMIC$probe[COSMIC$Gene.name == "40787"] = "unclear"
COSMIC$probe[COSMIC$Gene.name == "38777"] = "USP9YP33"  # pseudogene
COSMIC$probe[COSMIC$Gene.name == "38961"] = "MIR3622B"  # miRNA
COSMIC$probe[COSMIC$Gene.name == "40422"] = "unclear"
COSMIC$probe[COSMIC$Gene.name == "39508"] = "ATP5MC1P5" # pseudogene
COSMIC$probe[COSMIC$Gene.name == "38412"] = "unclear"
COSMIC$probe[COSMIC$Gene.name == "39142"] = "unclear"
COSMIC$probe[COSMIC$Gene.name == "37865"] = "GLYATL1B"  # coding
COSMIC$probe[COSMIC$Gene.name == "37500"] = "RNVU1-7"   # snRNA

# further corrections
COSMIC$probe[COSMIC$Gene.name == "MPP6"] = "PALS2"    
# alias for two official symbols, 
# but MPHOSPH6 exists in COSMIC_pre so we believe this is for PALS2
COSMIC$probe[COSMIC$Gene.name == "DUSP27"] = "unclear" 
# unclear if it must match to STYXL2 or DUSP29

ID_Map = ID_Map %>% dplyr::rename(probe = Gene.Symbol)
COSMIC_annot = COSMIC %>%
  left_join(ID_Map, by = "probe") %>%
  dplyr::filter(!probe == "unclear") %>%
  left_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(EntrezGene.ID, Gene.Symbol, Mutated.samples, Samples.tested, Ratio) %>%
  distinct()

# Cancer gene drivers
census = read.csv("COSMIC_census_09_02_2022.csv") %>%
  dplyr::rename(EntrezGene.ID = Entrez.GeneId, Gene.Symbol_COSMIC = Gene.Symbol) %>%
  dplyr::select(EntrezGene.ID, Gene.Symbol_COSMIC, Name, everything())
census$EntrezGene.ID = as.character(census$EntrezGene.ID)
census_annot = census %>% inner_join(official_df, by = "EntrezGene.ID")
census_annot_pancr = census_annot %>%
  dplyr::filter(grepl("pancr", Tumour.Types.Somatic.) | grepl("pancr", Tumour.Types.Germline.))
general_drivers = census_annot$Gene.Symbol
PDAC_drivers = census_annot_pancr$Gene.Symbol

# Basic Circos files Gene labels, Drug labels #####
# BC = BioCarta, KG = KEGG

for (i in 1:5){ # 4 stages + blood
  BC = melt(circos_path[[i]][["BioCarta"]]) %>%
    dplyr::rename(PathChr = L1, Gene.Symbol = value)
  BC$PathChr = gsub("BIOCARTA_", "", BC$PathChr)
  BC$PathChr = gsub("_PATHWAY", "", BC$PathChr)
  str_sub(BC$PathChr, 0, 0) = "hs"
  BC$GStart[1] = 0
  BC$GEnd[1] = 5
  
  for(k in 2:nrow(BC)){
    if(BC$PathChr[k] == BC$PathChr[k-1]){
      BC$GStart[k] = BC$GStart[k-1] + 5
      BC$GEnd[k]   = BC$GEnd[k-1] + 5
    } else {
      BC$GStart[k] = 0
      BC$GEnd[k]   = 5
    }
  }
  rm(k)
  
  BC = BC %>% left_join(sig_DGEA[[i]], by = c("Gene.Symbol" = "Gene.Symbol_pre")) %>%
    mutate(GVars = paste0("logfc=", logFC, ",adjpval=", adj.P.Val)) %>%
    dplyr::select(PathChr, GStart, GEnd, Gene.Symbol, GVars)
  
  KG = melt(circos_path[[i]][["KEGG"]]) %>%
    dplyr::rename(PathChr = L1, Gene.Symbol = value)
  KG$PathChr = gsub(" ", "_", KG$PathChr)
  KG$PathChr = gsub("-", "_", KG$PathChr)
  KG$PathChr = gsub("___", "_", KG$PathChr)
  KG$PathChr = gsub("__", "_", KG$PathChr)
  str_sub(KG$PathChr, 0, 0) = "hs"
  KG$GStart[1] = 0
  KG$GEnd[1] = 5
  
  for(k in 2:nrow(KG)){
    if(KG$PathChr[k] == KG$PathChr[k-1]){
      KG$GStart[k] = KG$GStart[k-1] + 5
      KG$GEnd[k]   = KG$GEnd[k-1] + 5
    } else {
      KG$GStart[k] = 0
      KG$GEnd[k]   = 5
    }
  }
  rm(k)
  
  KG = KG %>% left_join(sig_DGEA[[i]], by = c("Gene.Symbol" = "Gene.Symbol_pre")) %>%
    mutate(GVars = paste0("logfc=", logFC, ",adjpval=", adj.P.Val)) %>%
    dplyr::select(PathChr, GStart, GEnd, Gene.Symbol, GVars)
  
  
  # Drug Labels (only drugs that interact with our list of genes BC$Gene.Symbol and 
  # KG$Symbol). BioCarta:
  final_pharmacological_BC = BC %>% 
    inner_join(ultimate_pharmacological %>% dplyr::select(Gene.Symbol, API, Name0),  
               by = "Gene.Symbol")
  
  BC_Drugs = final_pharmacological_BC %>% dplyr::select(API, Name0) %>%
    dplyr::rename(DChr = Name0) %>%
    distinct()
  BC_Drugs = BC_Drugs[order(BC_Drugs$DChr, BC_Drugs$API),]
  
  str_sub(BC_Drugs$DChr, 0, 0) = "hs"
  BC_Drugs$DStart[1] = 0
  BC_Drugs$DEnd[1] = 5
  
  for(j in 2:nrow(BC_Drugs)){
    if(BC_Drugs$DChr[j] == BC_Drugs$DChr[j-1]){
      BC_Drugs$DStart[j] = BC_Drugs$DStart[j-1] + 5
      BC_Drugs$DEnd[j]   = BC_Drugs$DEnd[j-1] + 5
    } else {
      BC_Drugs$DStart[j] = 0
      BC_Drugs$DEnd[j]   = 5
    }
  }
  rm(j)
  BC_Drugs = BC_Drugs %>%
    dplyr::select(DChr, DStart, DEnd, API)
  
  # KEGG
  final_pharmacological_KG = KG %>% 
    inner_join(ultimate_pharmacological %>% dplyr::select(Gene.Symbol, API, Name0),  
               by = "Gene.Symbol")
  
  KG_Drugs = final_pharmacological_KG %>% dplyr::select(API, Name0) %>%
    dplyr::rename(DChr = Name0) %>%
    distinct()
  KG_Drugs = KG_Drugs[order(KG_Drugs$DChr, KG_Drugs$API),]
  
  str_sub(KG_Drugs$DChr, 0, 0) = "hs"
  KG_Drugs$DStart[1] = 0
  KG_Drugs$DEnd[1] = 5
  
  for(j in 2:nrow(KG_Drugs)){
    if(KG_Drugs$DChr[j] == KG_Drugs$DChr[j-1]){
      KG_Drugs$DStart[j] = KG_Drugs$DStart[j-1] + 5
      KG_Drugs$DEnd[j]   = KG_Drugs$DEnd[j-1] + 5
    } else {
      KG_Drugs$DStart[j] = 0
      KG_Drugs$DEnd[j]   = 5
    }
  }
  rm(j)
  KG_Drugs = KG_Drugs %>%
    dplyr::select(DChr, DStart, DEnd, API)
  
  # Join with gene_results_mienturnet to add miRNA degree to the genes
  BC = BC %>% left_join(gene_results_mienturnet[[i]], by = "Gene.Symbol")
  zero_degree_index = which(is.na(BC$Degree) == TRUE)
  BC$Degree[zero_degree_index] = "0"; rm(zero_degree_index)
  str_sub(BC$Degree, 0, 0) = "degree="
  BC$GVars_new = paste0(BC$GVars, ",", BC$Degree)
  BC = BC %>% dplyr::select(-GVars, -Degree) %>%
    dplyr::rename(GVars = GVars_new)
  
  KG = KG %>% left_join(gene_results_mienturnet[[i]], by = "Gene.Symbol")
  zero_degree_index = which(is.na(KG$Degree) == TRUE)
  KG$Degree[zero_degree_index] = "0"; rm(zero_degree_index)
  str_sub(KG$Degree, 0, 0) = "degree="
  KG$GVars_new = paste0(KG$GVars, ",", KG$Degree)
  KG = KG %>% dplyr::select(-GVars, -Degree) %>%
    dplyr::rename(GVars = GVars_new)
  
  # Preparing mutation histogram and driver glyph files
  # Mutations
  BC_mutation = BC %>% dplyr::select(-GVars) %>%
    left_join(COSMIC_annot, by = c("Gene.Symbol")) %>%
    dplyr::select(-EntrezGene.ID)
  zero_mut_index = which(is.na(BC_mutation$Ratio) == TRUE)
  BC_mutation$Ratio[zero_mut_index] = 0; rm(zero_mut_index)
  BC_mutation$Ratio = log(100*BC_mutation$Ratio + 1)
  BC_mutation = BC_mutation %>%
    mutate(MVars = paste0("gene=", Gene.Symbol, ",sampmut=", Mutated.samples,
                          ",samptest=", Samples.tested)) %>%
    dplyr::select(-Gene.Symbol, -Mutated.samples, -Samples.tested) %>%
    dplyr::select(PathChr, GStart, GEnd, Ratio, MVars) %>%
    dplyr::rename(logRatio = Ratio)
  
  KG_mutation = KG %>% dplyr::select(-GVars) %>%
    left_join(COSMIC_annot, by = c("Gene.Symbol")) %>%
    dplyr::select(-EntrezGene.ID)
  zero_mut_index = which(is.na(KG_mutation$Ratio) == TRUE)
  KG_mutation$Ratio[zero_mut_index] = 0; rm(zero_mut_index)
  KG_mutation$Ratio = log(100*KG_mutation$Ratio + 1)
  KG_mutation = KG_mutation %>%
    mutate(MVars = paste0("gene=", Gene.Symbol, ",sampmut=", Mutated.samples,
                          ",samptest=", Samples.tested)) %>%
    dplyr::select(-Gene.Symbol, -Mutated.samples, -Samples.tested) %>%
    dplyr::select(PathChr, GStart, GEnd, Ratio, MVars) %>%
    dplyr::rename(logRatio = Ratio)
  
  # Drivers
  BC_drivers = BC %>% dplyr::select(-GVars)
  gdriver_index = which(BC_drivers$Gene.Symbol %in% general_drivers)
  pdriver_index = which(BC_drivers$Gene.Symbol %in% PDAC_drivers)
  BC_drivers$gdriver = "no"; BC_drivers$pdriver = "no"
  BC_drivers$gdriver[gdriver_index] = "yes"; BC_drivers$pdriver[pdriver_index] = "yes"
  rm(gdriver_index, pdriver_index)
  BC_drivers$glyph_color = "white"
  BC_drivers$glyph_color[BC_drivers$gdriver=="yes"] = "vdred_a15"
  BC_drivers$glyph_color[BC_drivers$pdriver=="yes"] = "vdgreen_a15"
  BC_drivers = BC_drivers %>%
    mutate(DrVars = paste0("gdriver=", gdriver, ",pdriver=", pdriver, 
                           ",glyphcolor=", glyph_color)) %>%
    dplyr::select(-gdriver, -pdriver, -glyph_color)
  
  KG_drivers = KG %>% dplyr::select(-GVars)
  gdriver_index = which(KG_drivers$Gene.Symbol %in% general_drivers)
  pdriver_index = which(KG_drivers$Gene.Symbol %in% PDAC_drivers)
  KG_drivers$gdriver = "no"; KG_drivers$pdriver = "no"
  KG_drivers$gdriver[gdriver_index] = "yes"; KG_drivers$pdriver[pdriver_index] = "yes"
  rm(gdriver_index, pdriver_index)
  KG_drivers$glyph_color = "white"
  KG_drivers$glyph_color[KG_drivers$gdriver=="yes"] = "vdred_a15"
  KG_drivers$glyph_color[KG_drivers$pdriver=="yes"] = "vdgreen_a15"
  KG_drivers = KG_drivers %>%
    mutate(DrVars = paste0("gdriver=", gdriver, ",pdriver=", pdriver, 
                           ",glyphcolor=", glyph_color)) %>%
    dplyr::select(-gdriver, -pdriver, -glyph_color)
  
  # Create histogram files for the genes for -log10(adjpval), degree, logfc:
  # logfc
  BC_histogram_logfc = BC %>% separate(GVars, sep = ",", into = c("logfc", "adjpval", "degree")) %>%
    dplyr::select(-Gene.Symbol, -degree) %>%
    dplyr::rename(histvars = adjpval)
  BC_histogram_logfc$logfc = gsub("logfc=", "", BC_histogram_logfc$logfc)
  
  KG_histogram_logfc = KG %>% separate(GVars, sep = ",", into = c("logfc", "adjpval", "degree")) %>%
    dplyr::select(-Gene.Symbol, -degree) %>%
    dplyr::rename(histvars = adjpval)
  KG_histogram_logfc$logfc = gsub("logfc=", "", KG_histogram_logfc$logfc)
  
  # adjpval
  BC_histogram_adjpval = BC %>% separate(GVars, sep = ",", into = c("logfc", "adjpval", "degree")) %>%
    dplyr::select(-Gene.Symbol, -logfc, -degree)
  BC_histogram_adjpval$adjpval = gsub("adjpval=", "", BC_histogram_adjpval$adjpval)
  nas = which(BC_histogram_adjpval$adjpval == "NA")
  if (length(nas) > 0){
    BC_histogram_adjpval$adjpval[nas] = 1
  }
  rm(nas)
  BC_histogram_adjpval$adjpval = -log10(as.numeric(BC_histogram_adjpval$adjpval))
  BC_histogram_adjpval = BC_histogram_adjpval %>%
    dplyr::rename(`neglogadjpval` = adjpval)
  
  KG_histogram_adjpval = KG %>% separate(GVars, sep = ",", into = c("logfc", "adjpval", "degree")) %>%
    dplyr::select(-Gene.Symbol, -logfc, -degree)
  KG_histogram_adjpval$adjpval = gsub("adjpval=", "", KG_histogram_adjpval$adjpval)
  nas = which(KG_histogram_adjpval$adjpval == "NA")
  if (length(nas) > 0){
    KG_histogram_adjpval$adjpval[nas] = 1
  }
  rm(nas)
  KG_histogram_adjpval$adjpval = -log10(as.numeric(KG_histogram_adjpval$adjpval))
  KG_histogram_adjpval = KG_histogram_adjpval %>%
    dplyr::rename(`neglogadjpval` = adjpval)
  
  # degree (how many miRNAs a gene interacts with)
  BC_histogram_degree = BC %>% separate(GVars, sep = ",", into = c("logfc", "adjpval", "degree")) %>%
    dplyr::select(-Gene.Symbol) %>%
    dplyr::mutate(histvars = paste0(logfc, ",", adjpval)) %>%
    dplyr::select(-logfc, -adjpval)
  BC_histogram_degree$degree = gsub("degree=", "", BC_histogram_degree$degree)
  
  KG_histogram_degree = KG %>% separate(GVars, sep = ",", into = c("logfc", "adjpval", "degree")) %>%
    dplyr::select(-Gene.Symbol) %>%
    dplyr::mutate(histvars = paste0(logfc, ",", adjpval)) %>%
    dplyr::select(-logfc, -adjpval)
  KG_histogram_degree$degree = gsub("degree=", "", KG_histogram_degree$degree)
  
  # miRNA locations and labels
  miRNA = miRNA_results_mienturnet[[i]] %>%
    mutate(miRNA_Chr = "hsmiRNA", miRstart = 0, miRend = 5,
           miRvars = paste0("mirdegree=", degree, ",miradjpval=", adjpval)) %>%
    dplyr::select(-degree, -adjpval)
  
  for(j in 2:nrow(miRNA)){
    if(miRNA$miRNA_Chr[j] == miRNA$miRNA_Chr[j-1]){
      miRNA$miRstart[j] = miRNA$miRstart[j-1] + 5
      miRNA$miRend[j]   = miRNA$miRend[j-1] + 5
    } else {
      miRNA$miRstart[j] = 0
      miRNA$miRend[j]   = 5
    }
  }
  rm(j)
  
  miRNA$microRNA = gsub("-", "_", miRNA$microRNA)
  miRNA = miRNA %>% dplyr::select(miRNA_Chr, miRstart, miRend, microRNA, miRvars)
  
  # miRNA degree histogram
  miRNA_histogram_degree = miRNA %>% separate(miRvars, sep = ",", into = c("mirdegree", "adjpval")) %>%
    dplyr::select(-microRNA) %>%
    dplyr::rename(histvars = adjpval)
  miRNA_histogram_degree$mirdegree = gsub("mirdegree=", "", miRNA_histogram_degree$mirdegree)
  
  # Links
  miRNA_gene = as.data.frame(t(read.xlsx(paste0("Mienturnet/", mienturnet_dirs[i], 
                                               "/", mienturnet_dirs[i], "_genes_as_columns_miRTarBase.xlsx"),
                                        startRow = 2, colNames = TRUE)))
  miRNA_gene$Gene.Symbol = rownames(miRNA_gene)
  miRNA_gene_map = melt(miRNA_gene, na.rm = TRUE, id.vars = "Gene.Symbol") %>%
    dplyr::select(-variable) %>%
    dplyr::rename(microRNA = value) %>%
    distinct()
  miRNA_gene_map$microRNA = gsub("-", "_", miRNA_gene_map$microRNA)
  
  # BioCarta
  BC_miRNA_links = miRNA %>% inner_join(miRNA_gene_map, by = "microRNA") %>%
    inner_join(BC, by = "Gene.Symbol") %>%
    dplyr::select(miRNA_Chr, miRstart, miRend, PathChr, GStart, GEnd, miRvars, 
                  GVars, Gene.Symbol, microRNA) %>%
    mutate(LVars = paste0(miRvars, ",", GVars, ",gene=", Gene.Symbol, ",mirna=", microRNA, ",drug=none")) %>%
    dplyr::select(-miRvars, -GVars, -Gene.Symbol, -microRNA)
  
  BC_drug_links = final_pharmacological_BC %>% 
    dplyr::select(-PathChr, -GStart, -GEnd, -GVars) %>%
    inner_join(BC_Drugs, by = "API") %>%
    inner_join(BC, by = "Gene.Symbol") %>%
    dplyr::select(PathChr, GStart, GEnd, DChr, DStart, DEnd, GVars, API, Gene.Symbol) %>%
    mutate(LVars = paste0("mirdegree=0,miradjpval=0,", GVars, 
                          ",gene=", Gene.Symbol, ",mirna=none", ",drug=", API)) %>%
    dplyr::select(-GVars, -API, -Gene.Symbol)
  colnames(BC_drug_links) = colnames(BC_miRNA_links) # name of the columns doesn't matter for Circos input
  BC_links = rbind(BC_miRNA_links, BC_drug_links); rm(BC_miRNA_links, BC_drug_links)
  
  # KEGG
  KG_miRNA_links = miRNA %>% inner_join(miRNA_gene_map, by = "microRNA") %>%
    inner_join(KG, by = "Gene.Symbol") %>%
    dplyr::select(miRNA_Chr, miRstart, miRend, PathChr, GStart, GEnd, miRvars, 
                  GVars, Gene.Symbol, microRNA) %>%
    mutate(LVars = paste0(miRvars, ",", GVars, ",gene=", Gene.Symbol, ",mirna=", microRNA, ",drug=none")) %>%
    dplyr::select(-miRvars, -GVars, -Gene.Symbol, -microRNA)
  
  rm(miRNA_gene, miRNA_gene_map)
  
  KG_drug_links = final_pharmacological_KG %>%
    dplyr::select(-PathChr, -GStart, -GEnd, -GVars) %>%
    inner_join(KG_Drugs, by = "API") %>%
    inner_join(KG, by = "Gene.Symbol") %>%
    dplyr::select(PathChr, GStart, GEnd, DChr, DStart, DEnd, GVars, API, Gene.Symbol) %>%
    mutate(LVars = paste0("mirdegree=0,miradjpval=0,", GVars, 
                          ",gene=", Gene.Symbol, ",mirna=none", ",drug=", API)) %>%
    dplyr::select(-GVars, -API, -Gene.Symbol)
  colnames(KG_drug_links) = colnames(KG_miRNA_links) # name of the columns doesn't matter for Circos input
  KG_links = rbind(KG_miRNA_links, KG_drug_links); rm(KG_miRNA_links, KG_drug_links)
  
  # Karyotypes
  # BioCarta
  BC_karyotype = data.frame(matrix(ncol = 6, 0))
  for(l in unique(BC_Drugs$DChr)){
    one = "chr"; two = "-"; three = l; four = gsub("hs", "", l); five = "0";
    six = max(BC_Drugs$DEnd[BC_Drugs$DChr == l])
    BC_karyotype = rbind(BC_karyotype, c(one, two, three, four, five, six))
    rm(one, two, three, four, five, six)
  } 
  rm(l)
  for(l in unique(BC$PathChr)){
    one = "chr"; two = "-"; three = l; four = gsub("hs", "", l); five = "0";
    six = max(BC$GEnd[BC$PathChr == l])
    BC_karyotype = rbind(BC_karyotype, c(one, two, three, four, five, six))
    rm(one, two, three, four, five, six)
  }
  rm(l)
  for(l in unique(miRNA$miRNA_Chr)){
    one = "chr"; two = "-"; three = l; four = gsub("hs", "", l); five = "0";
    six = max(miRNA$miRend[miRNA$miRNA_Chr == l])
    BC_karyotype = rbind(BC_karyotype, c(one, two, three, four, five, six))
    rm(one, two, three, four, five, six)
  }
  rm(l)
  BC_karyotype = BC_karyotype[2:nrow(BC_karyotype),]
  BC_karyotype$X7 = "vvdblue"
  
  # KEGG
  KG_karyotype = data.frame(matrix(ncol = 6, 0))
  for(l in unique(KG_Drugs$DChr)){
    one = "chr"; two = "-"; three = l; four = gsub("hs", "", l); five = "0";
    six = max(KG_Drugs$DEnd[KG_Drugs$DChr == l])
    KG_karyotype = rbind(KG_karyotype, c(one, two, three, four, five, six))
    rm(one, two, three, four, five, six)
  } 
  rm(l)
  for(l in unique(KG$PathChr)){
    one = "chr"; two = "-"; three = l; four = gsub("hs", "", l); five = "0";
    six = max(KG$GEnd[KG$PathChr == l])
    KG_karyotype = rbind(KG_karyotype, c(one, two, three, four, five, six))
    rm(one, two, three, four, five, six)
  }
  rm(l)
  for(l in unique(miRNA$miRNA_Chr)){
    one = "chr"; two = "-"; three = l; four = gsub("hs", "", l); five = "0";
    six = max(miRNA$miRend[miRNA$miRNA_Chr == l])
    KG_karyotype = rbind(KG_karyotype, c(one, two, three, four, five, six))
    rm(one, two, three, four, five, six)
  }
  rm(l)
  KG_karyotype = KG_karyotype[2:nrow(KG_karyotype),]
  KG_karyotype$X7 = "vvdblue"
  
  # Final lists
  Circos_universe[[i]][["BioCarta"]][[1]]  = BC_karyotype
  Circos_universe[[i]][["BioCarta"]][[2]]  = BC
  Circos_universe[[i]][["BioCarta"]][[3]]  = BC_Drugs
  Circos_universe[[i]][["BioCarta"]][[4]]  = miRNA
  Circos_universe[[i]][["BioCarta"]][[5]]  = BC_links
  Circos_universe[[i]][["BioCarta"]][[6]]  = BC_histogram_adjpval
  Circos_universe[[i]][["BioCarta"]][[7]]  = BC_histogram_logfc
  Circos_universe[[i]][["BioCarta"]][[8]]  = BC_histogram_degree
  Circos_universe[[i]][["BioCarta"]][[9]]  = BC_mutation
  Circos_universe[[i]][["BioCarta"]][[10]] = BC_drivers
  Circos_universe[[i]][["BioCarta"]][[11]] = miRNA_histogram_degree
  names(Circos_universe[[i]][["BioCarta"]]) = c("Karyotype", "Gene_Labels", "Drug_Labels", "miRNA_Labels",
                             "Links", "adjpval_histogram", "logFC_histogram", 
                             "degree_histogram", "Mutations", "Drivers", "miRNA_histogram")
  
  rm(BC_karyotype, BC, BC_Drugs, BC_links, BC_histogram_adjpval, BC_histogram_degree,
     BC_histogram_logfc, BC_mutation, BC_drivers)
  
  Circos_universe[[i]][["KEGG"]][[1]]  = KG_karyotype
  Circos_universe[[i]][["KEGG"]][[2]]  = KG
  Circos_universe[[i]][["KEGG"]][[3]]  = KG_Drugs
  Circos_universe[[i]][["KEGG"]][[4]]  = miRNA
  Circos_universe[[i]][["KEGG"]][[5]]  = KG_links
  Circos_universe[[i]][["KEGG"]][[6]]  = KG_histogram_adjpval
  Circos_universe[[i]][["KEGG"]][[7]]  = KG_histogram_logfc
  Circos_universe[[i]][["KEGG"]][[8]]  = KG_histogram_degree
  Circos_universe[[i]][["KEGG"]][[9]]  = KG_mutation
  Circos_universe[[i]][["KEGG"]][[10]] = KG_drivers
  Circos_universe[[i]][["KEGG"]][[11]] = miRNA_histogram_degree
  names(Circos_universe[[i]][["KEGG"]]) = c("Karyotype", "Gene_Labels", "Drug_Labels", "miRNA_Labels",
                             "Links", "adjpval_histogram", "logFC_histogram", 
                             "degree_histogram", "Mutations", "Drivers", "miRNA_histogram")
  
  rm(KG_karyotype, KG, KG_Drugs, miRNA, KG_links, KG_histogram_adjpval, KG_histogram_degree,
     KG_histogram_logfc, KG_mutation, KG_drivers, miRNA_histogram_degree)
  
  # Write out the .txt files that are necessary for the plots
  for(r in 1:length(Circos_universe[[i]][["BioCarta"]])){ # we will use BioCarta for both
    write.table(Circos_universe[[i]][["BioCarta"]][[r]], paste0("Circos/", names(Circos_universe)[[i]], 
                                                                "/BioCarta/", 
                                                                names(Circos_universe[[i]][["BioCarta"]])[r],
                                                                ".txt"), 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(Circos_universe[[i]][["KEGG"]][[r]], paste0("Circos/", names(Circos_universe)[[i]], 
                                                                "/KEGG/", 
                                                                names(Circos_universe[[i]][["KEGG"]])[r],
                                                                ".txt"), 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  gc()
}