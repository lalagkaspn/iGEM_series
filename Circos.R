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
                        "clean_miRNA_venn", "miRNA_venn", "sig_DGEA")))

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
clustered_results = list(clustered_results_stage_1, clustered_results_stage_2,
                         clustered_results_stage_3, clustered_results_stage_4)
names(clustered_results) = c("Stage_1", "Stage_2", "Stage_3", "Stage_4")
rm(clustered_results_stage_1, clustered_results_stage_2,
   clustered_results_stage_3, clustered_results_stage_4); gc()

Stage_1_BioCarta = list(); Stage_1_KEGG = list()
Stage_2_BioCarta = list(); Stage_2_KEGG = list()
Stage_3_BioCarta = list(); Stage_3_KEGG = list()
Stage_4_BioCarta = list(); Stage_4_KEGG = list()

Stage_1_lists = list(Stage_1_BioCarta, Stage_1_KEGG); rm(Stage_1_BioCarta, Stage_1_KEGG)
Stage_2_lists = list(Stage_2_BioCarta, Stage_2_KEGG); rm(Stage_2_BioCarta, Stage_2_KEGG)
Stage_3_lists = list(Stage_3_BioCarta, Stage_3_KEGG); rm(Stage_3_BioCarta, Stage_3_KEGG)
Stage_4_lists = list(Stage_4_BioCarta, Stage_4_KEGG); rm(Stage_4_BioCarta, Stage_4_KEGG)
names(Stage_1_lists) = names(Stage_1_lists) = names(Stage_1_lists) = names(Stage_1_lists) =
  c("BioCarta", "KEGG")

circos_path = list(Stage_1_lists, Stage_2_lists, Stage_3_lists, Stage_4_lists)
names(circos_path) = c("Stage_1", "Stage_2", "Stage_3", "Stage_4")
rm(Stage_1_lists, Stage_2_lists, Stage_3_lists, Stage_4_lists); gc()

for (i in 1:4){ # 4 stages
  for (k in c("BioCarta", "KEGG")){ # databases of interest
    for (j in 1:10){ # top 10 representative terms
      up = scan(text = clustered_results[[i]][[k]]$Up_regulated[j], what = "")
      down = scan(text = clustered_results[[i]][[k]]$Down_regulated[j], what = "")
      
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
        circos_path[[i]][[k]][[j]] = c(up, down)
      } else if (checker_up == "yes" && checker_down == "not") {
        circos_path[[i]][[k]][[j]] = up
      } else if (checker_up == "not" && checker_down == "yes") {
        circos_path[[i]][[k]][[j]] = down
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

# The above will be generated four times: once for each stage

Circos_universe = list()
circos_stage_1 = list(); circos_stage_2 = list(); circos_stage_3 = list(); circos_stage_4 = list()
circos_stage_1_BioCarta = list(); circos_stage_1_KEGG = list()
circos_stage_2_BioCarta = list(); circos_stage_2_KEGG = list()
circos_stage_3_BioCarta = list(); circos_stage_3_KEGG = list()
circos_stage_4_BioCarta = list(); circos_stage_4_KEGG = list()
circos_stage_1 = list(circos_stage_1_BioCarta, circos_stage_1_KEGG)
circos_stage_2 = list(circos_stage_2_BioCarta, circos_stage_2_KEGG)
circos_stage_3 = list(circos_stage_3_BioCarta, circos_stage_3_KEGG)
circos_stage_4 = list(circos_stage_4_BioCarta, circos_stage_4_KEGG)
names(circos_stage_1) = names(circos_stage_2) = names(circos_stage_3) = 
  names(circos_stage_4) = c("BioCarta", "KEGG")
Circos_universe = list(circos_stage_1, circos_stage_2, circos_stage_3, circos_stage_4)

rm(circos_stage_1, circos_stage_2, circos_stage_3, circos_stage_4,
   circos_stage_1_BioCarta, circos_stage_1_KEGG, circos_stage_2_BioCarta,
   circos_stage_2_KEGG, circos_stage_3_BioCarta, circos_stage_3_KEGG,
   circos_stage_4_BioCarta, circos_stage_4_KEGG)

names(Circos_universe) = c("Stage_1", "Stage_2", "Stage_3", "Stage_4")

# miRNA predictions from Mienturnet #####
# One variable that we would like to enrich our genes with is the "degree", i.e.
# the number of miRNAs each gene is interacting with. We will do the same, but 
# in the opposite way for the miRNAs. We start with the genes:

results_mienturnet = list()

for (i in 1:4){ # 4 stages
  mienturnet = read.xlsx(paste0("Mienturnet/Stage_", i, 
                                "/Stage_", i, "_genes_as_columns_miRTarBase.xlsx"),
                         startRow = 2, colNames = TRUE)
  genes_mienturnet = data.frame(matrix(ncol = 2))
  colnames(genes_mienturnet) = c("Gene.Symbol", "Degree")
  for(j in 1:ncol(mienturnet)){
    degree = length(which(is.na(mienturnet[,j])==F))
    gene = colnames(mienturnet)[j]
    genes_mienturnet = rbind(genes_mienturnet, c(gene, degree))
  }
  rm(gene, degree, mienturnet)
  results_mienturnet[[i]] = genes_mienturnet[2:nrow(genes_mienturnet),]
  rm(genes_mienturnet)
}

names(results_mienturnet) = names(Circos_universe)

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

for (i in 1:4){
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
  
  # Join with results_mienturnet to add miRNA degree to the genes
  BC = BC %>% left_join(results_mienturnet[[i]], by = "Gene.Symbol")
  zero_degree_index = which(is.na(BC$Degree) == TRUE)
  BC$Degree[zero_degree_index] = "0"; rm(zero_degree_index)
  str_sub(BC$Degree, 0, 0) = "degree="
  BC$GVars_new = paste0(BC$GVars, ",", BC$Degree)
  BC = BC %>% dplyr::select(-GVars, -Degree) %>%
    dplyr::rename(GVars = GVars_new)
  
  KG = KG %>% left_join(results_mienturnet[[i]], by = "Gene.Symbol")
  zero_degree_index = which(is.na(KG$Degree) == TRUE)
  KG$Degree[zero_degree_index] = "0"; rm(zero_degree_index)
  str_sub(KG$Degree, 0, 0) = "degree="
  KG$GVars_new = paste0(KG$GVars, ",", KG$Degree)
  KG = KG %>% dplyr::select(-GVars, -Degree) %>%
    dplyr::rename(GVars = GVars_new)
  
  # Preparing mutation and driver glyph files
  # Mutations
  BC_mutation = BC %>% dplyr::select(-GVars) %>%
    left_join(COSMIC_annot, by = c("Gene.Symbol")) %>%
    dplyr::select(-EntrezGene.ID)
  zero_mut_index = which(is.na(BC_mutation$Ratio) == TRUE)
  BC_mutation$Ratio[zero_mut_index] = 0; rm(zero_mut_index)
  BC_mutation = BC_mutation %>%
    mutate(MVars = paste0("ratio=", Ratio, ",sampmut=", Mutated.samples,
                          ",samptest=", Samples.tested)) %>%
    dplyr::select(-Ratio, -Mutated.Samples, -Samples.tested)
  
  KG_mutation = KG %>% dplyr::select(-GVars) %>%
    left_join(COSMIC_annot, by = c("Gene.Symbol")) %>%
    dplyr::select(-EntrezGene.ID)
  zero_mut_index = which(is.na(BC_mutation$Ratio) == TRUE)
  BC_mutation$Ratio[zero_mut_index] = 0; rm(zero_mut_index)
  BC_mutation = BC_mutation %>%
    mutate(MVars = paste0("ratio=", Ratio, ",sampmut=", Mutated.samples,
                          ",samptest=", Samples.tested)) %>%
    dplyr::select(-Ratio, -Mutated.Samples, -Samples.tested)
  
  # Drivers
  BC_drivers = BC %>% dplyr::select(-GVars)
  gdriver_index = which(BC_drivers$Gene.Symbol %in% general_drivers)
  pdriver_index = which(BC_drivers$Gene.Symbol %in% PDAC_drivers)
  BC_drivers$gdriver = "no"; BC_drivers$pdriver = "no"
  BC_drivers$gdriver[gdriver_index] = "yes"; BC_drivers$pdriver[pdriver_index] = "yes"
  rm(gdriver_index, pdriver_index)
  BC_drivers = BC_drivers %>%
    mutate(DrVars = paste0("gdriver=", gdriver, ",pdriver=", pdriver)) %>%
    dplyr::select(-gdriver, -pdriver)
  
  KG_drivers = KG %>% dplyr::select(-GVars)
  gdriver_index = which(KG_drivers$Gene.Symbol %in% general_drivers)
  pdriver_index = which(KG_drivers$Gene.Symbol %in% PDAC_drivers)
  KG_drivers$gdriver = "no"; KG_drivers$pdriver = "no"
  KG_drivers$gdriver[gdriver_index] = "yes"; KG_drivers$pdriver[pdriver_index] = "yes"
  rm(gdriver_index, pdriver_index)
  KG_drivers = KG_drivers %>%
    mutate(DrVars = paste0("gdriver=", gdriver, ",pdriver=", pdriver)) %>%
    dplyr::select(-gdriver, -pdriver)
  
  # Create histogram files for the genes for -log10(adjpval), degree, logfc:
  BC_histogram_logfc = BC %>% separate(GVars)
  
}


bp_gene_ont = bp_gene_ont_pre %>% 
  group_by(BP_new) %>%
  dplyr::select(-BP) %>%
  dplyr::rename(GChr = BP_new) %>%
  unique()

bp_gene_ont = bp_gene_ont[order(bp_gene_ont$GChr, bp_gene_ont$Gene.Symbol),]
bp_gene_ont$GStart = NA
bp_gene_ont$GEnd = NA
bp_gene_ont$GVars = paste0("logfc=",bp_gene_ont$logFC,",adjpval=",bp_gene_ont$adj.P.Val,
                           ",fullname=",bp_gene_ont$Name)
bp_gene_ont = bp_gene_ont %>%
  dplyr::select(-logFC, -adj.P.Val, -Name) %>%
  unique()
bp_gene_ont$GStart[1] = 0
bp_gene_ont$GEnd[1] = 5

for(i in 2:nrow(bp_gene_ont)){
  if(bp_gene_ont$GChr[i] == bp_gene_ont$GChr[i-1]){
    bp_gene_ont$GStart[i] = bp_gene_ont$GStart[i-1] + 5
    bp_gene_ont$GEnd[i]   = bp_gene_ont$GEnd[i-1] + 5
  } else {
    bp_gene_ont$GStart[i] = 0
    bp_gene_ont$GEnd[i]   = 5
  }
}

bp_gene_ont$GChr = gsub(" ", "_", bp_gene_ont$GChr)
str_sub(bp_gene_ont$GChr, 0, 0) = "hs"
bp_gene_ont = bp_gene_ont[,c("GChr", "GStart", "GEnd", "Gene.Symbol", "GVars")]
rm(bp_gene_ont_pre)


# nonas_DE_mapped
significants = nonas_DE_mapped[nonas_DE_mapped$adj.P.Val < 0.05,]
gene_list = significants$EntrezGene.ID
names(gene_list) = significants$EntrezGene.ID
gene_list = gene_list[!is.na(gene_list)]

RNGversion("4.0.2")
set.seed(123)
MF_DE_GO = groupGO(gene = gene_list,
                   OrgDb    = org.Hs.eg.db,
                   ont      = "MF",
                   level    = 2,
                   readable = FALSE)

newcols = c()
for (i in 1:max(MF_DE_GO@result$Count)){
  newcols = c(newcols, paste0("geneID_", i))
}
MF_DE_frame = MF_DE_GO@result %>% separate(geneID, into = newcols,
                                           sep = "/") %>%
  dplyr::select(-Count, -GeneRatio, -ID) %>%
  dplyr::filter(!geneID_1=="") %>%
  melt(id.vars = "Description") %>%
  dplyr::select(-variable) %>%
  dplyr::rename(EntrezGene.ID = value, MF = Description) %>%
  distinct()

RNGversion("4.0.2")
set.seed(123)
BP_DE_GO = groupGO(gene     = gene_list,
                   OrgDb    = org.Hs.eg.db,
                   ont      = "BP",
                   level    = 2,
                   readable = FALSE)

newcols = c()
for (i in 1:max(BP_DE_GO@result$Count)){
  newcols = c(newcols, paste0("geneID_", i))
}
BP_DE_frame = BP_DE_GO@result %>% separate(geneID, into = newcols,
                                           sep = "/") %>%
  dplyr::select(-Count, -GeneRatio, -ID) %>%
  dplyr::filter(!geneID_1=="") %>%
  melt(id.vars = "Description") %>%
  dplyr::select(-variable) %>%
  dplyr::rename(EntrezGene.ID = value, BP = Description) %>%
  distinct()

# nonas_z_DE_mapped 
significants = nonas_z_DE_mapped[nonas_z_DE_mapped$adj.P.Val < 0.05,]
gene_list = significants$EntrezGene.ID
names(gene_list) = significants$EntrezGene.ID
gene_list = gene_list[!is.na(gene_list)]

RNGversion("4.0.2")
set.seed(123)
z_MF_DE_GO = groupGO(gene = gene_list,
                     OrgDb    = org.Hs.eg.db,
                     ont      = "MF",
                     level    = 2,
                     readable = FALSE)

newcols = c()
for (i in 1:max(z_MF_DE_GO@result$Count)){
  newcols = c(newcols, paste0("geneID_", i))
}
z_MF_DE_frame = z_MF_DE_GO@result %>% separate(geneID, into = newcols,
                                               sep = "/") %>%
  dplyr::select(-Count, -GeneRatio, -ID) %>%
  dplyr::filter(!geneID_1=="") %>%
  melt(id.vars = "Description") %>%
  dplyr::select(-variable) %>%
  dplyr::rename(EntrezGene.ID = value, MF = Description) %>%
  distinct()

RNGversion("4.0.2")
set.seed(123)
z_BP_DE_GO = groupGO(gene     = gene_list,
                     OrgDb    = org.Hs.eg.db,
                     ont      = "BP",
                     level    = 2,
                     readable = FALSE)

newcols = c()
for (i in 1:max(z_BP_DE_GO@result$Count)){
  newcols = c(newcols, paste0("geneID_", i))
}
z_BP_DE_frame = z_BP_DE_GO@result %>% separate(geneID, into = newcols,
                                               sep = "/") %>%
  dplyr::select(-Count, -GeneRatio, -ID) %>%
  dplyr::filter(!geneID_1=="") %>%
  melt(id.vars = "Description") %>%
  dplyr::select(-variable) %>%
  dplyr::rename(EntrezGene.ID = value, BP = Description) %>%
  distinct()

# These now need to be matched to the two gene expression matrices
# We are gonna do a full match and then filter based on our interest for Circos
enr_DE_mapped = nonas_DE_mapped %>%
  left_join(BP_DE_frame, by = "EntrezGene.ID") %>%
  left_join(MF_DE_frame, by = "EntrezGene.ID")

enr_z_DE_mapped = nonas_z_DE_mapped %>%
  left_join(BP_DE_frame, by = "EntrezGene.ID") %>%
  left_join(MF_DE_frame, by = "EntrezGene.ID")

##### Matching with drug ontologies - keeping only genes with adj.p-value < 0.05 #####

DE_list = list(enr_DE_mapped, enr_z_DE_mapped)
names(DE_list) = c("noNAs Limma", "noNAs KBZ")
dirs = list("noNA_DGEA/Simple_Limma/nonas_",
            "noNA_DGEA/KBZ/nonas_z_")
final_pharmacologicals = list()
final_all_list = list()

for(i in 1:length(DE_list)){
  DE_filt = DE_list[[i]] %>% dplyr::filter(adj.P.Val < 0.05) %>% 
    dplyr::rename(DE_Symbol = Gene.Symbol)
  final_all_list[[i]] = ultimate_all %>% inner_join(DE_filt, by = "EntrezGene.ID")
  final_pharmacologicals[[i]] = ultimate_pharmacological %>% 
    inner_join(DE_filt, by = "EntrezGene.ID")
  write.xlsx(final_pharmacologicals[[i]], paste0(dirs[[i]], "Drug-Gene_Interactions.xlsx"))
  write.xlsx(final_all_list[[i]], paste0(dirs[[i]], "all_Drug-Gene_Interactions.xlsx"))
}

names(final_pharmacologicals) = names(DE_list)
names(final_all_list) = names(DE_list)

##### Generating links for Circos. Only the z-score-transformed DEGs #####
# All the files we are going to generate will depend on the gene ontology we pick
# There are two such cases, so the file preparation will happen twice:
# 1) for BP gene ontology and 2) for MF gene ontology

##### BP #####
bp_circos = list()

# master data frame
bp_circos_pre = final_pharmacologicals[["noNAs KBZ"]] %>%
  dplyr::select(-MF, -Name0, -Name2, -Name3, -DE_Symbol, -Sub4, -NewName) %>%
  unique()

# Genes
bp_gene_ont_pre = bp_circos_pre %>%
  dplyr::select(Gene.Symbol, BP, logFC, adj.P.Val, Name) %>%
  group_by(BP) %>%
  unique()
bp_gene_ont_pre$Name = gsub(" ", "_", bp_gene_ont_pre$Name)
bp_gene_ont_pre$Name = gsub("-", "_", bp_gene_ont_pre$Name)
bp_gene_ont_pre$Name = gsub(",", "_", bp_gene_ont_pre$Name)
bp_gene_ont_pre$Name = gsub("'", "_", bp_gene_ont_pre$Name)
bp_gene_ont_pre = bp_gene_ont_pre[order(bp_gene_ont_pre$BP, bp_gene_ont_pre$Gene.Symbol),]
bp_gene_ont_pre$BP[which(is.na(bp_gene_ont_pre$BP)==TRUE)] = "Other"
bp_gene_ont_pre$BP_new = factor(bp_gene_ont_pre$BP, levels = unique(bp_gene_ont_pre$BP),
                                labels = c("Behavior", "Adhesion", "Interspecies",
                                           "Bioregulation", "Biomineralization", "Cellular_process",
                                           "Detoxification", "Developmental_process", "Growth",
                                           "Immune_system", "Localization", "Locomotion", "Metabolism",
                                           "Multi_organism", "Multicellular", "BP_downregulation",
                                           "BP_upregulation", "Bioregulation", "Reproduction", 
                                           "Reproduction", "Stimuli_response", "Rhythmic_process",
                                           "Signaling", "Other")) # two bioregulation-related and reproduction-related ontologies were merged

bp_gene_ont = bp_gene_ont_pre %>% 
  group_by(BP_new) %>%
  dplyr::select(-BP) %>%
  dplyr::rename(GChr = BP_new) %>%
  unique()

bp_gene_ont = bp_gene_ont[order(bp_gene_ont$GChr, bp_gene_ont$Gene.Symbol),]
bp_gene_ont$GStart = NA
bp_gene_ont$GEnd = NA
bp_gene_ont$GVars = paste0("logfc=",bp_gene_ont$logFC,",adjpval=",bp_gene_ont$adj.P.Val,
                           ",fullname=",bp_gene_ont$Name)
bp_gene_ont = bp_gene_ont %>%
  dplyr::select(-logFC, -adj.P.Val, -Name) %>%
  unique()
bp_gene_ont$GStart[1] = 0
bp_gene_ont$GEnd[1] = 5

for(i in 2:nrow(bp_gene_ont)){
  if(bp_gene_ont$GChr[i] == bp_gene_ont$GChr[i-1]){
    bp_gene_ont$GStart[i] = bp_gene_ont$GStart[i-1] + 5
    bp_gene_ont$GEnd[i]   = bp_gene_ont$GEnd[i-1] + 5
  } else {
    bp_gene_ont$GStart[i] = 0
    bp_gene_ont$GEnd[i]   = 5
  }
}

bp_gene_ont$GChr = gsub(" ", "_", bp_gene_ont$GChr)
str_sub(bp_gene_ont$GChr, 0, 0) = "hs"
bp_gene_ont = bp_gene_ont[,c("GChr", "GStart", "GEnd", "Gene.Symbol", "GVars")]
rm(bp_gene_ont_pre)

# Drugs
bp_drug_ont_pre = bp_circos_pre %>%
  dplyr::select(API, Name1) %>%
  dplyr::rename(Class = Name1) %>%
  unique()
bp_drug_ont_pre = bp_drug_ont_pre[order(bp_drug_ont_pre$Class, bp_drug_ont_pre$API),]
bp_drug_ont_pre$Class_new = factor(bp_drug_ont_pre$Class, levels = unique(bp_drug_ont_pre$Class),
                                   labels = c("Misc", "Analgesics", "Anesthetics", "Misc",
                                              "Anti_acne", "Antianemics", "Dermatologicals",
                                              "Gastrointestinals", "Antiepileptics", "Dermatologicals", "Antineoplastics",
                                              "Antihemorrhagics", "Hormonal_etc", "Cardiovascular",
                                              "Antiinflammatory", "Antineoplastics", "Misc",
                                              "Antipsoriatics", "Misc", "Antithrombotics", "Misc",
                                              "Cardiovascular", "Hormonal_etc", "Cardiovascular",
                                              "Corticosteroids", "Dermatologicals", "Radiopharmaceuticals",
                                              "Gastrointestinals", "Cardiovascular", "Gastrointestinals", 
                                              "Gastrointestinals", "Obstructive_airway", "Musculoskeletals",
                                              "Antidiabetics", "Hormonal_etc", "Gynecologicals",
                                              "Immunostimulants", "Immunosuppressants", "Lipid_modifying",
                                              "Musculoskeletals", "Nasal_preps", "Opth_and_Oto", 
                                              "Opth_and_Oto", "Gastrointestinals", "Dermatologicals",
                                              "Musculoskeletals", "Gynecologicals", "Opth_and_Oto",
                                              "Cardiovascular", "Dermatologicals", "Misc",
                                              "Hormonal_etc", "Stomatological_preps", "Throat_preps",
                                              "Hormonal_etc", "Musculoskeletals", "Misc",
                                              "Urologicals", "Cardiovascular", "Vitamins")) 

# Every category that had dermatological drugs was put under "Dermatologicals",
# Colchicine was put under "Antineoplastics" instead of "Antigout"
# All musculoskeletals were grouped together, same for all gastrointestinals
# Decalinium was put under "Other" instead of antiseptics
# Cardiac therapy, vasodilators, antihypertensives and vasoprotectives were grouped
# together under Cardiovascular
# Glycyrrhizic acid was put under "Other" instead of Bile and Liver
# Tamoxifen, Mifepristone, Ulipristal, Levothyroxine, Rupatadine  were put under "Hormonal_etc"
# instead of endocrine therapy, sex hormones, sex hormones, thyroid and antihistamines respectively
# Anthelminthics and antiprotozoals were put under "Other"
# Caffeine was put under "Other" instead of psychoanaleptics


bp_drug_ont = bp_drug_ont_pre %>% 
  group_by(Class_new) %>%
  dplyr::rename(DChr = Class_new) %>%
  unique()
bp_drug_ont$DStart = NA
bp_drug_ont$DEnd = NA
bp_drug_ont$DVars = paste0("realname=",bp_drug_ont$Class)
bp_drug_ont = bp_drug_ont %>%
  dplyr::select(-Class) %>%
  unique() %>%
  group_by(DChr, API)
bp_drug_ont = bp_drug_ont[order(bp_drug_ont$DChr, bp_drug_ont$API),]
bp_drug_ont$DStart[1] = 0
bp_drug_ont$DEnd[1] = 5

for(i in 2:nrow(bp_drug_ont)){
  if(bp_drug_ont$DChr[i] == bp_drug_ont$DChr[i-1]){
    bp_drug_ont$DStart[i] = bp_drug_ont$DStart[i-1] + 5
    bp_drug_ont$DEnd[i]   = bp_drug_ont$DEnd[i-1] + 5
  } else {
    bp_drug_ont$DStart[i] = 0
    bp_drug_ont$DEnd[i]   = 5
  }
}

bp_drug_ont = bp_drug_ont[order(bp_drug_ont$DChr, bp_drug_ont$API),]
bp_drug_ont$DChr = gsub(" ", "_", bp_drug_ont$DChr)
bp_drug_ont$API = gsub(" ", "_", bp_drug_ont$API)
bp_drug_ont$API = gsub("-", "_", bp_drug_ont$API)
bp_drug_ont$API = str_to_title(bp_drug_ont$API)
bp_drug_ont$DVars = gsub(" ", "_", bp_drug_ont$DVars)
bp_drug_ont$DVars = gsub("-", "_", bp_drug_ont$DVars)
bp_drug_ont$DVars = gsub(",", "_", bp_drug_ont$DVars)
str_sub(bp_drug_ont$DChr, 0, 0) = "hs"
bp_drug_ont = bp_drug_ont[,c("DChr", "DStart", "DEnd", "API", "DVars")]
rm(bp_drug_ont_pre)
gc()

# Links
bp_circos_pre$API = str_to_title(bp_circos_pre$API)
bp_links = bp_circos_pre %>%
  dplyr::select(Gene.Symbol, API) %>%
  unique() %>%
  inner_join(bp_gene_ont, by = "Gene.Symbol") %>%
  inner_join(bp_drug_ont, by = "API") %>%
  mutate("LVars" = paste0(GVars,",",DVars,",gene=",Gene.Symbol,",drug=",API)) %>%
  dplyr::select(-GVars, -DVars, -Gene.Symbol, -API) %>%
  dplyr::select(GChr, GStart, GEnd, DChr, DStart, DEnd, LVars) %>%
  unique()

# logFC histogram
bp_logfc_hist = bp_gene_ont %>%
  separate(GVars, into = c("logFC", "adjpval", "fullname"), sep = ",") %>%
  dplyr::select(-Gene.Symbol) %>%
  mutate("logfc" = gsub("logfc=", "", logFC),
         "logfc_histvars" = paste0("neglogadjpval=",
                                   -log10(as.numeric(gsub("adjpval=", "", adjpval))), 
                                   ",", fullname)) %>%
  dplyr::select(-adjpval, -fullname, -logFC) %>%
  unique()

# -log10(adjpval) histogram
bp_adjpval_hist = bp_gene_ont %>%
  separate(GVars, into = c("logFC", "adjpval", "fullname"), sep = ",") %>%
  dplyr::select(-Gene.Symbol) %>%
  mutate("neg_log_adj_p" = -log10(as.numeric(gsub("adjpval=", "", adjpval))),
         "adjpval_histvars" = paste0(logFC,",", fullname)) %>%
  dplyr::select(-logFC, -fullname, -adjpval) %>%
  unique()

# Karyotype
bp_karyotype = data.frame(matrix(ncol = 6, 0))
for(i in unique(bp_drug_ont$DChr)){
  one = "chr"; two = "-"; three = i; four = gsub("hs", "", i); five = "0";
  six = max(bp_drug_ont$DEnd[bp_drug_ont$DChr == i])
  bp_karyotype = rbind(bp_karyotype, c(one, two, three, four, five, six))
  rm(one, two, three, four, five, six)
}
for(i in unique(bp_gene_ont$GChr)){
  one = "chr"; two = "-"; three = i; four = gsub("hs", "", i); five = "0";
  six = max(bp_gene_ont$GEnd[bp_gene_ont$GChr == i])
  bp_karyotype = rbind(bp_karyotype, c(one, two, three, four, five, six))
  rm(one, two, three, four, five, six)
}
bp_karyotype = bp_karyotype[2:nrow(bp_karyotype),]
bp_karyotype$X7 = "vvdblue"

# Final list
bp_circos[[1]] = bp_karyotype
bp_circos[[2]] = bp_drug_ont
bp_circos[[3]] = bp_gene_ont
bp_circos[[4]] = bp_links
bp_circos[[5]] = bp_logfc_hist
bp_circos[[6]] = bp_adjpval_hist
names(bp_circos) = c("Karyotype", "Drug_Labels", "Gene_Labels", "Links",
                     "logFC_histogram", "adjpval_histogram")

rm(bp_karyotype, bp_drug_ont, bp_gene_ont, bp_links, bp_logfc_hist, 
   bp_adjpval_hist, bp_circos_pre)

# Write out the .txt files that are necessary for the plots
for(i in 1:length(bp_circos)){
  write.table(bp_circos[[i]], paste0("Circos/Circos_BP/", names(bp_circos)[[i]], ".txt"), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

gc()

##### MF #####
mf_circos = list()

# master data frame
mf_circos_pre = final_pharmacologicals[["noNAs KBZ"]] %>%
  dplyr::select(-BP, -Name0, -Name2, -Name3, -DE_Symbol, -Sub4, -NewName) %>%
  unique()

# Genes
mf_gene_ont_pre = mf_circos_pre %>%
  dplyr::select(Gene.Symbol, MF, logFC, adj.P.Val, Name) %>%
  group_by(MF) %>%
  unique()
mf_gene_ont_pre$Name = gsub(" ", "_", mf_gene_ont_pre$Name)
mf_gene_ont_pre$Name = gsub("-", "_", mf_gene_ont_pre$Name)
mf_gene_ont_pre$Name = gsub(",", "_", mf_gene_ont_pre$Name)
mf_gene_ont_pre$Name = gsub("'", "_", mf_gene_ont_pre$Name)
mf_gene_ont_pre = mf_gene_ont_pre[order(mf_gene_ont_pre$MF, mf_gene_ont_pre$Gene.Symbol),]
mf_gene_ont_pre$MF[which(is.na(mf_gene_ont_pre$MF)==TRUE)] = "Other"
mf_gene_ont_pre$MF_new = factor(mf_gene_ont_pre$MF, levels = unique(mf_gene_ont_pre$MF),
                                labels = c("Antioxidants", "Binding", "Cargo_receptor",
                                           "Catalysts", "Adaptors", "Carriers", 
                                           "MF_regulators", "Transducers", 
                                           "Structural_molecules", "Transcription_regulation",
                                           "Transporters", "Other"))

mf_gene_ont = mf_gene_ont_pre %>% 
  group_by(MF_new) %>%
  dplyr::select(-MF) %>%
  dplyr::rename(GChr = MF_new) %>%
  unique()

mf_gene_ont = mf_gene_ont[order(mf_gene_ont$GChr, mf_gene_ont$Gene.Symbol),]
mf_gene_ont$GStart = NA
mf_gene_ont$GEnd = NA
mf_gene_ont$GVars = paste0("logfc=",mf_gene_ont$logFC,",adjpval=",mf_gene_ont$adj.P.Val,
                           ",fullname=",mf_gene_ont$Name)
mf_gene_ont = mf_gene_ont %>%
  dplyr::select(-logFC, -adj.P.Val, -Name) %>%
  unique()
mf_gene_ont$GStart[1] = 0
mf_gene_ont$GEnd[1] = 5

for(i in 2:nrow(mf_gene_ont)){
  if(mf_gene_ont$GChr[i] == mf_gene_ont$GChr[i-1]){
    mf_gene_ont$GStart[i] = mf_gene_ont$GStart[i-1] + 5
    mf_gene_ont$GEnd[i]   = mf_gene_ont$GEnd[i-1] + 5
  } else {
    mf_gene_ont$GStart[i] = 0
    mf_gene_ont$GEnd[i]   = 5
  }
}

mf_gene_ont$GChr = gsub(" ", "_", mf_gene_ont$GChr)
str_sub(mf_gene_ont$GChr, 0, 0) = "hs"
mf_gene_ont = mf_gene_ont[,c("GChr", "GStart", "GEnd", "Gene.Symbol", "GVars")]
rm(mf_gene_ont_pre)

# Drugs
mf_drug_ont_pre = mf_circos_pre %>%
  dplyr::select(API, Name1) %>%
  dplyr::rename(Class = Name1) %>%
  unique()
mf_drug_ont_pre = mf_drug_ont_pre[order(mf_drug_ont_pre$Class, mf_drug_ont_pre$API),]
mf_drug_ont_pre$Class_new = factor(mf_drug_ont_pre$Class, levels = unique(mf_drug_ont_pre$Class),
                                   labels = c("Misc", "Analgesics", "Anesthetics", "Misc",
                                              "Anti_acne", "Antianemics", "Dermatologicals",
                                              "Gastrointestinals", "Antiepileptics", "Dermatologicals", "Antineoplastics",
                                              "Antihemorrhagics", "Hormonal_etc", "Cardiovascular",
                                              "Antiinflammatory", "Antineoplastics", "Misc",
                                              "Antipsoriatics", "Misc", "Antithrombotics", "Misc",
                                              "Cardiovascular", "Hormonal_etc", "Cardiovascular",
                                              "Corticosteroids", "Dermatologicals", "Radiopharmaceuticals",
                                              "Gastrointestinals", "Cardiovascular", "Gastrointestinals", 
                                              "Gastrointestinals", "Obstructive_airway", "Musculoskeletals",
                                              "Antidiabetics", "Hormonal_etc", "Gynecologicals",
                                              "Immunostimulants", "Immunosuppressants", "Lipid_modifying",
                                              "Musculoskeletals", "Nasal_preps", "Opth_and_Oto", 
                                              "Opth_and_Oto", "Gastrointestinals", "Dermatologicals",
                                              "Musculoskeletals", "Gynecologicals", "Opth_and_Oto",
                                              "Cardiovascular", "Dermatologicals", "Misc",
                                              "Hormonal_etc", "Stomatological_preps", "Throat_preps",
                                              "Hormonal_etc", "Musculoskeletals", "Misc",
                                              "Urologicals", "Cardiovascular", "Vitamins")) 

# Every category that had dermatological drugs was put under "Dermatologicals",
# Colchicine was put under "Antineoplastics" instead of "Antigout"
# All musculoskeletals were grouped together, same for all gastrointestinals
# Decalinium was put under "Other" instead of antiseptics
# Cardiac therapy, vasodilators, antihypertensives and vasoprotectives were grouped
# together under Cardiovascular
# Glycyrrhizic acid was put under "Other" instead of Bile and Liver
# Tamoxifen, Mifepristone, Ulipristal, Levothyroxine, Rupatadine  were put under "Hormonal_etc"
# instead of endocrine therapy, sex hormones, sex hormones, thyroid and antihistamines respectively
# Anthelminthics and antiprotozoals were put under "Other"
# Caffeine was put under "Other" instead of psychoanaleptics


mf_drug_ont = mf_drug_ont_pre %>% 
  group_by(Class_new) %>%
  dplyr::rename(DChr = Class_new) %>%
  unique()
mf_drug_ont$DStart = NA
mf_drug_ont$DEnd = NA
mf_drug_ont$DVars = paste0("realname=",mf_drug_ont$Class)
mf_drug_ont = mf_drug_ont %>%
  dplyr::select(-Class) %>%
  unique() %>%
  group_by(DChr, API)
mf_drug_ont = mf_drug_ont[order(mf_drug_ont$DChr, mf_drug_ont$API),]
mf_drug_ont$DStart[1] = 0
mf_drug_ont$DEnd[1] = 5

for(i in 2:nrow(mf_drug_ont)){
  if(mf_drug_ont$DChr[i] == mf_drug_ont$DChr[i-1]){
    mf_drug_ont$DStart[i] = mf_drug_ont$DStart[i-1] + 5
    mf_drug_ont$DEnd[i]   = mf_drug_ont$DEnd[i-1] + 5
  } else {
    mf_drug_ont$DStart[i] = 0
    mf_drug_ont$DEnd[i]   = 5
  }
}

mf_drug_ont = mf_drug_ont[order(mf_drug_ont$DChr, mf_drug_ont$API),]
mf_drug_ont$DChr = gsub(" ", "_", mf_drug_ont$DChr)
mf_drug_ont$API = gsub(" ", "_", mf_drug_ont$API)
mf_drug_ont$API = gsub("-", "_", mf_drug_ont$API)
mf_drug_ont$API = str_to_title(mf_drug_ont$API)
mf_drug_ont$DVars = gsub(" ", "_", mf_drug_ont$DVars)
mf_drug_ont$DVars = gsub("-", "_", mf_drug_ont$DVars)
mf_drug_ont$DVars = gsub(",", "_", mf_drug_ont$DVars)
str_sub(mf_drug_ont$DChr, 0, 0) = "hs"
mf_drug_ont = mf_drug_ont[,c("DChr", "DStart", "DEnd", "API", "DVars")]
rm(mf_drug_ont_pre)
gc()

# Links
mf_circos_pre$API = str_to_title(mf_circos_pre$API)
mf_links = mf_circos_pre %>%
  dplyr::select(Gene.Symbol, API) %>%
  unique() %>%
  inner_join(mf_gene_ont, by = "Gene.Symbol") %>%
  inner_join(mf_drug_ont, by = "API") %>%
  mutate("LVars" = paste0(GVars,",",DVars,",gene=",Gene.Symbol,",drug=",API)) %>%
  dplyr::select(-GVars, -DVars, -Gene.Symbol, -API) %>%
  dplyr::select(GChr, GStart, GEnd, DChr, DStart, DEnd, LVars) %>%
  unique()

# logFC histogram
mf_logfc_hist = mf_gene_ont %>%
  separate(GVars, into = c("logFC", "adjpval", "fullname"), sep = ",") %>%
  dplyr::select(-Gene.Symbol) %>%
  mutate("logfc" = gsub("logfc=", "", logFC),
         "logfc_histvars" = paste0("neglogadjpval=",
                                   -log10(as.numeric(gsub("adjpval=", "", adjpval))), 
                                   ",", fullname)) %>%
  dplyr::select(-adjpval, -fullname, -logFC) %>%
  unique()

# -log10(adjpval) histogram
mf_adjpval_hist = mf_gene_ont %>%
  separate(GVars, into = c("logFC", "adjpval", "fullname"), sep = ",") %>%
  dplyr::select(-Gene.Symbol) %>%
  mutate("neg_log_adj_p" = -log10(as.numeric(gsub("adjpval=", "", adjpval))),
         "adjpval_histvars" = paste0(logFC,",", fullname)) %>%
  dplyr::select(-logFC, -fullname, -adjpval) %>%
  unique()

# Karyotype
mf_karyotype = data.frame(matrix(ncol = 6, 0))
for(i in unique(mf_drug_ont$DChr)){
  one = "chr"; two = "-"; three = i; four = gsub("hs", "", i); five = "0";
  six = max(mf_drug_ont$DEnd[mf_drug_ont$DChr == i])
  mf_karyotype = rbind(mf_karyotype, c(one, two, three, four, five, six))
  rm(one, two, three, four, five, six)
}
for(i in unique(mf_gene_ont$GChr)){
  one = "chr"; two = "-"; three = i; four = gsub("hs", "", i); five = "0";
  six = max(mf_gene_ont$GEnd[mf_gene_ont$GChr == i])
  mf_karyotype = rbind(mf_karyotype, c(one, two, three, four, five, six))
  rm(one, two, three, four, five, six)
}
mf_karyotype = mf_karyotype[2:nrow(mf_karyotype),]
mf_karyotype$X7 = "vvdblue"

# Final list
mf_circos[[1]] = mf_karyotype
mf_circos[[2]] = mf_drug_ont
mf_circos[[3]] = mf_gene_ont
mf_circos[[4]] = mf_links
mf_circos[[5]] = mf_logfc_hist
mf_circos[[6]] = mf_adjpval_hist
names(mf_circos) = c("Karyotype", "Drug_Labels", "Gene_Labels", "Links",
                     "logFC_histogram", "adjpval_histogram")

rm(mf_karyotype, mf_drug_ont, mf_gene_ont, mf_links, mf_logfc_hist, 
   mf_adjpval_hist, mf_circos_pre)

# Write out the .txt files that are necessary for the plots
for(i in 1:length(mf_circos)){
  write.table(mf_circos[[i]], paste0("Circos/Circos_MF/", names(mf_circos)[[i]], ".txt"), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

gc()

##### Further processing #####
# Manually filtering for approved drugs for pancreatic cancer
PC_drugs = c("Paclitaxel", "Irinotecan", "Cisplatin", "Sunitinib")
BP_PC_drug_interactions_all = final_all_list[["noNAs KBZ"]] %>% 
  dplyr::filter(API %in% tolower(PC_drugs)) %>%
  dplyr::select(-MF, -NewName, -Name0, -Name2, -Name3, -EntrezGene.ID, -Sub4, -DE_Symbol) %>%
  dplyr::select(Name, BP, Gene.Symbol, API, everything()) %>%
  unique()
BP_PC_drug_interactions_pharmacological = final_pharmacologicals[["noNAs KBZ"]] %>% 
  dplyr::filter(API %in% tolower(PC_drugs)) %>%
  dplyr::select(-MF, -NewName, -Name0, -Name2, -Name3, -EntrezGene.ID, -Sub4, -DE_Symbol) %>%
  dplyr::select(Name, BP, Gene.Symbol, API, everything()) %>%
  unique()
MF_PC_drug_interactions_all = final_all_list[["noNAs KBZ"]] %>% 
  dplyr::filter(API %in% tolower(PC_drugs)) %>%
  dplyr::select(-BP, -NewName, -Name0, -Name2, -Name3, -EntrezGene.ID, -Sub4, -DE_Symbol) %>%
  dplyr::select(Name, MF, Gene.Symbol, API, everything()) %>%
  unique()
MF_PC_drug_interactions_pharmacological = final_pharmacologicals[["noNAs KBZ"]] %>% 
  dplyr::filter(API %in% tolower(PC_drugs)) %>%
  dplyr::select(-BP, -NewName, -Name0, -Name2, -Name3, -EntrezGene.ID, -Sub4, -DE_Symbol) %>%
  dplyr::select(Name, MF, Gene.Symbol, API, everything()) %>%
  unique()

# Write out interactions as an MS Excel file with multiple sheets
PC_drug_interactions = createWorkbook()
addWorksheet(PC_drug_interactions, "BP_all")
writeData(PC_drug_interactions, "BP_all", BP_PC_drug_interactions_all)
addWorksheet(PC_drug_interactions, "BP_pharm")
writeData(PC_drug_interactions, "BP_pharm", BP_PC_drug_interactions_pharmacological)
addWorksheet(PC_drug_interactions, "MF_all")
writeData(PC_drug_interactions, "MF_all", MF_PC_drug_interactions_all)
addWorksheet(PC_drug_interactions, "MF_pharm")
writeData(PC_drug_interactions, "MF_pharm", MF_PC_drug_interactions_pharmacological)
saveWorkbook(PC_drug_interactions, "PC_drug_interactions.xlsx", overwrite = TRUE)
