# Pathway and Gene Ontology analysis using pathfindR, clusterProfiler

# Libraries #####
library(pathfindR)
library(openxlsx)
library(dplyr)
library(ggplot2)

# Visualisations of pathways have been moved to a folder outside the repository
# due to large size and difficulties with uploading to GitHub

# Loading the input to pathfindR (the stage 1 vs normal topTable output):
Stage_1_tT = read.xlsx("DGEA/Union/Stage_1_vs_Normal_z_DE_topTable.xlsx") %>%
  dplyr::select(Gene.Symbol, logFC, adj.P.Val) %>%
  na.omit()

# Preparing a time-expensive pathfindR loop for enrichment analysis
gene_sets = c("BioCarta", "GO-All", "KEGG", "Reactome")
dirs = paste0("pathfindR/", gene_sets)
pathfindR_outputs = list()

for (i in 1:length(dirs)){
  pathfindR_outputs[[i]] = run_pathfindR(Stage_1_tT, gene_sets = gene_sets[i],
                                   p_val_threshold = 0.05, 
                                   visualize_enriched_terms = TRUE,
                                   output_dir = dirs[i], min_gset_size = 10,
                                   max_gset_size = 300, adj_method = 'fdr',
                                   enrichment_threshold = 0.05,
                                   pin_name_path = 'Biogrid', search_method = 'GR',
                                   grMaxDepth = 1, grSearchDepth = 1,
                                   iterations = 10, n_processes = 10)
  cat(paste0("Done with ", gene_sets[i], "\n"))
}
names(pathfindR_outputs) = gene_sets

# Perform hierarchical clustering on the results (average distance metric)
# Not run for GO-All and Reactome because the algorithm time complexity is O(n^3)

clustered_results = list()
for (i in c(1,3)){
  clustered_results[[i]] = cluster_enriched_terms(pathfindR_outputs[[i]],
                                                  method = "hierarchical")
}
names(clustered_results) = c("BioCarta", "NULL", "KEGG")

# BioCarta: The maximum average silhouette width was 0.17 for k = 66
# KEGG    : The maximum average silhouette width was 0.14 for k = 66

enrichment_dotplots = list()
cluster_enrichment_dotplots = list()

# Producing dotplots with the results
for (i in 1:length(pathfindR_outputs)){
  # unclustered results
  enrichment_dotplots[[i]] = enrichment_chart(result_df = pathfindR_outputs[[i]], top_terms = 10)
  png(paste0("pathfindR/", names(pathfindR_outputs)[i], "/",
             names(pathfindR_outputs)[i], "_top10_dotplot.png"), width = 1024, height = 768)
  print(enrichment_dotplots[[i]])
  dev.off()
  
  if (i == 1 | i == 3){
    # clustered results
    cluster_enrichment_dotplots[[i]] = enrichment_chart(result_df = clustered_results[[i]], top_terms = 10,
                            plot_by_cluster = TRUE)
    png(paste0("pathfindR/", names(pathfindR_outputs)[i], "/",
               names(pathfindR_outputs)[i], "_top10_dotplot_clustered.png"), width = 1024, 
        height = 768)
    print(cluster_enrichment_dotplots[[i]])
    dev.off()
  }
}

names(enrichment_dotplots) = names(pathfindR_outputs)
names(cluster_enrichment_dotplots) = names(clustered_results)

# Write out results in a comprehensive .xlsx file
wb = createWorkbook()
addWorksheet(wb, "BioCarta")
writeData(wb, "BioCarta", clustered_results[["BioCarta"]])
addWorksheet(wb, "GO-All")
writeData(wb, "GO-All", pathfindR_outputs[["GO-All"]])
addWorksheet(wb, "KEGG")
writeData(wb, "KEGG", clustered_results[["KEGG"]])
addWorksheet(wb, "Reactome")
writeData(wb, "Reactome", pathfindR_outputs[["Reactome"]])
saveWorkbook(wb, file = "pathfindR/Comprehensive_pathfindR_output.xlsx",
             overwrite = TRUE); rm(wb)

# Representative terms output file (BioCarta & KEGG)
wb2 = createWorkbook()
addWorksheet(wb2, "BioCarta - rep")
writeData(wb2, "BioCarta - rep", clustered_results[["BioCarta"]][clustered_results[["BioCarta"]]$Status
                                                                 == "Representative", ])
addWorksheet(wb2, "KEGG - rep")
writeData(wb2, "KEGG - rep", clustered_results[["KEGG"]][clustered_results[["KEGG"]]$Status
                                                                 == "Representative", ])
saveWorkbook(wb2, file = "pathfindR/Representative_terms.xlsx",
             overwrite = TRUE); rm(wb2)

# Term-gene heatmaps and term-gene graphs #####
term_gene_heatmaps = list()
term_gene_graphs = list()

for (i in 1:length(pathfindR_outputs)){
  # term-gene heatmaps
  term_gene_heatmaps[[i]] = term_gene_heatmap(result_df = pathfindR_outputs[[i]],
                                              genes_df = Stage_1_tT,
                                              num_terms = 5,
                                              use_description = TRUE,
                                              low = "midnightblue",
                                              high = "hotpink4",
                                              mid = "white")
  png(paste0("pathfindR/", names(pathfindR_outputs)[i], "/",
             names(pathfindR_outputs)[i], "_top5_term_gene_heatmap.png"), 
      width = 1920, height = 1080)
  print(term_gene_heatmaps[[i]])
  dev.off()
  
  # term-gene graphs
  term_gene_graphs[[i]] = term_gene_graph(result_df = pathfindR_outputs[[i]],
                                              num_terms = 3,
                                              use_description = TRUE,
                                          node_size = "p_val")
  png(paste0("pathfindR/", names(pathfindR_outputs)[i], "/",
             names(pathfindR_outputs)[i], "_top3_term_gene_graph.png"), 
      width = 1024, height = 768)
  print(term_gene_graphs[[i]])
  dev.off()
}

# UpSet plots #####
UpSet_plots = list()
for (i in 1:length(pathfindR_outputs)){
  # term-gene heatmaps
  UpSet_plots[[i]] = UpSet_plot(result_df = pathfindR_outputs[[i]],
                                genes_df = Stage_1_tT,
                                num_terms = 5,
                                use_description = TRUE,
                                low = "midnightblue",
                                high = "hotpink4",
                                mid = "white")
  png(paste0("pathfindR/", names(pathfindR_outputs)[i], "/",
             names(pathfindR_outputs)[i], "_top5_UpSet_plot.png"), 
      width = 1920, height = 1080)
  print(UpSet_plots[[i]])
  dev.off()
}

names(term_gene_graphs) = names(pathfindR_outputs)
names(term_gene_heatmaps) = names(pathfindR_outputs)
names(UpSet_plots) = names(pathfindR_outputs)

##### REQUIRES COMPLETE-CASE MATRIX #####
# Calculate agglomerated z-scores for representative terms and plot heatmap 
# using term descriptions
pdata = read.xlsx("DGEA/Pheno.xlsx")
pdata$Stage1 = 0
pdata$Stage1[pdata$AJCC_classification == "1a"] = 1
pdata$Stage1[pdata$AJCC_classification == "1b"] = 1

exp_mat = read.xlsx("DGEA/Union/Stage_1_vs_Normal_exp_mat.xlsx")
names = exp_mat$EntrezGene.ID
exp_mat = as.matrix(exp_mat %>% dplyr::select(-EntrezGene.ID))
rownames(exp_mat) = names; rm(names)

stage1 = intersect(colnames(exp_mat), pdata$GEO_accession[pdata$Stage1 == 1])
normal = intersect(colnames(exp_mat), pdata$GEO_accession[pdata$Tissue_type == "non_tumor"])

sample_score_plots = list()
for (i in c(1,3)){
  sample_score_plots[[i]] = score_terms(enrichment_table = 
                                          clustered_results[[i]][clustered_results[[i]]$Status
                                                                                  == "Representative", ],
                             exp_mat = exp_mat,
                             cases = stage1,
                             use_description = TRUE, # default FALSE
                             label_samples = FALSE, # default = TRUE
                             case_title = "Stage 1 PDAC",  # default = "Case"
                             control_title = "Non-tumor samples", # default = "Control"
                             low = "midnightblue", # default = "green"
                             mid = "white", # default = "black"
                             high = "hotpink4")
  png(paste0("pathfindR/", names(pathfindR_outputs)[i], "/",
             names(pathfindR_outputs)[i], "_sample_scores.png"), 
      width = 1920, height = 1080)
  print(sample_score_plots[[i]])
  dev.off()
}
#####

# Driver genes from CIRCOS census list #####
census = read.csv("COSMIC_census_09_02_2022.csv") %>%
  dplyr::rename(EntrezGene.ID = Entrez.GeneId, Gene.Symbol_COSMIC = Gene.Symbol) %>%
  dplyr::select(EntrezGene.ID, Gene.Symbol_COSMIC, Name, everything())
census$EntrezGene.ID = as.character(census$EntrezGene.ID)

sig_DGEA_subset = read.xlsx("DGEA/Union/Stage_1_vs_Normal_z_DE_topTable.xlsx") %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::rename(Gene.Symbol_pre = Gene.Symbol)

general_drivers = inner_join(census, sig_DGEA_subset, by = "EntrezGene.ID") %>%
  dplyr::select(EntrezGene.ID, Gene.Symbol_pre, Gene.Symbol_COSMIC, Name, logFC, 
                P.Value, B, HGNC_Official, adj.P.Val, everything())
PDAC_drivers = inner_join(census, sig_DGEA_subset, by = "EntrezGene.ID") %>%
  dplyr::filter(grepl("pancr", Tumour.Types.Somatic.) | grepl("pancr", Tumour.Types.Germline.)) %>%
  dplyr::select(EntrezGene.ID, Gene.Symbol_pre, Gene.Symbol_COSMIC, Name, logFC, 
                P.Value, B, HGNC_Official, adj.P.Val, everything())

# Writing out
drivers = createWorkbook()
addWorksheet(drivers, "General")
writeData(drivers, "General", general_drivers)
addWorksheet(drivers, "Pancreas")
writeData(drivers, "Pancreas", PDAC_drivers)
saveWorkbook(drivers, "Drivers.xlsx", overwrite = TRUE); rm(drivers); gc()
