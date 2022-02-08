# Pathway and Gene Ontology analysis using pathfindR, clusterProfiler

# Libraries #####
library(pathfindR)
library(openxlsx)
library(dplyr)

# Visualisations of pathways have been moved to a folder outside the repository
# due to large size and difficulties with uploading to GitHub

# Loading the input to pathfindR (the stage 1 vs normal topTable output):
Stage_1_tT = read.xlsx("DGEA/Union/Stage_1_vs_Normal_z_DE_topTable.xlsx") %>%
  dplyr::select(Gene.Symbol, logFC, adj.P.Val) %>%
  na.omit()

# Preparing a time-expensive pathfindR loop
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

# Producing dotplots with the results
for (i in 1:length(pathfindR_outputs)){
  # unclustered results
  plot = enrichment_chart(result_df = pathfindR_outputs[[i]], top_terms = 10)
  png(paste0("pathfindR/", names(pathfindR_outputs)[i], "/",
             names(pathfindR_outputs)[i], "_top10_dotplot.png"), width = 1024, height = 768)
  print(plot)
  dev.off(); rm(plot)
  
  if (i == 1 | i == 3){
    # clustered results
    plot = enrichment_chart(result_df = clustered_results[[i]], top_terms = 10,
                            plot_by_cluster = TRUE)
    png(paste0("pathfindR/", names(pathfindR_outputs)[i], "/",
               names(pathfindR_outputs)[i], "_top10_dotplot_clustered.png"), width = 1024, 
        height = 768)
    print(plot)
    dev.off(); rm(plot)
  }
}

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
