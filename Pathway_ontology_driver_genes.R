# Pathway and Gene Ontology analysis using pathfindR, clusterProfiler

# Libraries #####
library(pathfindR)
library(openxlsx)
library(dplyr)
library(ggplot2)

# Visualisations of pathways have been moved to a folder outside the repository
# due to large size and difficulties with uploading to GitHub
gene_sets = c("BioCarta", "GO-BP", "GO-CC", "GO-MF", "KEGG", "Reactome")
axis_text_size = c(10, 5, 5, 5, 5, 10)
names(axis_text_size) = gene_sets

# Stage 1 #####
# Loading the input to pathfindR (the stage 1 vs normal topTable output):
Stage_1_tT = read.xlsx("DGEA/Union/Stage_1_vs_Normal_z_DE_topTable.xlsx") %>%
  dplyr::select(Gene.Symbol, logFC, adj.P.Val) %>%
  na.omit()

# Preparing a pathfindR loop for enrichment analysis
dirs_stage_1 = paste0("pathfindR/Stage_1/", gene_sets)
pathfindR_outputs_stage_1 = list()

RNGversion("4.0.2")
set.seed(123)
for (i in 1:length(dirs_stage_1)){
  pathfindR_outputs_stage_1[[i]] = run_pathfindR(Stage_1_tT, gene_sets = gene_sets[i],
                                   p_val_threshold = 0.05, 
                                   visualize_enriched_terms = TRUE,
                                   output_dir = dirs_stage_1[i], min_gset_size = 10,
                                   max_gset_size = 300, adj_method = 'fdr',
                                   enrichment_threshold = 0.05,
                                   pin_name_path = 'Biogrid', search_method = 'GR',
                                   grMaxDepth = 1, grSearchDepth = 1,
                                   iterations = 10, n_processes = 10)
  cat(paste0("Done with ", gene_sets[i], "\n"))
}
names(pathfindR_outputs_stage_1) = gene_sets

# Perform hierarchical clustering on the results (average distance metric)
# Not run for GO-BP and Reactome because the algorithm time complexity is O(n^3)

clustered_results_stage_1 = list()
for (i in c(1, 3, 4, 5)){
  clustered_results_stage_1[[i]] = cluster_enriched_terms(pathfindR_outputs_stage_1[[i]],
                                                  method = "hierarchical")
}
names(clustered_results_stage_1) = c("BioCarta", "NULL", "GO-CC", "GO-MF", "KEGG")

# BioCarta : The maximum average silhouette width was 0.18 for k = 65
# GO-CC    : The maximum average silhouette width was 0.12 for k = 139
# GO-MF    : The maximum average silhouette width was 0.11 for k = 172
# KEGG     : The maximum average silhouette width was 0.14 for k = 68

enrichment_dotplots_stage_1 = list()
cluster_enrichment_dotplots_stage_1 = list()

# Producing dotplots with the results
for (i in 1:length(pathfindR_outputs_stage_1)){
  # unclustered results
  enrichment_dotplots_stage_1[[i]] = enrichment_chart(result_df = pathfindR_outputs_stage_1[[i]],
                                                      top_terms = 10)
  tiff(paste0("pathfindR/Stage_1/", names(pathfindR_outputs_stage_1)[i], "/",
             names(pathfindR_outputs_stage_1)[i], "_top10_dotplot.tif"), 
      width = 1920, height = 1080, res = 150)
  print(enrichment_dotplots_stage_1[[i]])
  dev.off()
  
  if (i == 1 | i == 3 | i == 4 | i == 5){
    # clustered results
    cluster_enrichment_dotplots_stage_1[[i]] = enrichment_chart(result_df = clustered_results_stage_1[[i]][clustered_results_stage_1[[i]]$Status
                                                                                                                    == "Representative", ][1:10,],
                                                                top_terms = NULL,
                            plot_by_cluster = TRUE)
    tiff(paste0("pathfindR/Stage_1/", names(pathfindR_outputs_stage_1)[i], "/",
               names(pathfindR_outputs_stage_1)[i], "_top10_dotplot_clustered.tif"), 
         width = 1920, height = 1080, res = 150)
    print(cluster_enrichment_dotplots_stage_1[[i]])
    dev.off()
  }
}

names(enrichment_dotplots_stage_1) = names(pathfindR_outputs_stage_1)
names(cluster_enrichment_dotplots_stage_1) = names(clustered_results_stage_1)

# Write out results in a comprehensive .xlsx file
wb = createWorkbook()
addWorksheet(wb, "BioCarta")
writeData(wb, "BioCarta", clustered_results_stage_1[["BioCarta"]])
addWorksheet(wb, "GO-BP")
writeData(wb, "GO-BP", pathfindR_outputs_stage_1[["GO-BP"]])
addWorksheet(wb, "GO-CC")
writeData(wb, "GO-CC", pathfindR_outputs_stage_1[["GO-CC"]])
addWorksheet(wb, "GO-MF")
writeData(wb, "GO-MF", pathfindR_outputs_stage_1[["GO-MF"]])
addWorksheet(wb, "KEGG")
writeData(wb, "KEGG", clustered_results_stage_1[["KEGG"]])
addWorksheet(wb, "Reactome")
writeData(wb, "Reactome", pathfindR_outputs_stage_1[["Reactome"]])
saveWorkbook(wb, file = "pathfindR/Stage_1/Comprehensive_pathfindR_output.xlsx",
             overwrite = TRUE); rm(wb)

# Representative terms output file (BioCarta, GO-CC, GO-MF & KEGG)
wb2 = createWorkbook()
addWorksheet(wb2, "BioCarta - rep")
writeData(wb2, "BioCarta - rep", clustered_results_stage_1[["BioCarta"]][clustered_results_stage_1[["BioCarta"]]$Status
                                                                 == "Representative", ])
addWorksheet(wb2, "GOCC - rep")
writeData(wb2, "GOCC - rep", clustered_results_stage_1[["GO-CC"]][clustered_results_stage_1[["GO-CC"]]$Status
                                                                         == "Representative", ])
addWorksheet(wb2, "GOMF - rep")
writeData(wb2, "GOMF - rep", clustered_results_stage_1[["GO-MF"]][clustered_results_stage_1[["GO-MF"]]$Status
                                                                  == "Representative", ])
addWorksheet(wb2, "KEGG - rep")
writeData(wb2, "KEGG - rep", clustered_results_stage_1[["KEGG"]][clustered_results_stage_1[["KEGG"]]$Status
                                                                 == "Representative", ])
saveWorkbook(wb2, file = "pathfindR/Stage_1/Representative_terms.xlsx",
             overwrite = TRUE); rm(wb2)

# Term-gene heatmaps and term-gene graphs #####
term_gene_heatmaps_stage_1 = list()
term_gene_graphs_stage_1 = list()

for (i in 1:length(pathfindR_outputs_stage_1)){
  # term-gene heatmaps
  term_gene_heatmaps_stage_1[[i]] = term_gene_heatmap(result_df = pathfindR_outputs_stage_1[[i]],
                                              genes_df = Stage_1_tT,
                                              num_terms = 5,
                                              use_description = TRUE,
                                              low = "darkgreen",
                                              high = "darkred",
                                              mid = "black")+
    theme(axis.text.x = element_text(size = axis_text_size[[i]]),
          axis.text.y = element_text(size = 12))
  tiff(paste0("pathfindR/Stage_1/", names(pathfindR_outputs_stage_1)[i], "/",
             names(pathfindR_outputs_stage_1)[i], "_top5_term_gene_heatmap.tif"), 
      width = 3840, height = 648, res = 150)
  print(term_gene_heatmaps_stage_1[[i]])
  dev.off()
  
  # term-gene graphs
  term_gene_graphs_stage_1[[i]] = term_gene_graph(result_df = pathfindR_outputs_stage_1[[i]],
                                              num_terms = 3,
                                              use_description = TRUE,
                                          node_size = "p_val")
  tiff(paste0("pathfindR/Stage_1/", names(pathfindR_outputs_stage_1)[i], "/",
             names(pathfindR_outputs_stage_1)[i], "_top3_term_gene_graph.tif"), 
      width = 1920, height = 1080, res = 100)
  print(term_gene_graphs_stage_1[[i]])
  dev.off()
}

# UpSet plots #####
UpSet_plots_stage_1 = list()
for (i in 1:length(pathfindR_outputs_stage_1)){
  # term-gene heatmaps
  UpSet_plots_stage_1[[i]] = UpSet_plot(result_df = pathfindR_outputs_stage_1[[i]],
                                genes_df = Stage_1_tT,
                                num_terms = 5,
                                use_description = TRUE,
                                low = "darkgreen",
                                high = "darkred",
                                mid = "black")+
    theme(axis.text.y = element_text(size = axis_text_size[[i]]))
  tiff(paste0("pathfindR/Stage_1/", names(pathfindR_outputs_stage_1)[i], "/",
             names(pathfindR_outputs_stage_1)[i], "_top5_UpSet_plot.tif"), 
      width = 1444, height = 3840, res = 150)
  print(UpSet_plots_stage_1[[i]])
  dev.off()
}

names(term_gene_graphs_stage_1) = names(pathfindR_outputs_stage_1)
names(term_gene_heatmaps_stage_1) = names(pathfindR_outputs_stage_1)
names(UpSet_plots_stage_1) = names(pathfindR_outputs_stage_1)

# Stage 2 #####
# Loading the input to pathfindR (the stage 2 vs normal topTable output):
Stage_2_tT = read.xlsx("DGEA/Union/Stage_2_vs_Normal_z_DE_topTable.xlsx") %>%
  dplyr::select(Gene.Symbol, logFC, adj.P.Val) %>%
  na.omit()

# Preparing a pathfindR loop for enrichment analysis
dirs_stage_2 = paste0("pathfindR/Stage_2/", gene_sets)
pathfindR_outputs_stage_2 = list()

RNGversion("4.0.2")
set.seed(123)
for (i in 1:length(dirs_stage_2)){
  pathfindR_outputs_stage_2[[i]] = run_pathfindR(Stage_2_tT, gene_sets = gene_sets[i],
                                                 p_val_threshold = 0.05, 
                                                 visualize_enriched_terms = TRUE,
                                                 output_dir = dirs_stage_2[i], min_gset_size = 10,
                                                 max_gset_size = 300, adj_method = 'fdr',
                                                 enrichment_threshold = 0.05,
                                                 pin_name_path = 'Biogrid', search_method = 'GR',
                                                 grMaxDepth = 1, grSearchDepth = 1,
                                                 iterations = 10, n_processes = 10)
  cat(paste0("Done with ", gene_sets[i], "\n"))
}
names(pathfindR_outputs_stage_2) = gene_sets

# Perform hierarchical clustering on the results (average distance metric)
# Not run for GO-All and Reactome because the algorithm time complexity is O(n^3)

clustered_results_stage_2 = list()
for (i in c(1, 3, 4, 5)){
  clustered_results_stage_2[[i]] = cluster_enriched_terms(pathfindR_outputs_stage_2[[i]],
                                                          method = "hierarchical")
}
names(clustered_results_stage_2) = c("BioCarta", "NULL", "GO-CC", "GO-MF", "KEGG")

# BioCarta : The maximum average silhouette width was 0.16 for k = 97 
# GO-CC    : The maximum average silhouette width was 0.11 for k = 147 
# GO-MF    : The maximum average silhouette width was 0.1  for k = 188 
# KEGG     : The maximum average silhouette width was 0.13 for k = 91 

enrichment_dotplots_stage_2 = list()
cluster_enrichment_dotplots_stage_2 = list()

# Producing dotplots with the results
for (i in 1:length(pathfindR_outputs_stage_2)){
  # unclustered results
  enrichment_dotplots_stage_2[[i]] = enrichment_chart(result_df = pathfindR_outputs_stage_2[[i]],
                                                      top_terms = 10)
  tiff(paste0("pathfindR/Stage_2/", names(pathfindR_outputs_stage_2)[i], "/",
              names(pathfindR_outputs_stage_2)[i], "_top10_dotplot.tif"), 
       width = 1920, height = 1080, res = 150)
  print(enrichment_dotplots_stage_2[[i]])
  dev.off()
  
  if (i == 1 | i == 3 | i == 4 | i == 5){
    # clustered results
    cluster_enrichment_dotplots_stage_2[[i]] = enrichment_chart(result_df = clustered_results_stage_2[[i]][clustered_results_stage_2[[i]]$Status
                                                                                                           == "Representative", ][1:10,],
                                                                top_terms = NULL,
                                                                plot_by_cluster = TRUE)
    tiff(paste0("pathfindR/Stage_2/", names(pathfindR_outputs_stage_2)[i], "/",
                names(pathfindR_outputs_stage_2)[i], "_top10_dotplot_clustered.tif"), 
         width = 1920, height = 1080, res = 150)
    print(cluster_enrichment_dotplots_stage_2[[i]])
    dev.off()
  }
}

names(enrichment_dotplots_stage_2) = names(pathfindR_outputs_stage_2)
names(cluster_enrichment_dotplots_stage_2) = names(clustered_results_stage_2)

# Write out results in a comprehensive .xlsx file
wb = createWorkbook()
addWorksheet(wb, "BioCarta")
writeData(wb, "BioCarta", clustered_results_stage_2[["BioCarta"]])
addWorksheet(wb, "GO-BP")
writeData(wb, "GO-BP", pathfindR_outputs_stage_2[["GO-BP"]])
addWorksheet(wb, "GO-CC")
writeData(wb, "GO-CC", pathfindR_outputs_stage_2[["GO-CC"]])
addWorksheet(wb, "GO-MF")
writeData(wb, "GO-MF", pathfindR_outputs_stage_2[["GO-MF"]])
addWorksheet(wb, "KEGG")
writeData(wb, "KEGG", clustered_results_stage_2[["KEGG"]])
addWorksheet(wb, "Reactome")
writeData(wb, "Reactome", pathfindR_outputs_stage_2[["Reactome"]])
saveWorkbook(wb, file = "pathfindR/Stage_2/Comprehensive_pathfindR_output.xlsx",
             overwrite = TRUE); rm(wb)

# Representative terms output file (BioCarta, GO-CC, GO-MF & KEGG)
wb2 = createWorkbook()
addWorksheet(wb2, "BioCarta - rep")
writeData(wb2, "BioCarta - rep", clustered_results_stage_2[["BioCarta"]][clustered_results_stage_2[["BioCarta"]]$Status
                                                                         == "Representative", ])
addWorksheet(wb2, "GOCC - rep")
writeData(wb2, "GOCC - rep", clustered_results_stage_2[["GO-CC"]][clustered_results_stage_2[["GO-CC"]]$Status
                                                                  == "Representative", ])
addWorksheet(wb2, "GOMF - rep")
writeData(wb2, "GOMF - rep", clustered_results_stage_2[["GO-MF"]][clustered_results_stage_2[["GO-MF"]]$Status
                                                                  == "Representative", ])
addWorksheet(wb2, "KEGG - rep")
writeData(wb2, "KEGG - rep", clustered_results_stage_2[["KEGG"]][clustered_results_stage_2[["KEGG"]]$Status
                                                                 == "Representative", ])
saveWorkbook(wb2, file = "pathfindR/Stage_2/Representative_terms.xlsx",
             overwrite = TRUE); rm(wb2)

# Term-gene heatmaps and term-gene graphs #####
term_gene_heatmaps_stage_2 = list()
term_gene_graphs_stage_2 = list()

for (i in 1:length(pathfindR_outputs_stage_2)){
  # term-gene heatmaps
  term_gene_heatmaps_stage_2[[i]] = term_gene_heatmap(result_df = pathfindR_outputs_stage_2[[i]],
                                                      genes_df = Stage_2_tT,
                                                      num_terms = 5,
                                                      use_description = TRUE,
                                                      low = "darkgreen",
                                                      high = "darkred",
                                                      mid = "black")+
    theme(axis.text.x = element_text(size = axis_text_size[[i]]),
          axis.text.y = element_text(size = 12))
  tiff(paste0("pathfindR/Stage_2/", names(pathfindR_outputs_stage_2)[i], "/",
              names(pathfindR_outputs_stage_2)[i], "_top5_term_gene_heatmap.tif"), 
       width = 3840, height = 648, res = 150)
  print(term_gene_heatmaps_stage_2[[i]])
  dev.off()
  
  # term-gene graphs
  term_gene_graphs_stage_2[[i]] = term_gene_graph(result_df = pathfindR_outputs_stage_2[[i]],
                                                  num_terms = 3,
                                                  use_description = TRUE,
                                                  node_size = "p_val")
  tiff(paste0("pathfindR/Stage_2/", names(pathfindR_outputs_stage_2)[i], "/",
              names(pathfindR_outputs_stage_2)[i], "_top3_term_gene_graph.tif"), 
       width = 1920, height = 1080, res = 100)
  print(term_gene_graphs_stage_2[[i]])
  dev.off()
}

# UpSet plots #####
UpSet_plots_stage_2 = list()
for (i in 1:length(pathfindR_outputs_stage_2)){
  # term-gene heatmaps
  UpSet_plots_stage_2[[i]] = UpSet_plot(result_df = pathfindR_outputs_stage_2[[i]],
                                        genes_df = Stage_2_tT,
                                        num_terms = 5,
                                        use_description = TRUE,
                                        low = "darkgreen",
                                        high = "darkred",
                                        mid = "black")+
    theme(axis.text.y = element_text(size = axis_text_size[[i]]))
  tiff(paste0("pathfindR/Stage_2/", names(pathfindR_outputs_stage_2)[i], "/",
              names(pathfindR_outputs_stage_2)[i], "_top5_UpSet_plot.tif"), 
       width = 1444, height = 3840, res = 150)
  print(UpSet_plots_stage_2[[i]])
  dev.off()
}

names(term_gene_graphs_stage_2) = names(pathfindR_outputs_stage_2)
names(term_gene_heatmaps_stage_2) = names(pathfindR_outputs_stage_2)
names(UpSet_plots_stage_2) = names(pathfindR_outputs_stage_2)


# Stage 3 #####
# Loading the input to pathfindR (the stage 3 vs normal topTable output):
Stage_3_tT = read.xlsx("DGEA/Union/Stage_3_vs_Normal_z_DE_topTable.xlsx") %>%
  dplyr::select(Gene.Symbol, logFC, adj.P.Val) %>%
  na.omit()

# Preparing a pathfindR loop for enrichment analysis
dirs_stage_3 = paste0("pathfindR/Stage_3/", gene_sets)
pathfindR_outputs_stage_3 = list()

RNGversion("4.0.2")
set.seed(123)
for (i in 1:length(dirs_stage_3)){
  pathfindR_outputs_stage_3[[i]] = run_pathfindR(Stage_3_tT, gene_sets = gene_sets[i],
                                                 p_val_threshold = 0.05, 
                                                 visualize_enriched_terms = TRUE,
                                                 output_dir = dirs_stage_3[i], min_gset_size = 10,
                                                 max_gset_size = 300, adj_method = 'fdr',
                                                 enrichment_threshold = 0.05,
                                                 pin_name_path = 'Biogrid', search_method = 'GR',
                                                 grMaxDepth = 1, grSearchDepth = 1,
                                                 iterations = 10, n_processes = 10)
  cat(paste0("Done with ", gene_sets[i], "\n"))
}
names(pathfindR_outputs_stage_3) = gene_sets

# Perform hierarchical clustering on the results (average distance metric)
# Not run for GO-All and Reactome because the algorithm time complexity is O(n^3)

clustered_results_stage_3 = list()
for (i in c(1, 3, 4, 5)){
  clustered_results_stage_3[[i]] = cluster_enriched_terms(pathfindR_outputs_stage_3[[i]],
                                                          method = "hierarchical")
}
names(clustered_results_stage_3) = c("BioCarta", "NULL", "GO-CC", "GO-MF", "KEGG")

# BioCarta: The maximum average silhouette width was 0.16 for k = 74 
# GO-CC   : The maximum average silhouette width was 0.11 for k = 156
# GO-MF   : The maximum average silhouette width was 0.1  for k = 165
# KEGG    : The maximum average silhouette width was 0.13 for k = 113

enrichment_dotplots_stage_3 = list()
cluster_enrichment_dotplots_stage_3 = list()

# Producing dotplots with the results
for (i in 1:length(pathfindR_outputs_stage_3)){
  # unclustered results
  enrichment_dotplots_stage_3[[i]] = enrichment_chart(result_df = pathfindR_outputs_stage_3[[i]],
                                                      top_terms = 10)
  tiff(paste0("pathfindR/Stage_3/", names(pathfindR_outputs_stage_3)[i], "/",
              names(pathfindR_outputs_stage_3)[i], "_top10_dotplot.tif"), 
       width = 1920, height = 1080, res = 150)
  print(enrichment_dotplots_stage_3[[i]])
  dev.off()
  
  if (i == 1 | i == 3 | i == 4 | i == 5){
    # clustered results
    cluster_enrichment_dotplots_stage_3[[i]] = enrichment_chart(result_df = clustered_results_stage_3[[i]][clustered_results_stage_3[[i]]$Status
                                                                                                           == "Representative", ][1:10,],
                                                                top_terms = NULL,
                                                                plot_by_cluster = TRUE)
    tiff(paste0("pathfindR/Stage_3/", names(pathfindR_outputs_stage_3)[i], "/",
                names(pathfindR_outputs_stage_3)[i], "_top10_dotplot_clustered.tif"), 
         width = 1920, height = 1080, res = 150)
    print(cluster_enrichment_dotplots_stage_3[[i]])
    dev.off()
  }
}

names(enrichment_dotplots_stage_3) = names(pathfindR_outputs_stage_3)
names(cluster_enrichment_dotplots_stage_3) = names(clustered_results_stage_3)

# Write out results in a comprehensive .xlsx file
wb = createWorkbook()
addWorksheet(wb, "BioCarta")
writeData(wb, "BioCarta", clustered_results_stage_3[["BioCarta"]])
addWorksheet(wb, "GO-BP")
writeData(wb, "GO-BP", pathfindR_outputs_stage_3[["GO-BP"]])
addWorksheet(wb, "GO-CC")
writeData(wb, "GO-CC", pathfindR_outputs_stage_3[["GO-CC"]])
addWorksheet(wb, "GO-MF")
writeData(wb, "GO-MF", pathfindR_outputs_stage_3[["GO-MF"]])
addWorksheet(wb, "KEGG")
writeData(wb, "KEGG", clustered_results_stage_3[["KEGG"]])
addWorksheet(wb, "Reactome")
writeData(wb, "Reactome", pathfindR_outputs_stage_3[["Reactome"]])
saveWorkbook(wb, file = "pathfindR/Stage_3/Comprehensive_pathfindR_output.xlsx",
             overwrite = TRUE); rm(wb)

# Representative terms output file (BioCarta, GO-CC, GO-MF & KEGG)
wb2 = createWorkbook()
addWorksheet(wb2, "BioCarta - rep")
writeData(wb2, "BioCarta - rep", clustered_results_stage_3[["BioCarta"]][clustered_results_stage_3[["BioCarta"]]$Status
                                                                         == "Representative", ])
addWorksheet(wb2, "GOCC - rep")
writeData(wb2, "GOCC - rep", clustered_results_stage_3[["GO-CC"]][clustered_results_stage_3[["GO-CC"]]$Status
                                                                  == "Representative", ])
addWorksheet(wb2, "GOMF - rep")
writeData(wb2, "GOMF - rep", clustered_results_stage_3[["GO-MF"]][clustered_results_stage_3[["GO-MF"]]$Status
                                                                  == "Representative", ])
addWorksheet(wb2, "KEGG - rep")
writeData(wb2, "KEGG - rep", clustered_results_stage_3[["KEGG"]][clustered_results_stage_3[["KEGG"]]$Status
                                                                 == "Representative", ])
saveWorkbook(wb2, file = "pathfindR/Stage_3/Representative_terms.xlsx",
             overwrite = TRUE); rm(wb2)

# Term-gene heatmaps and term-gene graphs #####
term_gene_heatmaps_stage_3 = list()
term_gene_graphs_stage_3 = list()

for (i in 1:length(pathfindR_outputs_stage_3)){
  # term-gene heatmaps
  term_gene_heatmaps_stage_3[[i]] = term_gene_heatmap(result_df = pathfindR_outputs_stage_3[[i]],
                                                      genes_df = Stage_3_tT,
                                                      num_terms = 5,
                                                      use_description = TRUE,
                                                      low = "darkgreen",
                                                      high = "darkred",
                                                      mid = "black")+
    theme(axis.text.x = element_text(size = axis_text_size[[i]]),
          axis.text.y = element_text(size = 12))
  tiff(paste0("pathfindR/Stage_3/", names(pathfindR_outputs_stage_3)[i], "/",
              names(pathfindR_outputs_stage_3)[i], "_top5_term_gene_heatmap.tif"), 
       width = 3840, height = 648, res = 150)
  print(term_gene_heatmaps_stage_3[[i]])
  dev.off()
  
  # term-gene graphs
  term_gene_graphs_stage_3[[i]] = term_gene_graph(result_df = pathfindR_outputs_stage_3[[i]],
                                                  num_terms = 3,
                                                  use_description = TRUE,
                                                  node_size = "p_val")
  tiff(paste0("pathfindR/Stage_3/", names(pathfindR_outputs_stage_3)[i], "/",
              names(pathfindR_outputs_stage_3)[i], "_top3_term_gene_graph.tif"), 
       width = 1920, height = 1080, res = 100)
  print(term_gene_graphs_stage_3[[i]])
  dev.off()
}

# UpSet plots #####
UpSet_plots_stage_3 = list()
for (i in 1:length(pathfindR_outputs_stage_3)){
  # term-gene heatmaps
  UpSet_plots_stage_3[[i]] = UpSet_plot(result_df = pathfindR_outputs_stage_3[[i]],
                                        genes_df = Stage_3_tT,
                                        num_terms = 5,
                                        use_description = TRUE,
                                        low = "darkgreen",
                                        high = "darkred",
                                        mid = "black")+
    theme(axis.text.y = element_text(size = axis_text_size[[i]]))
  tiff(paste0("pathfindR/Stage_3/", names(pathfindR_outputs_stage_3)[i], "/",
              names(pathfindR_outputs_stage_3)[i], "_top5_UpSet_plot.tif"), 
       width = 1444, height = 3840, res = 150)
  print(UpSet_plots_stage_3[[i]])
  dev.off()
}

names(term_gene_graphs_stage_3) = names(pathfindR_outputs_stage_3)
names(term_gene_heatmaps_stage_3) = names(pathfindR_outputs_stage_3)
names(UpSet_plots_stage_3) = names(pathfindR_outputs_stage_3)

# Stage 4 #####
# Loading the input to pathfindR (the stage 4 vs normal topTable output):
Stage_4_tT = read.xlsx("DGEA/Union/Stage_4_vs_Normal_z_DE_topTable.xlsx") %>%
  dplyr::select(Gene.Symbol, logFC, adj.P.Val) %>%
  na.omit()

# Preparing a pathfindR loop for enrichment analysis
dirs_stage_4 = paste0("pathfindR/Stage_4/", gene_sets)
pathfindR_outputs_stage_4 = list()

RNGversion("4.0.2")
set.seed(123)
for (i in 1:length(dirs_stage_4)){
  pathfindR_outputs_stage_4[[i]] = run_pathfindR(Stage_4_tT, gene_sets = gene_sets[i],
                                                 p_val_threshold = 0.05, 
                                                 visualize_enriched_terms = TRUE,
                                                 output_dir = dirs_stage_4[i], min_gset_size = 10,
                                                 max_gset_size = 300, adj_method = 'fdr',
                                                 enrichment_threshold = 0.05,
                                                 pin_name_path = 'Biogrid', search_method = 'GR',
                                                 grMaxDepth = 1, grSearchDepth = 1,
                                                 iterations = 10, n_processes = 10)
  cat(paste0("Done with ", gene_sets[i], "\n"))
}
names(pathfindR_outputs_stage_4) = gene_sets

# Perform hierarchical clustering on the results (average distance metric)
# Not run for GO-All and Reactome because the algorithm time complexity is O(n^3)

clustered_results_stage_4 = list()
for (i in c(1, 3, 4, 5)){
  clustered_results_stage_4[[i]] = cluster_enriched_terms(pathfindR_outputs_stage_4[[i]],
                                                          method = "hierarchical")
}
names(clustered_results_stage_4) = c("BioCarta", "NULL", "GO-CC", "GO-MF", "KEGG")

# BioCarta: The maximum average silhouette width was 0.15 for k = 43 
# GO-CC   : The maximum average silhouette width was 0.11 for k = 162
# GO-MF   : The maximum average silhouette width was 0.1  for k = 208
# KEGG    : The maximum average silhouette width was 0.14 for k = 112

enrichment_dotplots_stage_4 = list()
cluster_enrichment_dotplots_stage_4 = list()

# Producing dotplots with the results
for (i in 1:length(pathfindR_outputs_stage_4)){
  # unclustered results
  enrichment_dotplots_stage_4[[i]] = enrichment_chart(result_df = pathfindR_outputs_stage_4[[i]],
                                                      top_terms = 10)
  tiff(paste0("pathfindR/Stage_4/", names(pathfindR_outputs_stage_4)[i], "/",
              names(pathfindR_outputs_stage_4)[i], "_top10_dotplot.tif"), 
       width = 1920, height = 1080, res = 150)
  print(enrichment_dotplots_stage_4[[i]])
  dev.off()
  
  if (i == 1 | i == 3 | i == 4 | i == 5){
    # clustered results
    cluster_enrichment_dotplots_stage_4[[i]] = enrichment_chart(result_df = clustered_results_stage_4[[i]][clustered_results_stage_4[[i]]$Status
                                                                                                           == "Representative", ][1:10,],
                                                                top_terms = NULL,
                                                                plot_by_cluster = TRUE)
    tiff(paste0("pathfindR/Stage_4/", names(pathfindR_outputs_stage_4)[i], "/",
                names(pathfindR_outputs_stage_4)[i], "_top10_dotplot_clustered.tif"), 
         width = 1920, height = 1080, res = 150)
    print(cluster_enrichment_dotplots_stage_4[[i]])
    dev.off()
  }
}

names(enrichment_dotplots_stage_4) = names(pathfindR_outputs_stage_4)
names(cluster_enrichment_dotplots_stage_4) = names(clustered_results_stage_4)

# Write out results in a comprehensive .xlsx file
wb = createWorkbook()
addWorksheet(wb, "BioCarta")
writeData(wb, "BioCarta", clustered_results_stage_4[["BioCarta"]])
addWorksheet(wb, "GO-BP")
writeData(wb, "GO-BP", pathfindR_outputs_stage_4[["GO-BP"]])
addWorksheet(wb, "GO-CC")
writeData(wb, "GO-CC", pathfindR_outputs_stage_4[["GO-CC"]])
addWorksheet(wb, "GO-MF")
writeData(wb, "GO-MF", pathfindR_outputs_stage_4[["GO-MF"]])
addWorksheet(wb, "KEGG")
writeData(wb, "KEGG", clustered_results_stage_4[["KEGG"]])
addWorksheet(wb, "Reactome")
writeData(wb, "Reactome", pathfindR_outputs_stage_4[["Reactome"]])
saveWorkbook(wb, file = "pathfindR/Stage_4/Comprehensive_pathfindR_output.xlsx",
             overwrite = TRUE); rm(wb)

# Representative terms output file (BioCarta, GO-CC, GO-MF & KEGG)
wb2 = createWorkbook()
addWorksheet(wb2, "BioCarta - rep")
writeData(wb2, "BioCarta - rep", clustered_results_stage_4[["BioCarta"]][clustered_results_stage_4[["BioCarta"]]$Status
                                                                         == "Representative", ])
addWorksheet(wb2, "GOCC - rep")
writeData(wb2, "GOCC - rep", clustered_results_stage_4[["GO-CC"]][clustered_results_stage_4[["GO-CC"]]$Status
                                                                  == "Representative", ])
addWorksheet(wb2, "GOMF - rep")
writeData(wb2, "GOMF - rep", clustered_results_stage_4[["GO-MF"]][clustered_results_stage_4[["GO-MF"]]$Status
                                                                  == "Representative", ])
addWorksheet(wb2, "KEGG - rep")
writeData(wb2, "KEGG - rep", clustered_results_stage_4[["KEGG"]][clustered_results_stage_4[["KEGG"]]$Status
                                                                 == "Representative", ])
saveWorkbook(wb2, file = "pathfindR/Stage_4/Representative_terms.xlsx",
             overwrite = TRUE); rm(wb2)

# Term-gene heatmaps and term-gene graphs #####
term_gene_heatmaps_stage_4 = list()
term_gene_graphs_stage_4 = list()

for (i in 1:length(pathfindR_outputs_stage_4)){
  # term-gene heatmaps
  term_gene_heatmaps_stage_4[[i]] = term_gene_heatmap(result_df = pathfindR_outputs_stage_4[[i]],
                                                      genes_df = Stage_4_tT,
                                                      num_terms = 5,
                                                      use_description = TRUE,
                                                      low = "darkgreen",
                                                      high = "darkred",
                                                      mid = "black")+
    theme(axis.text.x = element_text(size = axis_text_size[[i]]),
          axis.text.y = element_text(size = 12))
  tiff(paste0("pathfindR/Stage_4/", names(pathfindR_outputs_stage_4)[i], "/",
              names(pathfindR_outputs_stage_4)[i], "_top5_term_gene_heatmap.tif"), 
       width = 3840, height = 648, res = 150)
  print(term_gene_heatmaps_stage_4[[i]])
  dev.off()
  
  # term-gene graphs
  term_gene_graphs_stage_4[[i]] = term_gene_graph(result_df = pathfindR_outputs_stage_4[[i]],
                                                  num_terms = 3,
                                                  use_description = TRUE,
                                                  node_size = "p_val")
  tiff(paste0("pathfindR/Stage_4/", names(pathfindR_outputs_stage_4)[i], "/",
              names(pathfindR_outputs_stage_4)[i], "_top3_term_gene_graph.tif"), 
       width = 1920, height = 1080, res = 100)
  print(term_gene_graphs_stage_4[[i]])
  dev.off()
}

# UpSet plots #####
UpSet_plots_stage_4 = list()
for (i in 1:length(pathfindR_outputs_stage_4)){
  # term-gene heatmaps
  UpSet_plots_stage_4[[i]] = UpSet_plot(result_df = pathfindR_outputs_stage_4[[i]],
                                        genes_df = Stage_4_tT,
                                        num_terms = 5,
                                        use_description = TRUE,
                                        low = "darkgreen",
                                        high = "darkred",
                                        mid = "black")+
    theme(axis.text.y = element_text(size = axis_text_size[[i]]))
  tiff(paste0("pathfindR/Stage_4/", names(pathfindR_outputs_stage_4)[i], "/",
              names(pathfindR_outputs_stage_4)[i], "_top5_UpSet_plot.tif"), 
       width = 1444, height = 3840, res = 150)
  print(UpSet_plots_stage_4[[i]])
  dev.off()
}

names(term_gene_graphs_stage_4) = names(pathfindR_outputs_stage_4)
names(term_gene_heatmaps_stage_4) = names(pathfindR_outputs_stage_4)
names(UpSet_plots_stage_4) = names(pathfindR_outputs_stage_4)
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

sig_DGEA_stage_1 = read.xlsx("DGEA/Union/Stage_1_vs_Normal_z_DE_topTable.xlsx") %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::rename(Gene.Symbol_pre = Gene.Symbol)

sig_DGEA_stage_2 = read.xlsx("DGEA/Union/Stage_2_vs_Normal_z_DE_topTable.xlsx") %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::rename(Gene.Symbol_pre = Gene.Symbol)

sig_DGEA_stage_3 = read.xlsx("DGEA/Union/Stage_3_vs_Normal_z_DE_topTable.xlsx") %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::rename(Gene.Symbol_pre = Gene.Symbol)

sig_DGEA_stage_4 = read.xlsx("DGEA/Union/Stage_4_vs_Normal_z_DE_topTable.xlsx") %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::rename(Gene.Symbol_pre = Gene.Symbol)

sig_DGEA = list(sig_DGEA_stage_1, sig_DGEA_stage_2, sig_DGEA_stage_3, sig_DGEA_stage_4)
names(sig_DGEA) = c("Stage_1", "Stage_2", "Stage_3", "Stage_4")
rm(sig_DGEA_stage_1, sig_DGEA_stage_2, sig_DGEA_stage_3, sig_DGEA_stage_4); gc()

# Drivers lists
general_drivers = list()
PDAC_drivers = list()

for (i in 1:4){
  general_drivers[[i]] = inner_join(census, sig_DGEA[[i]], by = "EntrezGene.ID") %>%
    dplyr::select(EntrezGene.ID, Gene.Symbol_pre, Gene.Symbol_COSMIC, Name, logFC, 
                  P.Value, B, HGNC_Official, adj.P.Val, everything())
  PDAC_drivers[[i]] = inner_join(census, sig_DGEA[[i]], by = "EntrezGene.ID") %>%
    dplyr::filter(grepl("pancr", Tumour.Types.Somatic.) | grepl("pancr", Tumour.Types.Germline.)) %>%
    dplyr::select(EntrezGene.ID, Gene.Symbol_pre, Gene.Symbol_COSMIC, Name, logFC, 
                  P.Value, B, HGNC_Official, adj.P.Val, everything())
}

names(general_drivers) = names(sig_DGEA)
names(PDAC_drivers) = names(sig_DGEA)

# Writing out
for (i in 1:4){
  drivers = createWorkbook()
  addWorksheet(drivers, "General")
  writeData(drivers, "General", general_drivers[[i]])
  addWorksheet(drivers, "Pancreas")
  writeData(drivers, "Pancreas", PDAC_drivers[[i]])
  saveWorkbook(drivers, paste0("DGEA/Drivers/", names(sig_DGEA)[i], "_Drivers.xlsx"),
               overwrite = TRUE); rm(drivers); gc()
}

rm(i)

##### Plotting #####
library(ggVennDiagram)
library(ggvenn)
library(EnhancedVolcano)

# ggVennDiagram
# Generating a suitable object with stat. sig (p.adj < 0.05) DEGs from all stages
venn_DEG = list(`Stage 1` = sig_DGEA[["Stage_1"]]$EntrezGene.ID[sig_DGEA[["Stage_1"]]$adj.P.Val < 0.05],
            `Stage 2` = sig_DGEA[["Stage_2"]]$EntrezGene.ID[sig_DGEA[["Stage_2"]]$adj.P.Val < 0.05],
            `Stage 3` = sig_DGEA[["Stage_3"]]$EntrezGene.ID[sig_DGEA[["Stage_3"]]$adj.P.Val < 0.05],
            `Stage 4` = sig_DGEA[["Stage_4"]]$EntrezGene.ID[sig_DGEA[["Stage_4"]]$adj.P.Val < 0.05])

diagram1 = ggVennDiagram(
  venn_DEG, label_alpha = 0,
  category.names = names(venn_DEG)
) +
  scale_fill_gradient(low = "white", high = "darkred") +
  theme(plot.title = element_text(face = "bold"))+
  labs(title = "Venn diagram of significantly differentially expressed genes (DEGs) between stages")
tiff("Additional_plots/Venn/ggVennDiagram_Stages_DEG_Venn.tif", 
     width = 1920, height = 1080, res = 200)
diagram1
dev.off(); rm(diagram1)

# ggvenn
diagram2 = ggvenn(
  venn_DEG, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)+
  theme(plot.title = element_text(face = "bold"))+
  labs(title = "Venn diagram of significantly differentially expressed genes (DEGs) between stages")
tiff("Additional_plots/Venn/ggvenn_Stages_DEG_Venn.tif", 
     width = 1920, height = 1080, res = 150)
diagram2
dev.off(); rm(diagram2)

# Defining a multiplot function #####
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Create a 2x3 Venn diagram multiplot that compares the subnetwork enrichment
# between 4 stages with respect to BioCarta, GO-BP, GO-CC, GO-MF, KEGG, Reactome
# 3 versions: all results, Top 100, clustered results (representative terms)
# All will be produced by ggVennDiagram

# Producing the necessary lists:
venn_BioCarta = list(pathfindR_outputs_stage_1[["BioCarta"]]$Term_Description,
                     pathfindR_outputs_stage_2[["BioCarta"]]$Term_Description,
                     pathfindR_outputs_stage_3[["BioCarta"]]$Term_Description,
                     pathfindR_outputs_stage_4[["BioCarta"]]$Term_Description)
names(venn_BioCarta) = names(venn_DEG)

venn_GO_BP = list(pathfindR_outputs_stage_1[["GO-BP"]]$Term_Description,
                  pathfindR_outputs_stage_2[["GO-BP"]]$Term_Description,
                  pathfindR_outputs_stage_3[["GO-BP"]]$Term_Description,
                  pathfindR_outputs_stage_4[["GO-BP"]]$Term_Description)
names(venn_GO_BP) = names(venn_DEG)

venn_GO_CC = list(pathfindR_outputs_stage_1[["GO-CC"]]$Term_Description,
                  pathfindR_outputs_stage_2[["GO-CC"]]$Term_Description,
                  pathfindR_outputs_stage_3[["GO-CC"]]$Term_Description,
                  pathfindR_outputs_stage_4[["GO-CC"]]$Term_Description)
names(venn_GO_CC) = names(venn_DEG)

venn_GO_MF = list(pathfindR_outputs_stage_1[["GO-MF"]]$Term_Description,
                  pathfindR_outputs_stage_2[["GO-MF"]]$Term_Description,
                  pathfindR_outputs_stage_3[["GO-MF"]]$Term_Description,
                  pathfindR_outputs_stage_4[["GO-MF"]]$Term_Description)
names(venn_GO_MF) = names(venn_DEG)

venn_KEGG = list(pathfindR_outputs_stage_1[["KEGG"]]$Term_Description,
                  pathfindR_outputs_stage_2[["KEGG"]]$Term_Description,
                  pathfindR_outputs_stage_3[["KEGG"]]$Term_Description,
                  pathfindR_outputs_stage_4[["KEGG"]]$Term_Description)
names(venn_KEGG) = names(venn_DEG)

venn_Reactome = list(pathfindR_outputs_stage_1[["Reactome"]]$Term_Description,
                 pathfindR_outputs_stage_2[["Reactome"]]$Term_Description,
                 pathfindR_outputs_stage_3[["Reactome"]]$Term_Description,
                 pathfindR_outputs_stage_4[["Reactome"]]$Term_Description)
names(venn_Reactome) = names(venn_DEG)

venn_path = list(venn_BioCarta, venn_GO_BP, venn_GO_CC, venn_GO_MF, venn_KEGG,
                 venn_Reactome)
names(venn_path) = names(pathfindR_outputs_stage_1)
rm(venn_BioCarta, venn_GO_BP, venn_GO_CC, venn_GO_MF, venn_KEGG,
   venn_Reactome)

# top100 list
top100_venn_BioCarta = list(pathfindR_outputs_stage_1[["BioCarta"]]$Term_Description[1:100],
                     pathfindR_outputs_stage_2[["BioCarta"]]$Term_Description[1:100],
                     pathfindR_outputs_stage_3[["BioCarta"]]$Term_Description[1:100],
                     pathfindR_outputs_stage_4[["BioCarta"]]$Term_Description[1:100])
names(top100_venn_BioCarta) = names(venn_DEG)

top100_venn_GO_BP = list(pathfindR_outputs_stage_1[["GO-BP"]]$Term_Description[1:100],
                  pathfindR_outputs_stage_2[["GO-BP"]]$Term_Description[1:100],
                  pathfindR_outputs_stage_3[["GO-BP"]]$Term_Description[1:100],
                  pathfindR_outputs_stage_4[["GO-BP"]]$Term_Description[1:100])
names(top100_venn_GO_BP) = names(venn_DEG)

top100_venn_GO_CC = list(pathfindR_outputs_stage_1[["GO-CC"]]$Term_Description[1:100],
                  pathfindR_outputs_stage_2[["GO-CC"]]$Term_Description[1:100],
                  pathfindR_outputs_stage_3[["GO-CC"]]$Term_Description[1:100],
                  pathfindR_outputs_stage_4[["GO-CC"]]$Term_Description[1:100])
names(top100_venn_GO_CC) = names(venn_DEG)

top100_venn_GO_MF = list(pathfindR_outputs_stage_1[["GO-MF"]]$Term_Description[1:100],
                  pathfindR_outputs_stage_2[["GO-MF"]]$Term_Description[1:100],
                  pathfindR_outputs_stage_3[["GO-MF"]]$Term_Description[1:100],
                  pathfindR_outputs_stage_4[["GO-MF"]]$Term_Description[1:100])
names(top100_venn_GO_MF) = names(venn_DEG)

top100_venn_KEGG = list(pathfindR_outputs_stage_1[["KEGG"]]$Term_Description[1:100],
                 pathfindR_outputs_stage_2[["KEGG"]]$Term_Description[1:100],
                 pathfindR_outputs_stage_3[["KEGG"]]$Term_Description[1:100],
                 pathfindR_outputs_stage_4[["KEGG"]]$Term_Description[1:100])
names(top100_venn_KEGG) = names(venn_DEG)

top100_venn_Reactome = list(pathfindR_outputs_stage_1[["Reactome"]]$Term_Description[1:100],
                     pathfindR_outputs_stage_2[["Reactome"]]$Term_Description[1:100],
                     pathfindR_outputs_stage_3[["Reactome"]]$Term_Description[1:100],
                     pathfindR_outputs_stage_4[["Reactome"]]$Term_Description[1:100])
names(top100_venn_Reactome) = names(venn_DEG)

top100_venn_path = list(top100_venn_BioCarta, top100_venn_GO_BP, top100_venn_GO_CC, 
                        top100_venn_GO_MF, top100_venn_KEGG, top100_venn_Reactome)
names(top100_venn_path) = names(pathfindR_outputs_stage_1)
rm(top100_venn_BioCarta, top100_venn_GO_BP, top100_venn_GO_CC, top100_venn_GO_MF,
   top100_venn_KEGG, top100_venn_Reactome)

# clustered list
clustered_venn_BioCarta = list(clustered_results_stage_1[["BioCarta"]]$Term_Description[clustered_results_stage_1[["BioCarta"]]$Status == "Representative"],
                               clustered_results_stage_2[["BioCarta"]]$Term_Description[clustered_results_stage_2[["BioCarta"]]$Status == "Representative"],
                               clustered_results_stage_3[["BioCarta"]]$Term_Description[clustered_results_stage_3[["BioCarta"]]$Status == "Representative"],
                               clustered_results_stage_4[["BioCarta"]]$Term_Description[clustered_results_stage_4[["BioCarta"]]$Status == "Representative"])
names(clustered_venn_BioCarta) = names(venn_DEG)

clustered_venn_GO_CC = list(clustered_results_stage_1[["GO-CC"]]$Term_Description[clustered_results_stage_1[["GO-CC"]]$Status == "Representative"],
                            clustered_results_stage_2[["GO-CC"]]$Term_Description[clustered_results_stage_2[["GO-CC"]]$Status == "Representative"],
                            clustered_results_stage_3[["GO-CC"]]$Term_Description[clustered_results_stage_3[["GO-CC"]]$Status == "Representative"],
                            clustered_results_stage_4[["GO-CC"]]$Term_Description[clustered_results_stage_4[["GO-CC"]]$Status == "Representative"])
names(clustered_venn_GO_CC) = names(venn_DEG)

clustered_venn_GO_MF = list(clustered_results_stage_1[["GO-MF"]]$Term_Description[clustered_results_stage_1[["GO-MF"]]$Status == "Representative"],
                            clustered_results_stage_2[["GO-MF"]]$Term_Description[clustered_results_stage_2[["GO-MF"]]$Status == "Representative"],
                            clustered_results_stage_3[["GO-MF"]]$Term_Description[clustered_results_stage_3[["GO-MF"]]$Status == "Representative"],
                            clustered_results_stage_4[["GO-MF"]]$Term_Description[clustered_results_stage_4[["GO-MF"]]$Status == "Representative"])
names(clustered_venn_GO_MF) = names(venn_DEG)

clustered_venn_KEGG = list(clustered_results_stage_1[["KEGG"]]$Term_Description[clustered_results_stage_1[["KEGG"]]$Status == "Representative"],
                           clustered_results_stage_2[["KEGG"]]$Term_Description[clustered_results_stage_2[["KEGG"]]$Status == "Representative"],
                           clustered_results_stage_3[["KEGG"]]$Term_Description[clustered_results_stage_3[["KEGG"]]$Status == "Representative"],
                           clustered_results_stage_4[["KEGG"]]$Term_Description[clustered_results_stage_4[["KEGG"]]$Status == "Representative"])
names(clustered_venn_KEGG) = names(venn_DEG)

clustered_venn_path = list(clustered_venn_BioCarta, clustered_venn_GO_CC, 
                           clustered_venn_GO_MF, clustered_venn_KEGG)
names(clustered_venn_path) = c("BioCarta", "GO-CC", "GO-MF", "KEGG")
rm(clustered_venn_BioCarta, clustered_venn_GO_CC, clustered_venn_GO_MF,
   clustered_venn_KEGG)

# Preparing lists for "for" loops
all_results_path_venn = list()
top100_results_path_venn = list()
clustered_results_path_venn = list()

# All results
for (i in 1:length(venn_path)){
  all_results_path_venn[[i]] = ggVennDiagram(
    venn_path[[i]], label_alpha = 0,
    category.names = names(venn_path[[i]])
  ) +
    scale_fill_gradient(low = "white", high = "darkred") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))+
    labs(title = paste0("Venn diagram of ", names(venn_path)[i], " enriched networks"))
}

tiff("Additional_plots/Venn/all_path_results_ggVennDiagram_multiplot.tif", 
     width = 2880, height = 1620, res = 150)
multiplot(all_results_path_venn[[1]], all_results_path_venn[[2]],
              all_results_path_venn[[3]], all_results_path_venn[[4]],
              all_results_path_venn[[5]], all_results_path_venn[[6]], cols = 3)
m = ggplot(multiplot(all_results_path_venn[[1]], all_results_path_venn[[2]],
                     all_results_path_venn[[3]], all_results_path_venn[[4]],
                     all_results_path_venn[[5]], all_results_path_venn[[6]], cols = 3))
dev.off(); rm(m)

# Top 100
for (i in 1:length(top100_venn_path)){
  top100_results_path_venn[[i]] = ggVennDiagram(
    top100_venn_path[[i]], label_alpha = 0,
    category.names = names(top100_venn_path[[i]])
  ) +
    scale_fill_gradient(low = "white", high = "darkred") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))+
    labs(title = paste0("Venn diagram of top 100 ", names(top100_venn_path)[i], 
                        " enriched networks"))
}

tiff("Additional_plots/Venn/top100_path_results_ggVennDiagram_multiplot.tif", 
     width = 2880, height = 1620, res = 150)
multiplot(top100_results_path_venn[[1]], top100_results_path_venn[[2]],
          top100_results_path_venn[[3]], top100_results_path_venn[[4]],
          top100_results_path_venn[[5]], top100_results_path_venn[[6]], cols = 3)
m = ggplot(multiplot(top100_results_path_venn[[1]], top100_results_path_venn[[2]],
                     top100_results_path_venn[[3]], top100_results_path_venn[[4]],
                     top100_results_path_venn[[5]], top100_results_path_venn[[6]], cols = 3))
dev.off(); rm(m)

# clustered enriched networks multiplot (not for GO-BP, Reactome)
# Top 100
for (i in 1:length(clustered_venn_path)){
  clustered_results_path_venn[[i]] = ggVennDiagram(
    clustered_venn_path[[i]], label_alpha = 0,
    category.names = names(clustered_venn_path[[i]])
  ) +
    scale_fill_gradient(low = "white", high = "darkred") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))+
    labs(title = paste0("Venn diagram of clustered ", names(top100_venn_path)[i], 
                        " enriched networks"))
}

tiff("Additional_plots/Venn/clustered_path_results_ggVennDiagram_multiplot.tif", 
     width = 2880, height = 1620, res = 150)
multiplot(clustered_results_path_venn[[1]], clustered_results_path_venn[[2]],
          clustered_results_path_venn[[3]], clustered_results_path_venn[[4]], cols = 2)
m = ggplot(multiplot(clustered_results_path_venn[[1]], clustered_results_path_venn[[2]],
                     clustered_results_path_venn[[3]], clustered_results_path_venn[[4]], cols = 2))
dev.off(); rm(m)

#####
# Load volcano plots from tumor stage DGEA as an .RData file:
# load("/Your/path/your_plots.RData")

# All-stage-volcano-dotplot (BioCarta)-pairs
tiff("Additional_plots/all_Volcano_CE_Dotplot_multiplot.tif", 
     width = 7680, height = 4320, res = 100)
multiplot(union_one_normal_volcano, union_two_normal_volcano, 
              union_three_normal_volcano, union_four_normal_volcano, 
          cluster_enrichment_dotplots_stage_1[["BioCarta"]],
          cluster_enrichment_dotplots_stage_2[["BioCarta"]],
          cluster_enrichment_dotplots_stage_3[["BioCarta"]],
          cluster_enrichment_dotplots_stage_4[["BioCarta"]], cols = 2)
m = ggplot(multiplot(union_one_normal_volcano, union_two_normal_volcano, 
                     union_three_normal_volcano, union_four_normal_volcano, 
                     cluster_enrichment_dotplots_stage_1[["BioCarta"]],
                     cluster_enrichment_dotplots_stage_2[["BioCarta"]],
                     cluster_enrichment_dotplots_stage_3[["BioCarta"]],
                     cluster_enrichment_dotplots_stage_4[["BioCarta"]], cols = 2))
dev.off(); rm(m)

# 2 pairs: Stages 1-2, Stages 3-4
# Stages 1-2
tiff("Additional_plots/Stages_1_2_Volcano_CE_Dotplot_multiplot.tif", 
     width = 3840, height = 2160, res = 150)
multiplot(union_one_normal_volcano, union_two_normal_volcano, 
          cluster_enrichment_dotplots_stage_1[["BioCarta"]],
          cluster_enrichment_dotplots_stage_2[["BioCarta"]], cols = 2)
m = ggplot(multiplot(union_one_normal_volcano, union_two_normal_volcano, 
                     cluster_enrichment_dotplots_stage_1[["BioCarta"]],
                     cluster_enrichment_dotplots_stage_2[["BioCarta"]], cols = 2))
dev.off(); rm(m)

# Stages 3-4
tiff("Additional_plots/Stages_3_4_Volcano_CE_Dotplot_multiplot.tif", 
     width = 3840, height = 2160, res = 150)
multiplot(union_three_normal_volcano, union_four_normal_volcano, 
          cluster_enrichment_dotplots_stage_3[["BioCarta"]],
          cluster_enrichment_dotplots_stage_4[["BioCarta"]], cols = 2)
m = ggplot(multiplot(union_three_normal_volcano, union_four_normal_volcano, 
                     cluster_enrichment_dotplots_stage_3[["BioCarta"]],
                     cluster_enrichment_dotplots_stage_4[["BioCarta"]], cols = 2))
dev.off(); rm(m)

# Unique pairs
# Stage 1
tiff("Additional_plots/Stage_1_Volcano_CE_Dotplot.tif", 
     width = 1920, height = 1080, res = 150)
multiplot(union_one_normal_volcano, 
          cluster_enrichment_dotplots_stage_1[["BioCarta"]], cols = 2)
m = ggplot(multiplot(union_one_normal_volcano, 
                     cluster_enrichment_dotplots_stage_1[["BioCarta"]], cols = 2))
dev.off(); rm(m)

# Stage 2
tiff("Additional_plots/Stage_2_Volcano_CE_Dotplot.tif", 
     width = 1920, height = 1080, res = 150)
multiplot(union_two_normal_volcano, 
          cluster_enrichment_dotplots_stage_2[["BioCarta"]], cols = 2)
m = ggplot(multiplot(union_two_normal_volcano, 
                     cluster_enrichment_dotplots_stage_2[["BioCarta"]], cols = 2))
dev.off(); rm(m)

# Stage 3
tiff("Additional_plots/Stage_3_Volcano_CE_Dotplot.tif", 
     width = 1920, height = 1080, res = 150)
multiplot(union_three_normal_volcano, 
          cluster_enrichment_dotplots_stage_3[["BioCarta"]], cols = 2)
m = ggplot(multiplot(union_three_normal_volcano, 
                     cluster_enrichment_dotplots_stage_3[["BioCarta"]], cols = 2))
dev.off(); rm(m)

# Stage 4
tiff("Additional_plots/Stage_4_Volcano_CE_Dotplot.tif", 
     width = 1920, height = 1080, res = 150)
multiplot(union_four_normal_volcano, 
          cluster_enrichment_dotplots_stage_4[["BioCarta"]], cols = 2)
m = ggplot(multiplot(union_four_normal_volcano, 
                     cluster_enrichment_dotplots_stage_4[["BioCarta"]], cols = 2))
dev.off(); rm(m)

# DEmiRNAs #####
miRNA_enrichment = read.xlsx("Mienturnet/Stage_1/Stage_1_Mienturnet_input_DEG_list.xlsx")
colnames(miRNA_enrichment) = miRNA_enrichment[1,]
miRNA_enrichment = miRNA_enrichment[-1,]

miRNA_targets = read.xlsx("Mienturnet/Stage_1/Stage_1_Mienturnet_Enrichment_results_miRTarBase.xlsx")
colnames(miRNA_targets) = miRNA_targets[1,]
miRNA_targets = miRNA_targets[-1,]
checklist = miRNA_targets[miRNA_targets$microRNA=="hsa-miR-124-3p", 6:747]
