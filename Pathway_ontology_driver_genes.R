# Pathway and Gene Ontology analysis using pathfindR

# Libraries #####
library(pathfindR)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(cowplot)

# Visualisations of pathways have been moved to a folder outside the repository
# due to large size and difficulties with uploading to GitHub
gene_sets = c("BioCarta", "GO-BP", "GO-CC", "GO-MF", "KEGG", "Reactome")
axis_text_size = c(8, 5, 5, 5, 6, 8)

# Stage 1 #####
# Loading the input to pathfindR (the stage 1 vs normal topTable output):
Stage_1_tT = read.xlsx("DGEA/Union/DGEA_results.xlsx", sheet = 2) %>%
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

cluster_names = c("BioCarta", "GO-CC", "GO-MF", "KEGG")
RNGversion("4.0.2")
set.seed(123)
clustered_results_stage_1 = list()
for (i in 1:length(cluster_names)){
  clustered_results_stage_1[[i]] = cluster_enriched_terms(pathfindR_outputs_stage_1[[cluster_names[[i]]]],
                                                  method = "hierarchical")
}
names(clustered_results_stage_1) = cluster_names

# BioCarta : The maximum average silhouette width was 0.17 for k = 75
# GO-CC    : The maximum average silhouette width was 0.11 for k = 155 
# GO-MF    : The maximum average silhouette width was 0.12 for k = 172
# KEGG     : The maximum average silhouette width was 0.14 for k = 108 

# Wrapping the text of terms with too many characters in their description
wrapped_pathfindR_outputs_stage_1 = pathfindR_outputs_stage_1
for (i in 1:length(wrapped_pathfindR_outputs_stage_1)){
  wrapped_pathfindR_outputs_stage_1[[i]]$Term_Description = stringr::str_wrap(pathfindR_outputs_stage_1[[i]]$Term_Description, 
                                                                      width = 41)
}
rm(i)

wrapped_clustered_pathfindR_outputs_stage_1 = clustered_results_stage_1
for (i in 1:length(wrapped_clustered_pathfindR_outputs_stage_1)){
  wrapped_clustered_pathfindR_outputs_stage_1[[i]]$Term_Description = stringr::str_wrap(clustered_results_stage_1[[i]]$Term_Description, 
                                                                                width = 41)
}
rm(i)

stages = c("Stage 1", "Stage 2", "Stage 3", "Stage 4", "Blood samples")
names(stages) = c("stage_1", "stage_2", "stage_3", "stage_4", "blood")

enrichment_dotplots_stage_1 = list()
cluster_enrichment_dotplots_stage_1 = list()

# Producing dotplots with the results
for (i in 1:length(pathfindR_outputs_stage_1)){
  # unclustered results
  enrichment_dotplots_stage_1[[i]] = enrichment_chart(result_df = wrapped_pathfindR_outputs_stage_1[[i]],
                                                      top_terms = 10)+
    scale_color_gradient(low = "#fca4a4", high = "#fc0303")+
    theme(plot.title = element_text(size = 15, face = "bold", vjust = 1),
          axis.text.y = element_text(color = "black", size = 14),
          axis.text.x = element_text(color = "black", size = 14),
          axis.title.x = element_text(size = 15, face = "bold"))+
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_stage_1)[i],
                        " terms enrichment dotplot - (", stages["stage_1"], ")"))
  tiff(paste0("pathfindR/Stage_1/", names(wrapped_pathfindR_outputs_stage_1)[i], "/",
             names(pathfindR_outputs_stage_1)[i], "_top10_dotplot.tif"), 
       width = 2880, height = 1620, res = 210)
  print(enrichment_dotplots_stage_1[[i]])
  dev.off()
  
  if (names(pathfindR_outputs_stage_1)[i] == "BioCarta" |
      names(pathfindR_outputs_stage_1)[i] == "GO-CC" |
      names(pathfindR_outputs_stage_1)[i] == "GO-MF" |
      names(pathfindR_outputs_stage_1)[i] == "KEGG"){
    # clustered results
    cluster_enrichment_dotplots_stage_1[[i]] = enrichment_chart(result_df = wrapped_clustered_pathfindR_outputs_stage_1[[names(pathfindR_outputs_stage_1)[i]]][clustered_results_stage_1[[names(pathfindR_outputs_stage_1)[i]]]$Status
                                                                                                                    == "Representative", ][1:10,],
                                                                top_terms = NULL,
                            plot_by_cluster = TRUE)+
      scale_color_gradient(low = "#fca4a4", high = "#fc0303")+
      theme(plot.title = element_text(size = 15, face = "bold", vjust = 1),
            axis.text.y = element_text(color = "black", size = 12),
            axis.text.x = element_text(color = "black", size = 14),
            axis.title.x = element_text(size = 15, face = "bold"))+
      labs(title = paste0("Top 10 clustered ", names(wrapped_pathfindR_outputs_stage_1)[i],
                          " terms enrichment dotplot - (", stages["stage_1"], ")"))
    tiff(paste0("pathfindR/Stage_1/", names(wrapped_pathfindR_outputs_stage_1)[i], "/",
               names(pathfindR_outputs_stage_1)[i], "_top10_dotplot_clustered.tif"), 
         width = 2880, height = 1620, res = 210)
    print(cluster_enrichment_dotplots_stage_1[[i]])
    dev.off()
  }
}

names(enrichment_dotplots_stage_1) = names(pathfindR_outputs_stage_1)
names(cluster_enrichment_dotplots_stage_1) = c("BioCarta", "NULL", "GO-CC", "GO-MF", "KEGG")

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

# Defining a legend alignment function 
align_legend <- function(p, hjust = 0.5)
{
  # extract legend
  g <- cowplot::plot_to_gtable(p)
  grobs <- g$grobs
  legend_index <- which(sapply(grobs, function(x) x$name) == "guide-box")
  legend <- grobs[[legend_index]]
  
  # extract guides table
  guides_index <- which(sapply(legend$grobs, function(x) x$name) == "layout")
  
  # there can be multiple guides within one legend box  
  for (gi in guides_index) {
    guides <- legend$grobs[[gi]]
    
    # add extra column for spacing
    # guides$width[5] is the extra spacing from the end of the legend text
    # to the end of the legend title. If we instead distribute it by `hjust:(1-hjust)` on
    # both sides, we get an aligned legend
    spacing <- guides$width[5]
    guides <- gtable::gtable_add_cols(guides, hjust*spacing, 1)
    guides$widths[6] <- (1-hjust)*spacing
    title_index <- guides$layout$name == "title"
    guides$layout$l[title_index] <- 2
    
    # reconstruct guides and write back
    legend$grobs[[gi]] <- guides
  }
  
  # reconstruct legend and write back
  g$grobs[[legend_index]] <- legend
  g
}

term_gene_heatmaps_stage_1 = list()
term_gene_graphs_stage_1 = list()

for (i in 1:length(pathfindR_outputs_stage_1)){
  # term-gene heatmaps
  term_gene_heatmaps_stage_1[[i]] = term_gene_heatmap(result_df = wrapped_pathfindR_outputs_stage_1[[i]],
                                              genes_df = Stage_1_tT,
                                              num_terms = 5,
                                              use_description = TRUE,
                                              low = "darkgreen",
                                              high = "darkred",
                                              mid = "black")+
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0.5),
          axis.text.x = element_text(size = axis_text_size[i], vjust = 0.5, color = "black"),
          axis.text.y = element_text(size = 13),
          legend.title = element_text(size = 13),
          legend.title.align = 0.5,
          legend.direction = "vertical") +
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_stage_1)[i], 
                        " terms - differentially expressed genes heatmap (",
                        stages["stage_1"], ")"),
         fill = "Differential\nExpression\nunits: sd")
  tiff(paste0("pathfindR/Stage_1/", names(wrapped_pathfindR_outputs_stage_1)[i], "/",
             names(wrapped_pathfindR_outputs_stage_1)[i], "_top5_term_gene_heatmap.tif"), 
      width = 3840, height = 648, res = 150)
  print(ggdraw(align_legend(term_gene_heatmaps_stage_1[[i]], hjust = 0.5)))
  dev.off()
  
  # term-gene graphs
  term_gene_graphs_stage_1[[i]] = term_gene_graph(result_df = pathfindR_outputs_stage_1[[i]],
                                              num_terms = 3,
                                              use_description = TRUE,
                                          node_size = "p_val")+
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.subtitle = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 0.5),
          legend.title = element_text(size = 15),
          legend.title.align = 0.5,
          legend.direction = "vertical",
          legend.text = element_text(size = 13))+
    labs(title = paste0("Top 3 ", names(wrapped_pathfindR_outputs_stage_1)[i], 
                        " term - gene graph (",
                        stages["stage_1"], ")"),
         fill = "Differential\nExpression\nunits: sd")
  tiff(paste0("pathfindR/Stage_1/", names(wrapped_pathfindR_outputs_stage_1)[i], "/",
             names(wrapped_pathfindR_outputs_stage_1)[i], "_top3_term_gene_graph.tif"), 
       width = 1920, height = 1080, res = 100)
  print(term_gene_graphs_stage_1[[i]])
  dev.off()
}

# UpSet plots #####
UpSet_plots_stage_1 = list()
for (i in 1:length(pathfindR_outputs_stage_1)){
  # UpSet plot
  UpSet_plots_stage_1[[i]] = UpSet_plot(result_df = wrapped_pathfindR_outputs_stage_1[[i]],
                                genes_df = Stage_1_tT,
                                num_terms = 5,
                                use_description = TRUE,
                                low = "darkgreen",
                                high = "darkred",
                                mid = "black")+
    theme(axis.text.y = element_text(size = axis_text_size[[i]]))
  tiff(paste0("pathfindR/Stage_1/", names(wrapped_pathfindR_outputs_stage_1)[i], "/",
             names(wrapped_pathfindR_outputs_stage_1)[i], "_top5_UpSet_plot.tif"), 
      width = 1444, height = 3840, res = 150)
  print(UpSet_plots_stage_1[[i]])
  dev.off()
}

names(term_gene_graphs_stage_1) = names(pathfindR_outputs_stage_1)
names(term_gene_heatmaps_stage_1) = names(pathfindR_outputs_stage_1)
names(UpSet_plots_stage_1) = names(pathfindR_outputs_stage_1)

# Stage 2 #####
# Loading the input to pathfindR (the stage 2 vs normal topTable output):
Stage_2_tT = read.xlsx("DGEA/Union/DGEA_results.xlsx", sheet = 3) %>%
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

cluster_names = c("BioCarta", "GO-CC", "GO-MF", "KEGG")
RNGversion("4.0.2")
set.seed(123)
clustered_results_stage_2 = list()
for (i in 1:length(cluster_names)){
  clustered_results_stage_2[[i]] = cluster_enriched_terms(pathfindR_outputs_stage_2[[cluster_names[[i]]]],
                                                          method = "hierarchical")
}
names(clustered_results_stage_2) = cluster_names

# BioCarta : The maximum average silhouette width was 0.17 for k = 58
# GO-CC    : The maximum average silhouette width was 0.11 for k = 152 
# GO-MF    : The maximum average silhouette width was 0.09 for k = 168 
# KEGG     : The maximum average silhouette width was 0.13 for k = 78 

# Wrapping the text of terms with too many characters in their description
wrapped_pathfindR_outputs_stage_2 = pathfindR_outputs_stage_2
for (i in 1:length(wrapped_pathfindR_outputs_stage_2)){
  wrapped_pathfindR_outputs_stage_2[[i]]$Term_Description = stringr::str_wrap(pathfindR_outputs_stage_2[[i]]$Term_Description, 
                                                                              width = 41)
}
rm(i)

wrapped_clustered_pathfindR_outputs_stage_2 = clustered_results_stage_2
for (i in 1:length(wrapped_clustered_pathfindR_outputs_stage_2)){
  wrapped_clustered_pathfindR_outputs_stage_2[[i]]$Term_Description = stringr::str_wrap(clustered_results_stage_2[[i]]$Term_Description, 
                                                                                        width = 41)
}
rm(i)

enrichment_dotplots_stage_2 = list()
cluster_enrichment_dotplots_stage_2 = list()

# Producing dotplots with the results
for (i in 1:length(pathfindR_outputs_stage_2)){
  # unclustered results
  enrichment_dotplots_stage_2[[i]] = enrichment_chart(result_df = wrapped_pathfindR_outputs_stage_2[[i]],
                                                      top_terms = 10)+
    scale_color_gradient(low = "#fca4a4", high = "#fc0303")+
    theme(plot.title = element_text(size = 15, face = "bold", vjust = 1),
          axis.text.y = element_text(color = "black", size = 14),
          axis.text.x = element_text(color = "black", size = 14),
          axis.title.x = element_text(size = 15, face = "bold"))+
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_stage_2)[i],
                        " terms enrichment dotplot - (", stages["stage_2"], ")"))
  tiff(paste0("pathfindR/Stage_2/", names(wrapped_pathfindR_outputs_stage_2)[i], "/",
              names(pathfindR_outputs_stage_2)[i], "_top10_dotplot.tif"), 
       width = 2880, height = 1620, res = 210)
  print(enrichment_dotplots_stage_2[[i]])
  dev.off()
  
  if (names(pathfindR_outputs_stage_2)[i] == "BioCarta" |
      names(pathfindR_outputs_stage_2)[i] == "GO-CC" |
      names(pathfindR_outputs_stage_2)[i] == "GO-MF" |
      names(pathfindR_outputs_stage_2)[i] == "KEGG"){
    # clustered results
    cluster_enrichment_dotplots_stage_2[[i]] = enrichment_chart(result_df = wrapped_clustered_pathfindR_outputs_stage_2[[names(pathfindR_outputs_stage_2)[i]]][clustered_results_stage_2[[names(pathfindR_outputs_stage_2)[i]]]$Status
                                                                                                                                                               == "Representative", ][1:10,],
                                                                top_terms = NULL,
                                                                plot_by_cluster = TRUE)+
      scale_color_gradient(low = "#fca4a4", high = "#fc0303")+
      theme(plot.title = element_text(size = 15, face = "bold", vjust = 1),
            axis.text.y = element_text(color = "black", size = 12),
            axis.text.x = element_text(color = "black", size = 14),
            axis.title.x = element_text(size = 15, face = "bold"))+
      labs(title = paste0("Top 10 clustered ", names(wrapped_pathfindR_outputs_stage_2)[i],
                          " terms enrichment dotplot - (", stages["stage_2"], ")"))
    tiff(paste0("pathfindR/Stage_2/", names(wrapped_pathfindR_outputs_stage_2)[i], "/",
                names(pathfindR_outputs_stage_2)[i], "_top10_dotplot_clustered.tif"), 
         width = 2880, height = 1620, res = 210)
    print(cluster_enrichment_dotplots_stage_2[[i]])
    dev.off()
  }
}

names(enrichment_dotplots_stage_2) = names(pathfindR_outputs_stage_2)
names(cluster_enrichment_dotplots_stage_2) = c("BioCarta", "NULL", "GO-CC", "GO-MF", "KEGG")

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

# Loop not executable for the term-gene graph for i = 2 #
for (i in 1:length(pathfindR_outputs_stage_2)){
  # term-gene heatmaps
  term_gene_heatmaps_stage_2[[i]] = term_gene_heatmap(result_df = wrapped_pathfindR_outputs_stage_2[[i]],
                                                      genes_df = Stage_2_tT,
                                                      num_terms = 5,
                                                      use_description = TRUE,
                                                      low = "darkgreen",
                                                      high = "darkred",
                                                      mid = "black")+
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0.5),
          axis.text.x = element_text(size = axis_text_size[i], vjust = 0.5, color = "black"),
          axis.text.y = element_text(size = 13),
          legend.title = element_text(size = 13),
          legend.title.align = 0.5,
          legend.direction = "vertical") +
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_stage_2)[i], 
                        " terms - differentially expressed genes heatmap (",
                        stages["stage_2"], ")"),
         fill = "Differential\nExpression\nunits: sd")
  tiff(paste0("pathfindR/Stage_2/", names(wrapped_pathfindR_outputs_stage_2)[i], "/",
              names(wrapped_pathfindR_outputs_stage_2)[i], "_top5_term_gene_heatmap.tif"), 
       width = 3840, height = 648, res = 150)
  print(ggdraw(align_legend(term_gene_heatmaps_stage_2[[i]], hjust = 0.5)))
  dev.off()

  # term-gene graphs
  term_gene_graphs_stage_2[[i]] = term_gene_graph(result_df = pathfindR_outputs_stage_2[[i]],
                                                  num_terms = 3,
                                                  use_description = TRUE,
                                                  node_size = "p_val")+
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.subtitle = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 0.5),
          legend.title = element_text(size = 15),
          legend.title.align = 0.5,
          legend.direction = "vertical",
          legend.text = element_text(size = 13))+
    labs(title = paste0("Top 3 ", names(wrapped_pathfindR_outputs_stage_2)[i], 
                        " term - gene graph (",
                        stages["stage_2"], ")"),
         fill = "Differential\nExpression\nunits: sd")
  tiff(paste0("pathfindR/Stage_2/", names(wrapped_pathfindR_outputs_stage_2)[i], "/",
              names(wrapped_pathfindR_outputs_stage_2)[i], "_top3_term_gene_graph.tif"), 
       width = 1920, height = 1080, res = 100)
  print(term_gene_graphs_stage_2[[i]])
  dev.off()
}

# UpSet plots #####
UpSet_plots_stage_2 = list()
for (i in 1:length(pathfindR_outputs_stage_2)){
  # UpSet plot
  UpSet_plots_stage_2[[i]] = UpSet_plot(result_df = wrapped_pathfindR_outputs_stage_2[[i]],
                                        genes_df = Stage_2_tT,
                                        num_terms = 5,
                                        use_description = TRUE,
                                        low = "darkgreen",
                                        high = "darkred",
                                        mid = "black")+
    theme(axis.text.y = element_text(size = axis_text_size[[i]]))
  tiff(paste0("pathfindR/Stage_2/", names(wrapped_pathfindR_outputs_stage_2)[i], "/",
              names(wrapped_pathfindR_outputs_stage_2)[i], "_top5_UpSet_plot.tif"), 
       width = 1444, height = 3840, res = 150)
  print(UpSet_plots_stage_2[[i]])
  dev.off()
}

names(term_gene_graphs_stage_2) = names(pathfindR_outputs_stage_2)
names(term_gene_heatmaps_stage_2) = names(pathfindR_outputs_stage_2)
names(UpSet_plots_stage_2) = names(pathfindR_outputs_stage_2)

# Stage 3 #####
# Loading the input to pathfindR (the stage 3 vs normal topTable output):
Stage_3_tT = read.xlsx("DGEA/Union/DGEA_results.xlsx", sheet = 4) %>%
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

cluster_names = c("BioCarta", "GO-CC", "GO-MF", "KEGG")
RNGversion("4.0.2")
set.seed(123)
clustered_results_stage_3 = list()
for (i in 1:length(cluster_names)){
  clustered_results_stage_3[[i]] = cluster_enriched_terms(pathfindR_outputs_stage_3[[cluster_names[[i]]]],
                                                          method = "hierarchical")
}
names(clustered_results_stage_3) = cluster_names

# BioCarta: The maximum average silhouette width was 0.15 for k = 64 
# GO-CC   : The maximum average silhouette width was 0.11 for k = 197
# GO-MF   : The maximum average silhouette width was 0.1  for k = 186
# KEGG    : The maximum average silhouette width was 0.13 for k = 75

# Wrapping the text of terms with too many characters in their description
wrapped_pathfindR_outputs_stage_3 = pathfindR_outputs_stage_3
for (i in 1:length(wrapped_pathfindR_outputs_stage_3)){
  wrapped_pathfindR_outputs_stage_3[[i]]$Term_Description = stringr::str_wrap(pathfindR_outputs_stage_3[[i]]$Term_Description, 
                                                                              width = 41)
}
rm(i)

wrapped_clustered_pathfindR_outputs_stage_3 = clustered_results_stage_3
for (i in 1:length(wrapped_clustered_pathfindR_outputs_stage_3)){
  wrapped_clustered_pathfindR_outputs_stage_3[[i]]$Term_Description = stringr::str_wrap(clustered_results_stage_3[[i]]$Term_Description, 
                                                                                        width = 41)
}
rm(i)

enrichment_dotplots_stage_3 = list()
cluster_enrichment_dotplots_stage_3 = list()

# Producing dotplots with the results
for (i in 1:length(pathfindR_outputs_stage_3)){
  # unclustered results
  enrichment_dotplots_stage_3[[i]] = enrichment_chart(result_df = wrapped_pathfindR_outputs_stage_3[[i]],
                                                      top_terms = 10)+
    scale_color_gradient(low = "#fca4a4", high = "#fc0303")+
    theme(plot.title = element_text(size = 15, face = "bold", vjust = 1),
          axis.text.y = element_text(color = "black", size = 14),
          axis.text.x = element_text(color = "black", size = 14),
          axis.title.x = element_text(size = 15, face = "bold"))+
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_stage_3)[i],
                        " terms enrichment dotplot - (", stages["stage_3"], ")"))
  tiff(paste0("pathfindR/Stage_3/", names(wrapped_pathfindR_outputs_stage_3)[i], "/",
              names(pathfindR_outputs_stage_3)[i], "_top10_dotplot.tif"), 
       width = 2880, height = 1620, res = 210)
  print(enrichment_dotplots_stage_3[[i]])
  dev.off()
  
  if (names(pathfindR_outputs_stage_3)[i] == "BioCarta" |
      names(pathfindR_outputs_stage_3)[i] == "GO-CC" |
      names(pathfindR_outputs_stage_3)[i] == "GO-MF" |
      names(pathfindR_outputs_stage_3)[i] == "KEGG"){
    # clustered results
    cluster_enrichment_dotplots_stage_3[[i]] = enrichment_chart(result_df = wrapped_clustered_pathfindR_outputs_stage_3[[names(pathfindR_outputs_stage_3)[i]]][clustered_results_stage_3[[names(pathfindR_outputs_stage_3)[i]]]$Status
                                                                                                                                                               == "Representative", ][1:10,],
                                                                top_terms = NULL,
                                                                plot_by_cluster = TRUE)+
      scale_color_gradient(low = "#fca4a4", high = "#fc0303")+
      theme(plot.title = element_text(size = 15, face = "bold", vjust = 1),
            axis.text.y = element_text(color = "black", size = 12),
            axis.text.x = element_text(color = "black", size = 14),
            axis.title.x = element_text(size = 15, face = "bold"))+
      labs(title = paste0("Top 10 clustered ", names(wrapped_pathfindR_outputs_stage_3)[i],
                          " terms enrichment dotplot - (", stages["stage_3"], ")"))
    tiff(paste0("pathfindR/Stage_3/", names(wrapped_pathfindR_outputs_stage_3)[i], "/",
                names(pathfindR_outputs_stage_3)[i], "_top10_dotplot_clustered.tif"), 
         width = 2880, height = 1620, res = 210)
    print(cluster_enrichment_dotplots_stage_3[[i]])
    dev.off()
  }
}

names(enrichment_dotplots_stage_3) = names(pathfindR_outputs_stage_3)
names(cluster_enrichment_dotplots_stage_3) = c("BioCarta", "NULL", "GO-CC", "GO-MF", "KEGG")

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
  term_gene_heatmaps_stage_3[[i]] = term_gene_heatmap(result_df = wrapped_pathfindR_outputs_stage_3[[i]],
                                                      genes_df = Stage_3_tT,
                                                      num_terms = 5,
                                                      use_description = TRUE,
                                                      low = "darkgreen",
                                                      high = "darkred",
                                                      mid = "black")+
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0.5),
          axis.text.x = element_text(size = axis_text_size[i], vjust = 0.5, color = "black"),
          axis.text.y = element_text(size = 13),
          legend.title = element_text(size = 13),
          legend.title.align = 0.5,
          legend.direction = "vertical") +
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_stage_3)[i], 
                        " terms - differentially expressed genes heatmap (",
                        stages["stage_3"], ")"),
         fill = "Differential\nExpression\nunits: sd")
  tiff(paste0("pathfindR/Stage_3/", names(wrapped_pathfindR_outputs_stage_3)[i], "/",
              names(wrapped_pathfindR_outputs_stage_3)[i], "_top5_term_gene_heatmap.tif"), 
       width = 3840, height = 648, res = 150)
  print(ggdraw(align_legend(term_gene_heatmaps_stage_3[[i]], hjust = 0.5)))
  dev.off()
  
  # term-gene graphs
  term_gene_graphs_stage_3[[i]] = term_gene_graph(result_df = pathfindR_outputs_stage_3[[i]],
                                                  num_terms = 3,
                                                  use_description = TRUE,
                                                  node_size = "p_val")+
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.subtitle = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 0.5),
          legend.title = element_text(size = 15),
          legend.title.align = 0.5,
          legend.direction = "vertical",
          legend.text = element_text(size = 13))+
    labs(title = paste0("Top 3 ", names(wrapped_pathfindR_outputs_stage_3)[i], 
                        " term - gene graph (",
                        stages["stage_3"], ")"),
         fill = "Differential\nExpression\nunits: sd")
  tiff(paste0("pathfindR/Stage_3/", names(wrapped_pathfindR_outputs_stage_3)[i], "/",
              names(wrapped_pathfindR_outputs_stage_3)[i], "_top3_term_gene_graph.tif"), 
       width = 1920, height = 1080, res = 100)
  print(term_gene_graphs_stage_3[[i]])
  dev.off()
}

# UpSet plots #####
UpSet_plots_stage_3 = list()
for (i in 1:length(pathfindR_outputs_stage_3)){
  # UpSet plot
  UpSet_plots_stage_3[[i]] = UpSet_plot(result_df = wrapped_pathfindR_outputs_stage_3[[i]],
                                        genes_df = Stage_3_tT,
                                        num_terms = 5,
                                        use_description = TRUE,
                                        low = "darkgreen",
                                        high = "darkred",
                                        mid = "black")+
    theme(axis.text.y = element_text(size = axis_text_size[[i]]))
  tiff(paste0("pathfindR/Stage_3/", names(wrapped_pathfindR_outputs_stage_3)[i], "/",
              names(wrapped_pathfindR_outputs_stage_3)[i], "_top5_UpSet_plot.tif"), 
       width = 1444, height = 3840, res = 150)
  print(UpSet_plots_stage_3[[i]])
  dev.off()
}

names(term_gene_graphs_stage_3) = names(pathfindR_outputs_stage_3)
names(term_gene_heatmaps_stage_3) = names(pathfindR_outputs_stage_3)
names(UpSet_plots_stage_3) = names(pathfindR_outputs_stage_3)

# Stage 4 #####
# Loading the input to pathfindR (the stage 4 vs normal topTable output):
Stage_4_tT = read.xlsx("DGEA/Union/DGEA_results.xlsx", sheet = 5) %>%
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

cluster_names = c("BioCarta", "GO-CC", "GO-MF", "KEGG")
RNGversion("4.0.2")
set.seed(123)
clustered_results_stage_4 = list()
for (i in 1:length(cluster_names)){
  clustered_results_stage_4[[i]] = cluster_enriched_terms(pathfindR_outputs_stage_4[[cluster_names[[i]]]],
                                                          method = "hierarchical")
}
names(clustered_results_stage_4) = cluster_names

# BioCarta: The maximum average silhouette width was 0.16 for k = 45 
# GO-CC   : The maximum average silhouette width was 0.12 for k = 144
# GO-MF   : The maximum average silhouette width was 0.11 for k = 188
# KEGG    : The maximum average silhouette width was 0.14 for k = 109

# Wrapping the text of terms with too many characters in their description
wrapped_pathfindR_outputs_stage_4 = pathfindR_outputs_stage_4
for (i in 1:length(wrapped_pathfindR_outputs_stage_4)){
  wrapped_pathfindR_outputs_stage_4[[i]]$Term_Description = stringr::str_wrap(pathfindR_outputs_stage_4[[i]]$Term_Description, 
                                                                              width = 41)
}
rm(i)

wrapped_clustered_pathfindR_outputs_stage_4 = clustered_results_stage_4
for (i in 1:length(wrapped_clustered_pathfindR_outputs_stage_4)){
  wrapped_clustered_pathfindR_outputs_stage_4[[i]]$Term_Description = stringr::str_wrap(clustered_results_stage_4[[i]]$Term_Description, 
                                                                                        width = 41)
}
rm(i)

enrichment_dotplots_stage_4 = list()
cluster_enrichment_dotplots_stage_4 = list()

# Producing dotplots with the results
for (i in 1:length(pathfindR_outputs_stage_4)){
  # unclustered results
  enrichment_dotplots_stage_4[[i]] = enrichment_chart(result_df = wrapped_pathfindR_outputs_stage_4[[i]],
                                                      top_terms = 10)+
    scale_color_gradient(low = "#fca4a4", high = "#fc0303")+
    theme(plot.title = element_text(size = 15, face = "bold", vjust = 1),
          axis.text.y = element_text(color = "black", size = 14),
          axis.text.x = element_text(color = "black", size = 14),
          axis.title.x = element_text(size = 15, face = "bold"))+
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_stage_4)[i],
                        " terms enrichment dotplot - (", stages["stage_4"], ")"))
  tiff(paste0("pathfindR/Stage_4/", names(wrapped_pathfindR_outputs_stage_4)[i], "/",
              names(pathfindR_outputs_stage_4)[i], "_top10_dotplot.tif"), 
       width = 2880, height = 1620, res = 210)
  print(enrichment_dotplots_stage_4[[i]])
  dev.off()
  
  if (names(pathfindR_outputs_stage_4)[i] == "BioCarta" |
      names(pathfindR_outputs_stage_4)[i] == "GO-CC" |
      names(pathfindR_outputs_stage_4)[i] == "GO-MF" |
      names(pathfindR_outputs_stage_4)[i] == "KEGG"){
    # clustered results
    cluster_enrichment_dotplots_stage_4[[i]] = enrichment_chart(result_df = wrapped_clustered_pathfindR_outputs_stage_4[[names(pathfindR_outputs_stage_4)[i]]][clustered_results_stage_4[[names(pathfindR_outputs_stage_4)[i]]]$Status
                                                                                                                                                               == "Representative", ][1:10,],
                                                                top_terms = NULL,
                                                                plot_by_cluster = TRUE)+
      scale_color_gradient(low = "#fca4a4", high = "#fc0303")+
      theme(plot.title = element_text(size = 15, face = "bold", vjust = 1),
            axis.text.y = element_text(color = "black", size = 12),
            axis.text.x = element_text(color = "black", size = 14),
            axis.title.x = element_text(size = 15, face = "bold"))+
      labs(title = paste0("Top 10 clustered ", names(wrapped_pathfindR_outputs_stage_4)[i],
                          " terms enrichment dotplot - (", stages["stage_4"], ")"))
    tiff(paste0("pathfindR/Stage_4/", names(wrapped_pathfindR_outputs_stage_4)[i], "/",
                names(pathfindR_outputs_stage_4)[i], "_top10_dotplot_clustered.tif"), 
         width = 2880, height = 1620, res = 210)
    print(cluster_enrichment_dotplots_stage_4[[i]])
    dev.off()
  }
}

names(enrichment_dotplots_stage_4) = names(pathfindR_outputs_stage_4)
names(cluster_enrichment_dotplots_stage_4) = c("BioCarta", "NULL", "GO-CC", "GO-MF", "KEGG")

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
  term_gene_heatmaps_stage_4[[i]] = term_gene_heatmap(result_df = wrapped_pathfindR_outputs_stage_4[[i]],
                                                      genes_df = Stage_4_tT,
                                                      num_terms = 5,
                                                      use_description = TRUE,
                                                      low = "darkgreen",
                                                      high = "darkred",
                                                      mid = "black")+
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0.5),
          axis.text.x = element_text(size = axis_text_size[i], vjust = 0.5, color = "black"),
          axis.text.y = element_text(size = 13),
          legend.title = element_text(size = 13),
          legend.title.align = 0.5,
          legend.direction = "vertical") +
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_stage_4)[i], 
                        " terms - differentially expressed genes heatmap (",
                        stages["stage_4"], ")"),
         fill = "Differential\nExpression\nunits: sd")
  tiff(paste0("pathfindR/Stage_4/", names(wrapped_pathfindR_outputs_stage_4)[i], "/",
              names(wrapped_pathfindR_outputs_stage_4)[i], "_top5_term_gene_heatmap.tif"), 
       width = 3840, height = 648, res = 150)
  print(ggdraw(align_legend(term_gene_heatmaps_stage_4[[i]], hjust = 0.5)))
  dev.off()
  
  # term-gene graphs
  term_gene_graphs_stage_4[[i]] = term_gene_graph(result_df = pathfindR_outputs_stage_4[[i]],
                                                  num_terms = 3,
                                                  use_description = TRUE,
                                                  node_size = "p_val")+
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.subtitle = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 0.5),
          legend.title = element_text(size = 15),
          legend.title.align = 0.5,
          legend.direction = "vertical",
          legend.text = element_text(size = 13))+
    labs(title = paste0("Top 3 ", names(wrapped_pathfindR_outputs_stage_4)[i], 
                        " term - gene graph (",
                        stages["stage_4"], ")"),
         fill = "Differential\nExpression\nunits: sd")
  tiff(paste0("pathfindR/Stage_4/", names(wrapped_pathfindR_outputs_stage_4)[i], "/",
              names(wrapped_pathfindR_outputs_stage_4)[i], "_top3_term_gene_graph.tif"), 
       width = 1920, height = 1080, res = 100)
  print(term_gene_graphs_stage_4[[i]])
  dev.off()
}

# UpSet plots #####
UpSet_plots_stage_4 = list()
for (i in 1:length(pathfindR_outputs_stage_4)){
  # UpSet plot
  UpSet_plots_stage_4[[i]] = UpSet_plot(result_df = wrapped_pathfindR_outputs_stage_4[[i]],
                                        genes_df = Stage_4_tT,
                                        num_terms = 5,
                                        use_description = TRUE,
                                        low = "darkgreen",
                                        high = "darkred",
                                        mid = "black")+
    theme(axis.text.y = element_text(size = axis_text_size[[i]]))
  tiff(paste0("pathfindR/Stage_4/", names(wrapped_pathfindR_outputs_stage_4)[i], "/",
              names(wrapped_pathfindR_outputs_stage_4)[i], "_top5_UpSet_plot.tif"), 
       width = 1444, height = 3840, res = 150)
  print(UpSet_plots_stage_4[[i]])
  dev.off()
}

names(term_gene_graphs_stage_4) = names(pathfindR_outputs_stage_4)
names(term_gene_heatmaps_stage_4) = names(pathfindR_outputs_stage_4)
names(UpSet_plots_stage_4) = names(pathfindR_outputs_stage_4)

# Blood samples #####
# Loading the input to pathfindR (the blood tumor vs normal topTable output):
Blood_tT = read.xlsx("DGEA/Union/Blood_samples_analysis/Blood_TN_z_DE_topTable.xlsx") %>%
  dplyr::select(Gene.Symbol, logFC, adj.P.Val) %>%
  na.omit()

# Preparing a pathfindR loop for enrichment analysis
dirs_blood = paste0("pathfindR/Blood/", gene_sets)
pathfindR_outputs_blood = list()

RNGversion("4.0.2")
set.seed(123)
for (i in 1:length(dirs_blood)){
  pathfindR_outputs_blood[[i]] = run_pathfindR(Blood_tT, gene_sets = gene_sets[i],
                                               p_val_threshold = 0.05, 
                                               visualize_enriched_terms = TRUE,
                                               output_dir = dirs_blood[i], min_gset_size = 10,
                                               max_gset_size = 300, adj_method = 'fdr',
                                               enrichment_threshold = 0.05,
                                               pin_name_path = 'Biogrid', search_method = 'GR',
                                               grMaxDepth = 1, grSearchDepth = 1,
                                               iterations = 10, n_processes = 10)
  cat(paste0("Done with ", gene_sets[i], "\n"))
}
names(pathfindR_outputs_blood) = gene_sets

# Perform hierarchical clustering on the results (average distance metric)
# Not run for GO-All and Reactome because the algorithm time complexity is O(n^3)

cluster_names = c("BioCarta", "GO-CC", "GO-MF", "KEGG")
RNGversion("4.0.2")
set.seed(123)
clustered_results_blood = list()
for (i in 1:length(cluster_names)){
  clustered_results_blood[[i]] = cluster_enriched_terms(pathfindR_outputs_blood[[cluster_names[[i]]]],
                                                          method = "hierarchical")
}
names(clustered_results_blood) = cluster_names

# BioCarta: The maximum average silhouette width was 0.18 for k = 77 
# GO-CC   : The maximum average silhouette width was 0.14 for k = 135
# GO-MF   : The maximum average silhouette width was 0.12 for k = 135
# KEGG    : The maximum average silhouette width was 0.13 for k = 70

# Wrapping the text of terms with too many characters in their description
wrapped_pathfindR_outputs_blood = pathfindR_outputs_blood
for (i in 1:length(wrapped_pathfindR_outputs_blood)){
  wrapped_pathfindR_outputs_blood[[i]]$Term_Description = stringr::str_wrap(pathfindR_outputs_blood[[i]]$Term_Description, 
                                                                            width = 41)
}
rm(i)

wrapped_clustered_pathfindR_outputs_blood = clustered_results_blood
for (i in 1:length(wrapped_clustered_pathfindR_outputs_blood)){
  wrapped_clustered_pathfindR_outputs_blood[[i]]$Term_Description = stringr::str_wrap(clustered_results_blood[[i]]$Term_Description, 
                                                                                      width = 41)
}
rm(i)

enrichment_dotplots_blood = list()
cluster_enrichment_dotplots_blood = list()

# Producing dotplots with the results
for (i in 1:length(pathfindR_outputs_blood)){
  # unclustered results
  enrichment_dotplots_blood[[i]] = enrichment_chart(result_df = wrapped_pathfindR_outputs_blood[[i]],
                                                    top_terms = 10)+
    scale_color_gradient(low = "#fca4a4", high = "#fc0303")+
    theme(plot.title = element_text(size = 15, face = "bold", vjust = 1),
          axis.text.y = element_text(color = "black", size = 14),
          axis.text.x = element_text(color = "black", size = 14),
          axis.title.x = element_text(size = 15, face = "bold"))+
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_blood)[i],
                        " terms enrichment dotplot - (", stages["blood"], ")"))
  tiff(paste0("pathfindR/Blood/", names(wrapped_pathfindR_outputs_blood)[i], "/",
              names(pathfindR_outputs_blood)[i], "_top10_dotplot.tif"), 
       width = 2880, height = 1620, res = 210)
  print(enrichment_dotplots_blood[[i]])
  dev.off()
  
  if (names(pathfindR_outputs_blood)[i] == "BioCarta" |
      names(pathfindR_outputs_blood)[i] == "GO-CC" |
      names(pathfindR_outputs_blood)[i] == "GO-MF" |
      names(pathfindR_outputs_blood)[i] == "KEGG"){
    # clustered results
    cluster_enrichment_dotplots_blood[[i]] = enrichment_chart(result_df = wrapped_clustered_pathfindR_outputs_blood[[names(pathfindR_outputs_blood)[i]]][clustered_results_blood[[names(pathfindR_outputs_blood)[i]]]$Status
                                                                                                                                                         == "Representative", ][1:10,],
                                                              top_terms = NULL,
                                                              plot_by_cluster = TRUE)+
      scale_color_gradient(low = "#fca4a4", high = "#fc0303")+
      theme(plot.title = element_text(size = 15, face = "bold", vjust = 1),
            axis.text.y = element_text(color = "black", size = 12),
            axis.text.x = element_text(color = "black", size = 14),
            axis.title.x = element_text(size = 15, face = "bold"))+
      labs(title = paste0("Top 10 clustered ", names(wrapped_pathfindR_outputs_blood)[i],
                          " terms enrichment dotplot - (", stages["blood"], ")"))
    tiff(paste0("pathfindR/Blood/", names(wrapped_pathfindR_outputs_blood)[i], "/",
                names(pathfindR_outputs_blood)[i], "_top10_dotplot_clustered.tif"), 
         width = 2880, height = 1620, res = 210)
    print(cluster_enrichment_dotplots_blood[[i]])
    dev.off()
  }
}

names(enrichment_dotplots_blood) = names(pathfindR_outputs_blood)
names(cluster_enrichment_dotplots_blood) = c("BioCarta", "NULL", "GO-CC", "GO-MF", "KEGG")

# Write out results in a comprehensive .xlsx file
wb = createWorkbook()
addWorksheet(wb, "BioCarta")
writeData(wb, "BioCarta", clustered_results_blood[["BioCarta"]])
addWorksheet(wb, "GO-BP")
writeData(wb, "GO-BP", pathfindR_outputs_blood[["GO-BP"]])
addWorksheet(wb, "GO-CC")
writeData(wb, "GO-CC", pathfindR_outputs_blood[["GO-CC"]])
addWorksheet(wb, "GO-MF")
writeData(wb, "GO-MF", pathfindR_outputs_blood[["GO-MF"]])
addWorksheet(wb, "KEGG")
writeData(wb, "KEGG", clustered_results_blood[["KEGG"]])
addWorksheet(wb, "Reactome")
writeData(wb, "Reactome", pathfindR_outputs_blood[["Reactome"]])
saveWorkbook(wb, file = "pathfindR/Blood/Comprehensive_pathfindR_output.xlsx",
             overwrite = TRUE); rm(wb)

# Representative terms output file (BioCarta, GO-CC, GO-MF & KEGG)
wb2 = createWorkbook()
addWorksheet(wb2, "BioCarta - rep")
writeData(wb2, "BioCarta - rep", clustered_results_blood[["BioCarta"]][clustered_results_blood[["BioCarta"]]$Status
                                                                       == "Representative", ])
addWorksheet(wb2, "GOCC - rep")
writeData(wb2, "GOCC - rep", clustered_results_blood[["GO-CC"]][clustered_results_blood[["GO-CC"]]$Status
                                                                == "Representative", ])
addWorksheet(wb2, "GOMF - rep")
writeData(wb2, "GOMF - rep", clustered_results_blood[["GO-MF"]][clustered_results_blood[["GO-MF"]]$Status
                                                                == "Representative", ])
addWorksheet(wb2, "KEGG - rep")
writeData(wb2, "KEGG - rep", clustered_results_blood[["KEGG"]][clustered_results_blood[["KEGG"]]$Status
                                                               == "Representative", ])
saveWorkbook(wb2, file = "pathfindR/Blood/Representative_terms.xlsx",
             overwrite = TRUE); rm(wb2)

# Term-gene heatmaps and term-gene graphs #####
term_gene_heatmaps_blood = list()
term_gene_graphs_blood = list()

for (i in 1:length(pathfindR_outputs_blood)){
  # term-gene heatmaps
  term_gene_heatmaps_blood[[i]] = term_gene_heatmap(result_df = wrapped_pathfindR_outputs_blood[[i]],
                                                    genes_df = Blood_tT,
                                                    num_terms = 5,
                                                    use_description = TRUE,
                                                    low = "darkgreen",
                                                    high = "darkred",
                                                    mid = "black")+
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0.5),
          axis.text.x = element_text(size = rel(1), vjust = 0.5, color = "black"),
          axis.text.y = element_text(size = 13),
          legend.title = element_text(size = 13),
          legend.title.align = 0.5,
          legend.direction = "vertical") +
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_blood)[i], 
                        " terms - differentially expressed genes heatmap (",
                        stages["blood"], ")"),
         fill = "Differential\nExpression\nunits: sd")
  tiff(paste0("pathfindR/Blood/", names(wrapped_pathfindR_outputs_blood)[i], "/",
              names(wrapped_pathfindR_outputs_blood)[i], "_top5_term_gene_heatmap.tif"), 
       width = 3840, height = 648, res = 150)
  print(ggdraw(align_legend(term_gene_heatmaps_blood[[i]], hjust = 0.5)))
  dev.off()
  
  # term-gene graphs
  term_gene_graphs_blood[[i]] = term_gene_graph(result_df = pathfindR_outputs_blood[[i]],
                                                num_terms = 3,
                                                use_description = TRUE,
                                                node_size = "p_val")+
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.subtitle = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 0.5),
          legend.title = element_text(size = 15),
          legend.title.align = 0.5,
          legend.direction = "vertical",
          legend.text = element_text(size = 13))+
    labs(title = paste0("Top 3 ", names(wrapped_pathfindR_outputs_blood)[i], 
                        " term - gene graph (",
                        stages["blood"], ")"),
         fill = "Differential\nExpression\nunits: sd")
  tiff(paste0("pathfindR/Blood/", names(wrapped_pathfindR_outputs_blood)[i], "/",
              names(wrapped_pathfindR_outputs_blood)[i], "_top3_term_gene_graph.tif"), 
       width = 1920, height = 1080, res = 100)
  print(term_gene_graphs_blood[[i]])
  dev.off()
}

# UpSet plots #####
UpSet_plots_blood = list()
for (i in 1:length(pathfindR_outputs_blood)){
  # UpSet plot
  UpSet_plots_blood[[i]] = UpSet_plot(result_df = wrapped_pathfindR_outputs_blood[[i]],
                                      genes_df = Blood_tT,
                                      num_terms = 5,
                                      use_description = TRUE,
                                      low = "darkgreen",
                                      high = "darkred",
                                      mid = "black")+
    theme(axis.text.y = element_text(size = axis_text_size[[i]]))
  tiff(paste0("pathfindR/Blood/", names(wrapped_pathfindR_outputs_blood)[i], "/",
              names(wrapped_pathfindR_outputs_blood)[i], "_top5_UpSet_plot.tif"), 
       width = 1444, height = 3840, res = 150)
  print(UpSet_plots_blood[[i]])
  dev.off()
}

names(term_gene_graphs_blood) = names(pathfindR_outputs_blood)
names(term_gene_heatmaps_blood) = names(pathfindR_outputs_blood)
names(UpSet_plots_blood) = names(pathfindR_outputs_blood)

##### REQUIRES COMPLETE-CASE MATRIX - CANNOT BE RUN IN OUR CASE #####
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

sig_DGEA_stage_1 = read.xlsx("DGEA/Union/DGEA_results.xlsx", sheet = 2) %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::rename(Gene.Symbol_pre = Gene.Symbol)

sig_DGEA_stage_2 = read.xlsx("DGEA/Union/DGEA_results.xlsx", sheet = 3) %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::rename(Gene.Symbol_pre = Gene.Symbol)

sig_DGEA_stage_3 = read.xlsx("DGEA/Union/DGEA_results.xlsx", sheet = 4) %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::rename(Gene.Symbol_pre = Gene.Symbol)

sig_DGEA_stage_4 = read.xlsx("DGEA/Union/DGEA_results.xlsx", sheet = 5) %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::rename(Gene.Symbol_pre = Gene.Symbol)

sig_DGEA_blood = read.xlsx("DGEA/Union/Blood_samples_analysis/Blood_TN_z_DE_topTable.xlsx") %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::rename(Gene.Symbol_pre = Gene.Symbol)

sig_DGEA = list(sig_DGEA_stage_1, sig_DGEA_stage_2, sig_DGEA_stage_3, sig_DGEA_stage_4,
                sig_DGEA_blood)
names(sig_DGEA) = c("Stage_1", "Stage_2", "Stage_3", "Stage_4", "Blood samples")
rm(sig_DGEA_stage_1, sig_DGEA_stage_2, sig_DGEA_stage_3, sig_DGEA_stage_4, sig_DGEA_blood); gc()

# Drivers lists
general_drivers = list()
PDAC_drivers = list()

for (i in 1:length(sig_DGEA)){
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
for (i in 1:length(general_drivers)){
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
library(EnhancedVolcano)

# ggVennDiagram
# Generating a suitable object with stat. sig (p.adj < 0.05) DEGs from all stages
venn_DEG = list(`Stage 1` = sig_DGEA[["Stage_1"]]$EntrezGene.ID[sig_DGEA[["Stage_1"]]$adj.P.Val < 0.05],
            `Stage 2` = sig_DGEA[["Stage_2"]]$EntrezGene.ID[sig_DGEA[["Stage_2"]]$adj.P.Val < 0.05],
            `Stage 3` = sig_DGEA[["Stage_3"]]$EntrezGene.ID[sig_DGEA[["Stage_3"]]$adj.P.Val < 0.05],
            `Stage 4` = sig_DGEA[["Stage_4"]]$EntrezGene.ID[sig_DGEA[["Stage_4"]]$adj.P.Val < 0.05],
            `Blood samples` = sig_DGEA[["Blood samples"]]$EntrezGene.ID[sig_DGEA[["Blood samples"]]$adj.P.Val < 0.05])

diagram1 = ggVennDiagram(
  venn_DEG, label_alpha = 0,
  category.names = names(venn_DEG)
) +
  scale_fill_gradient(low = "white", high = "darkred") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size =18))+
  labs(title = "Venn diagram of significantly differentially expressed genes (DEGs)")
tiff("Additional_plots/Venn/ggVennDiagram_DEG_Venn_tissue_and_blood.tif", 
     width = 1920, height = 1080, res = 130)
diagram1
dev.off()

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

# A Venn diagram just for tumor samples from pancreatic tissue
diagram2 = ggVennDiagram(
  venn_DEG[1:4], label_alpha = 0,
  category.names = names(venn_DEG[1:4])
) +
  scale_fill_gradient(low = "white", high = "darkred") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size  = 18))+
  labs(title = "Venn diagram of sig. diff. expressed genes (DEGs): tumor tissue")
tiff("Additional_plots/Venn/ggVennDiagram_stages_DEG_Venn.tif", 
     width = 1920, height = 1080, res = 200)
diagram2
dev.off()

# DEG venn with blood and tumors multiplot
tiff("Additional_plots/Venn/ggVennDiagram_DEG_multiplot.tif", 
     width = 2880, height = 1620, res = 150)
multiplot(diagram1, diagram2, cols = 2)
m = ggplot(multiplot(diagram1, diagram2, cols = 2))
dev.off(); rm(m, diagram1, diagram2)

# Create a 2x3 Venn diagram multiplot that compares the subnetwork enrichment
# between 4 stages with respect to BioCarta, GO-BP, GO-CC, GO-MF, KEGG, Reactome
# 3 versions: all results, Top 100, clustered results (representative terms)
# All will be produced by ggVennDiagram

# Producing the necessary lists:
venn_BioCarta = list(pathfindR_outputs_stage_1[["BioCarta"]]$Term_Description,
                     pathfindR_outputs_stage_2[["BioCarta"]]$Term_Description,
                     pathfindR_outputs_stage_3[["BioCarta"]]$Term_Description,
                     pathfindR_outputs_stage_4[["BioCarta"]]$Term_Description,
                     pathfindR_outputs_blood[["BioCarta"]]$Term_Description)
names(venn_BioCarta) = names(venn_DEG)

venn_GO_BP = list(pathfindR_outputs_stage_1[["GO-BP"]]$Term_Description,
                  pathfindR_outputs_stage_2[["GO-BP"]]$Term_Description,
                  pathfindR_outputs_stage_3[["GO-BP"]]$Term_Description,
                  pathfindR_outputs_stage_4[["GO-BP"]]$Term_Description,
                  pathfindR_outputs_blood[["GO-BP"]]$Term_Description)
names(venn_GO_BP) = names(venn_DEG)

venn_GO_CC = list(pathfindR_outputs_stage_1[["GO-CC"]]$Term_Description,
                  pathfindR_outputs_stage_2[["GO-CC"]]$Term_Description,
                  pathfindR_outputs_stage_3[["GO-CC"]]$Term_Description,
                  pathfindR_outputs_stage_4[["GO-CC"]]$Term_Description,
                  pathfindR_outputs_blood[["GO-CC"]]$Term_Description)
names(venn_GO_CC) = names(venn_DEG)

venn_GO_MF = list(pathfindR_outputs_stage_1[["GO-MF"]]$Term_Description,
                  pathfindR_outputs_stage_2[["GO-MF"]]$Term_Description,
                  pathfindR_outputs_stage_3[["GO-MF"]]$Term_Description,
                  pathfindR_outputs_stage_4[["GO-MF"]]$Term_Description,
                  pathfindR_outputs_blood[["GO-MF"]]$Term_Description)
names(venn_GO_MF) = names(venn_DEG)

venn_KEGG = list(pathfindR_outputs_stage_1[["KEGG"]]$Term_Description,
                  pathfindR_outputs_stage_2[["KEGG"]]$Term_Description,
                  pathfindR_outputs_stage_3[["KEGG"]]$Term_Description,
                  pathfindR_outputs_stage_4[["KEGG"]]$Term_Description,
                 pathfindR_outputs_blood[["KEGG"]]$Term_Description)
names(venn_KEGG) = names(venn_DEG)

venn_Reactome = list(pathfindR_outputs_stage_1[["Reactome"]]$Term_Description,
                 pathfindR_outputs_stage_2[["Reactome"]]$Term_Description,
                 pathfindR_outputs_stage_3[["Reactome"]]$Term_Description,
                 pathfindR_outputs_stage_4[["Reactome"]]$Term_Description,
                 pathfindR_outputs_blood[["Reactome"]]$Term_Description)
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
                     pathfindR_outputs_stage_4[["BioCarta"]]$Term_Description[1:100],
                     pathfindR_outputs_blood[["BioCarta"]]$Term_Description[1:100])
names(top100_venn_BioCarta) = names(venn_DEG)

top100_venn_GO_BP = list(pathfindR_outputs_stage_1[["GO-BP"]]$Term_Description[1:100],
                  pathfindR_outputs_stage_2[["GO-BP"]]$Term_Description[1:100],
                  pathfindR_outputs_stage_3[["GO-BP"]]$Term_Description[1:100],
                  pathfindR_outputs_stage_4[["GO-BP"]]$Term_Description[1:100],
                  pathfindR_outputs_blood[["GO-BP"]]$Term_Description[1:100])
names(top100_venn_GO_BP) = names(venn_DEG)

top100_venn_GO_CC = list(pathfindR_outputs_stage_1[["GO-CC"]]$Term_Description[1:100],
                  pathfindR_outputs_stage_2[["GO-CC"]]$Term_Description[1:100],
                  pathfindR_outputs_stage_3[["GO-CC"]]$Term_Description[1:100],
                  pathfindR_outputs_stage_4[["GO-CC"]]$Term_Description[1:100],
                  pathfindR_outputs_blood[["GO-CC"]]$Term_Description[1:100])
names(top100_venn_GO_CC) = names(venn_DEG)

top100_venn_GO_MF = list(pathfindR_outputs_stage_1[["GO-MF"]]$Term_Description[1:100],
                  pathfindR_outputs_stage_2[["GO-MF"]]$Term_Description[1:100],
                  pathfindR_outputs_stage_3[["GO-MF"]]$Term_Description[1:100],
                  pathfindR_outputs_stage_4[["GO-MF"]]$Term_Description[1:100],
                  pathfindR_outputs_blood[["GO-MF"]]$Term_Description[1:100])
names(top100_venn_GO_MF) = names(venn_DEG)

top100_venn_KEGG = list(pathfindR_outputs_stage_1[["KEGG"]]$Term_Description[1:100],
                 pathfindR_outputs_stage_2[["KEGG"]]$Term_Description[1:100],
                 pathfindR_outputs_stage_3[["KEGG"]]$Term_Description[1:100],
                 pathfindR_outputs_stage_4[["KEGG"]]$Term_Description[1:100],
                 pathfindR_outputs_blood[["KEGG"]]$Term_Description[1:100])
names(top100_venn_KEGG) = names(venn_DEG)

top100_venn_Reactome = list(pathfindR_outputs_stage_1[["Reactome"]]$Term_Description[1:100],
                     pathfindR_outputs_stage_2[["Reactome"]]$Term_Description[1:100],
                     pathfindR_outputs_stage_3[["Reactome"]]$Term_Description[1:100],
                     pathfindR_outputs_stage_4[["Reactome"]]$Term_Description[1:100],
                     pathfindR_outputs_blood[["Reactome"]]$Term_Description[1:100])
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
                               clustered_results_stage_4[["BioCarta"]]$Term_Description[clustered_results_stage_4[["BioCarta"]]$Status == "Representative"],
                               clustered_results_blood[["BioCarta"]]$Term_Description[clustered_results_blood[["BioCarta"]]$Status == "Representative"])
names(clustered_venn_BioCarta) = names(venn_DEG)

clustered_venn_GO_CC = list(clustered_results_stage_1[["GO-CC"]]$Term_Description[clustered_results_stage_1[["GO-CC"]]$Status == "Representative"],
                            clustered_results_stage_2[["GO-CC"]]$Term_Description[clustered_results_stage_2[["GO-CC"]]$Status == "Representative"],
                            clustered_results_stage_3[["GO-CC"]]$Term_Description[clustered_results_stage_3[["GO-CC"]]$Status == "Representative"],
                            clustered_results_stage_4[["GO-CC"]]$Term_Description[clustered_results_stage_4[["GO-CC"]]$Status == "Representative"],
                            clustered_results_blood[["GO-CC"]]$Term_Description[clustered_results_blood[["GO-CC"]]$Status == "Representative"])
names(clustered_venn_GO_CC) = names(venn_DEG)

clustered_venn_GO_MF = list(clustered_results_stage_1[["GO-MF"]]$Term_Description[clustered_results_stage_1[["GO-MF"]]$Status == "Representative"],
                            clustered_results_stage_2[["GO-MF"]]$Term_Description[clustered_results_stage_2[["GO-MF"]]$Status == "Representative"],
                            clustered_results_stage_3[["GO-MF"]]$Term_Description[clustered_results_stage_3[["GO-MF"]]$Status == "Representative"],
                            clustered_results_stage_4[["GO-MF"]]$Term_Description[clustered_results_stage_4[["GO-MF"]]$Status == "Representative"],
                            clustered_results_blood[["GO-MF"]]$Term_Description[clustered_results_blood[["GO-MF"]]$Status == "Representative"])
names(clustered_venn_GO_MF) = names(venn_DEG)

clustered_venn_KEGG = list(clustered_results_stage_1[["KEGG"]]$Term_Description[clustered_results_stage_1[["KEGG"]]$Status == "Representative"],
                           clustered_results_stage_2[["KEGG"]]$Term_Description[clustered_results_stage_2[["KEGG"]]$Status == "Representative"],
                           clustered_results_stage_3[["KEGG"]]$Term_Description[clustered_results_stage_3[["KEGG"]]$Status == "Representative"],
                           clustered_results_stage_4[["KEGG"]]$Term_Description[clustered_results_stage_4[["KEGG"]]$Status == "Representative"],
                           clustered_results_blood[["KEGG"]]$Term_Description[clustered_results_blood[["KEGG"]]$Status == "Representative"])
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
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18))+
    labs(title = paste0("Venn diagram of ", names(venn_path)[i], " enriched networks"))
}

tiff("Additional_plots/Venn/all_path_results_ggVennDiagram_multiplot.tif", 
     width = 3840, height = 2160, res = 135)
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
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18))+
    labs(title = paste0("Venn diagram of top 100 ", names(top100_venn_path)[i], 
                        " enriched networks"))
}

tiff("Additional_plots/Venn/top100_path_results_ggVennDiagram_multiplot.tif", 
     width = 3840, height = 2160, res = 135)
multiplot(top100_results_path_venn[[1]], top100_results_path_venn[[2]],
          top100_results_path_venn[[3]], top100_results_path_venn[[4]],
          top100_results_path_venn[[5]], top100_results_path_venn[[6]], cols = 3)
m = ggplot(multiplot(top100_results_path_venn[[1]], top100_results_path_venn[[2]],
                     top100_results_path_venn[[3]], top100_results_path_venn[[4]],
                     top100_results_path_venn[[5]], top100_results_path_venn[[6]], cols = 3))
dev.off(); rm(m)

# clustered enriched networks multiplot (not for GO-BP, Reactome)
for (i in 1:length(clustered_venn_path)){
  clustered_results_path_venn[[i]] = ggVennDiagram(
    clustered_venn_path[[i]], label_alpha = 0,
    category.names = names(clustered_venn_path[[i]])
  ) +
    scale_fill_gradient(low = "white", high = "darkred") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18))+
    labs(title = paste0("Venn diagram of clustered ", names(clustered_venn_path)[i], 
                        " enriched networks"))
}

tiff("Additional_plots/Venn/clustered_path_results_ggVennDiagram_multiplot.tif", 
     width = 2560, height = 2160, res = 140)
multiplot(clustered_results_path_venn[[1]], clustered_results_path_venn[[2]],
          clustered_results_path_venn[[3]], clustered_results_path_venn[[4]], cols = 2)
m = ggplot(multiplot(clustered_results_path_venn[[1]], clustered_results_path_venn[[2]],
                     clustered_results_path_venn[[3]], clustered_results_path_venn[[4]], cols = 2))
dev.off(); rm(m)

# Venn mapping data frames
# BioCarta
BioCarta_map = as.data.frame(matrix(nrow=length(unique(c(pathfindR_outputs_blood$BioCarta$Term_Description, 
                                                         pathfindR_outputs_stage_1$BioCarta$Term_Description,
                                                         pathfindR_outputs_stage_2$BioCarta$Term_Description, 
                                                         pathfindR_outputs_stage_3$BioCarta$Term_Description,
                                                         pathfindR_outputs_stage_4$BioCarta$Term_Description))), ncol = 6))
colnames(BioCarta_map) = c("Term", "Stage_1", "Stage_2", "Stage_3", "Stage_4", "Blood")
BioCarta_map$Term = sort(unique(c(pathfindR_outputs_blood$BioCarta$Term_Description, 
                                  pathfindR_outputs_stage_1$BioCarta$Term_Description,
                                  pathfindR_outputs_stage_2$BioCarta$Term_Description, 
                                  pathfindR_outputs_stage_3$BioCarta$Term_Description,
                                  pathfindR_outputs_stage_4$BioCarta$Term_Description)))
BioCarta_map$Stage_1 = ifelse(BioCarta_map$Term %in% pathfindR_outputs_stage_1$BioCarta$Term_Description,
                              "Present", "Absent")
BioCarta_map$Stage_2 = ifelse(BioCarta_map$Term %in% pathfindR_outputs_stage_2$BioCarta$Term_Description,
                              "Present", "Absent")
BioCarta_map$Stage_3 = ifelse(BioCarta_map$Term %in% pathfindR_outputs_stage_3$BioCarta$Term_Description,
                              "Present", "Absent")
BioCarta_map$Stage_4 = ifelse(BioCarta_map$Term %in% pathfindR_outputs_stage_4$BioCarta$Term_Description,
                              "Present", "Absent")
BioCarta_map$Blood = ifelse(BioCarta_map$Term %in% pathfindR_outputs_blood$BioCarta$Term_Description,
                            "Present", "Absent")
BioCarta_map[which(BioCarta_map$Stage_1 == "Present" & BioCarta_map$Stage_2 == "Absent" &
                     BioCarta_map$Stage_3 == "Absent" & BioCarta_map$Stage_4 == "Absent" &
                     BioCarta_map$Blood == "Present"),]

# GO-BP
GO_BP_map = as.data.frame(matrix(nrow=length(unique(c(pathfindR_outputs_blood$`GO-BP`$Term_Description, 
                                                      pathfindR_outputs_stage_1$`GO-BP`$Term_Description,
                                                      pathfindR_outputs_stage_2$`GO-BP`$Term_Description, 
                                                      pathfindR_outputs_stage_3$`GO-BP`$Term_Description,
                                                      pathfindR_outputs_stage_4$`GO-BP`$Term_Description))), ncol = 6))
colnames(GO_BP_map) = c("Term", "Stage_1", "Stage_2", "Stage_3", "Stage_4", "Blood")
GO_BP_map$Term = sort(unique(c(pathfindR_outputs_blood$`GO-BP`$Term_Description, 
                               pathfindR_outputs_stage_1$`GO-BP`$Term_Description,
                               pathfindR_outputs_stage_2$`GO-BP`$Term_Description, 
                               pathfindR_outputs_stage_3$`GO-BP`$Term_Description,
                               pathfindR_outputs_stage_4$`GO-BP`$Term_Description)))
GO_BP_map$Stage_1 = ifelse(GO_BP_map$Term %in% pathfindR_outputs_stage_1$`GO-BP`$Term_Description,
                           "Present", "Absent")
GO_BP_map$Stage_2 = ifelse(GO_BP_map$Term %in% pathfindR_outputs_stage_2$`GO-BP`$Term_Description,
                           "Present", "Absent")
GO_BP_map$Stage_3 = ifelse(GO_BP_map$Term %in% pathfindR_outputs_stage_3$`GO-BP`$Term_Description,
                           "Present", "Absent")
GO_BP_map$Stage_4 = ifelse(GO_BP_map$Term %in% pathfindR_outputs_stage_4$`GO-BP`$Term_Description,
                           "Present", "Absent")
GO_BP_map$Blood = ifelse(GO_BP_map$Term %in% pathfindR_outputs_blood$`GO-BP`$Term_Description,
                         "Present", "Absent")
GO_BP_map[which(GO_BP_map$Stage_1 == "Present" & GO_BP_map$Stage_2 == "Absent" &
                  GO_BP_map$Stage_3 == "Absent" & GO_BP_map$Stage_4 == "Absent" &
                  GO_BP_map$Blood == "Present"),]

# GO-CC
GO_CC_map = as.data.frame(matrix(nrow=length(unique(c(pathfindR_outputs_blood$`GO-CC`$Term_Description, 
                                                      pathfindR_outputs_stage_1$`GO-CC`$Term_Description,
                                                      pathfindR_outputs_stage_2$`GO-CC`$Term_Description, 
                                                      pathfindR_outputs_stage_3$`GO-CC`$Term_Description,
                                                      pathfindR_outputs_stage_4$`GO-CC`$Term_Description))), ncol = 6))
colnames(GO_CC_map) = c("Term", "Stage_1", "Stage_2", "Stage_3", "Stage_4", "Blood")
GO_CC_map$Term = sort(unique(c(pathfindR_outputs_blood$`GO-CC`$Term_Description, 
                               pathfindR_outputs_stage_1$`GO-CC`$Term_Description,
                               pathfindR_outputs_stage_2$`GO-CC`$Term_Description, 
                               pathfindR_outputs_stage_3$`GO-CC`$Term_Description,
                               pathfindR_outputs_stage_4$`GO-CC`$Term_Description)))
GO_CC_map$Stage_1 = ifelse(GO_CC_map$Term %in% pathfindR_outputs_stage_1$`GO-CC`$Term_Description,
                           "Present", "Absent")
GO_CC_map$Stage_2 = ifelse(GO_CC_map$Term %in% pathfindR_outputs_stage_2$`GO-CC`$Term_Description,
                           "Present", "Absent")
GO_CC_map$Stage_3 = ifelse(GO_CC_map$Term %in% pathfindR_outputs_stage_3$`GO-CC`$Term_Description,
                           "Present", "Absent")
GO_CC_map$Stage_4 = ifelse(GO_CC_map$Term %in% pathfindR_outputs_stage_4$`GO-CC`$Term_Description,
                           "Present", "Absent")
GO_CC_map$Blood = ifelse(GO_CC_map$Term %in% pathfindR_outputs_blood$`GO-CC`$Term_Description,
                         "Present", "Absent")
GO_CC_map[which(GO_CC_map$Stage_1 == "Present" & GO_CC_map$Stage_2 == "Absent" &
                  GO_CC_map$Stage_3 == "Absent" & GO_CC_map$Stage_4 == "Absent" &
                  GO_CC_map$Blood == "Present"),]

# GO-MF
GO_MF_map = as.data.frame(matrix(nrow=length(unique(c(pathfindR_outputs_blood$`GO-MF`$Term_Description, 
                                                      pathfindR_outputs_stage_1$`GO-MF`$Term_Description,
                                                      pathfindR_outputs_stage_2$`GO-MF`$Term_Description, 
                                                      pathfindR_outputs_stage_3$`GO-MF`$Term_Description,
                                                      pathfindR_outputs_stage_4$`GO-MF`$Term_Description))), ncol = 6))
colnames(GO_MF_map) = c("Term", "Stage_1", "Stage_2", "Stage_3", "Stage_4", "Blood")
GO_MF_map$Term = sort(unique(c(pathfindR_outputs_blood$`GO-MF`$Term_Description, 
                               pathfindR_outputs_stage_1$`GO-MF`$Term_Description,
                               pathfindR_outputs_stage_2$`GO-MF`$Term_Description, 
                               pathfindR_outputs_stage_3$`GO-MF`$Term_Description,
                               pathfindR_outputs_stage_4$`GO-MF`$Term_Description)))
GO_MF_map$Stage_1 = ifelse(GO_MF_map$Term %in% pathfindR_outputs_stage_1$`GO-MF`$Term_Description,
                           "Present", "Absent")
GO_MF_map$Stage_2 = ifelse(GO_MF_map$Term %in% pathfindR_outputs_stage_2$`GO-MF`$Term_Description,
                           "Present", "Absent")
GO_MF_map$Stage_3 = ifelse(GO_MF_map$Term %in% pathfindR_outputs_stage_3$`GO-MF`$Term_Description,
                           "Present", "Absent")
GO_MF_map$Stage_4 = ifelse(GO_MF_map$Term %in% pathfindR_outputs_stage_4$`GO-MF`$Term_Description,
                           "Present", "Absent")
GO_MF_map$Blood = ifelse(GO_MF_map$Term %in% pathfindR_outputs_blood$`GO-MF`$Term_Description,
                         "Present", "Absent")
GO_MF_map[which(GO_MF_map$Stage_1 == "Present" & GO_MF_map$Stage_2 == "Absent" &
                  GO_MF_map$Stage_3 == "Absent" & GO_MF_map$Stage_4 == "Absent" &
                  GO_MF_map$Blood == "Present"),]

# KEGG
KEGG_map = as.data.frame(matrix(nrow=length(unique(c(pathfindR_outputs_blood$KEGG$Term_Description, 
                                                     pathfindR_outputs_stage_1$KEGG$Term_Description,
                                                     pathfindR_outputs_stage_2$KEGG$Term_Description, 
                                                     pathfindR_outputs_stage_3$KEGG$Term_Description,
                                                     pathfindR_outputs_stage_4$KEGG$Term_Description))), ncol = 6))
colnames(KEGG_map) = c("Term", "Stage_1", "Stage_2", "Stage_3", "Stage_4", "Blood")
KEGG_map$Term = sort(unique(c(pathfindR_outputs_blood$KEGG$Term_Description, 
                              pathfindR_outputs_stage_1$KEGG$Term_Description,
                              pathfindR_outputs_stage_2$KEGG$Term_Description, 
                              pathfindR_outputs_stage_3$KEGG$Term_Description,
                              pathfindR_outputs_stage_4$KEGG$Term_Description)))
KEGG_map$Stage_1 = ifelse(KEGG_map$Term %in% pathfindR_outputs_stage_1$KEGG$Term_Description,
                          "Present", "Absent")
KEGG_map$Stage_2 = ifelse(KEGG_map$Term %in% pathfindR_outputs_stage_2$KEGG$Term_Description,
                          "Present", "Absent")
KEGG_map$Stage_3 = ifelse(KEGG_map$Term %in% pathfindR_outputs_stage_3$KEGG$Term_Description,
                          "Present", "Absent")
KEGG_map$Stage_4 = ifelse(KEGG_map$Term %in% pathfindR_outputs_stage_4$KEGG$Term_Description,
                          "Present", "Absent")
KEGG_map$Blood = ifelse(KEGG_map$Term %in% pathfindR_outputs_blood$KEGG$Term_Description,
                        "Present", "Absent")
KEGG_map[which(KEGG_map$Stage_1 == "Present" & KEGG_map$Stage_2 == "Absent" &
                 KEGG_map$Stage_3 == "Absent" & KEGG_map$Stage_4 == "Absent" &
                 KEGG_map$Blood == "Present"),]

# Reactome
Reactome_map = as.data.frame(matrix(nrow=length(unique(c(pathfindR_outputs_blood$Reactome$Term_Description, 
                                                         pathfindR_outputs_stage_1$Reactome$Term_Description,
                                                         pathfindR_outputs_stage_2$Reactome$Term_Description, 
                                                         pathfindR_outputs_stage_3$Reactome$Term_Description,
                                                         pathfindR_outputs_stage_4$Reactome$Term_Description))), ncol = 6))
colnames(Reactome_map) = c("Term", "Stage_1", "Stage_2", "Stage_3", "Stage_4", "Blood")
Reactome_map$Term = sort(unique(c(pathfindR_outputs_blood$Reactome$Term_Description, 
                                  pathfindR_outputs_stage_1$Reactome$Term_Description,
                                  pathfindR_outputs_stage_2$Reactome$Term_Description, 
                                  pathfindR_outputs_stage_3$Reactome$Term_Description,
                                  pathfindR_outputs_stage_4$Reactome$Term_Description)))
Reactome_map$Stage_1 = ifelse(Reactome_map$Term %in% pathfindR_outputs_stage_1$Reactome$Term_Description,
                              "Present", "Absent")
Reactome_map$Stage_2 = ifelse(Reactome_map$Term %in% pathfindR_outputs_stage_2$Reactome$Term_Description,
                              "Present", "Absent")
Reactome_map$Stage_3 = ifelse(Reactome_map$Term %in% pathfindR_outputs_stage_3$Reactome$Term_Description,
                              "Present", "Absent")
Reactome_map$Stage_4 = ifelse(Reactome_map$Term %in% pathfindR_outputs_stage_4$Reactome$Term_Description,
                              "Present", "Absent")
Reactome_map$Blood = ifelse(Reactome_map$Term %in% pathfindR_outputs_blood$Reactome$Term_Description,
                            "Present", "Absent")
Reactome_map[which(Reactome_map$Stage_1 == "Present" & Reactome_map$Stage_2 == "Absent" &
                     Reactome_map$Stage_3 == "Absent" & Reactome_map$Stage_4 == "Absent" &
                     Reactome_map$Blood == "Present"),]
# Top100
# BioCarta
Top100_BioCarta_map = as.data.frame(matrix(nrow=length(unique(c(pathfindR_outputs_blood$BioCarta$Term_Description[1:100], 
                                                                pathfindR_outputs_stage_1$BioCarta$Term_Description[1:100],
                                                                pathfindR_outputs_stage_2$BioCarta$Term_Description[1:100], 
                                                                pathfindR_outputs_stage_3$BioCarta$Term_Description[1:100],
                                                                pathfindR_outputs_stage_4$BioCarta$Term_Description[1:100]))), ncol = 6))
colnames(Top100_BioCarta_map) = c("Term", "Stage_1", "Stage_2", "Stage_3", "Stage_4", "Blood")
Top100_BioCarta_map$Term = sort(unique(c(pathfindR_outputs_blood$BioCarta$Term_Description[1:100], 
                                         pathfindR_outputs_stage_1$BioCarta$Term_Description[1:100],
                                         pathfindR_outputs_stage_2$BioCarta$Term_Description[1:100], 
                                         pathfindR_outputs_stage_3$BioCarta$Term_Description[1:100],
                                         pathfindR_outputs_stage_4$BioCarta$Term_Description[1:100])))
Top100_BioCarta_map$Stage_1 = ifelse(Top100_BioCarta_map$Term %in% pathfindR_outputs_stage_1$BioCarta$Term_Description[1:100],
                                     "Present", "Absent")
Top100_BioCarta_map$Stage_2 = ifelse(Top100_BioCarta_map$Term %in% pathfindR_outputs_stage_2$BioCarta$Term_Description[1:100],
                                     "Present", "Absent")
Top100_BioCarta_map$Stage_3 = ifelse(Top100_BioCarta_map$Term %in% pathfindR_outputs_stage_3$BioCarta$Term_Description[1:100],
                                     "Present", "Absent")
Top100_BioCarta_map$Stage_4 = ifelse(Top100_BioCarta_map$Term %in% pathfindR_outputs_stage_4$BioCarta$Term_Description[1:100],
                                     "Present", "Absent")
Top100_BioCarta_map$Blood = ifelse(Top100_BioCarta_map$Term %in% pathfindR_outputs_blood$BioCarta$Term_Description[1:100],
                                   "Present", "Absent")
Top100_BioCarta_map[which(Top100_BioCarta_map$Stage_1 == "Present" & Top100_BioCarta_map$Stage_2 == "Absent" &
                            Top100_BioCarta_map$Stage_3 == "Absent" & Top100_BioCarta_map$Stage_4 == "Absent" &
                            Top100_BioCarta_map$Blood == "Present"),]

# GO-BP
Top100_GO_BP_map = as.data.frame(matrix(nrow=length(unique(c(pathfindR_outputs_blood$`GO-BP`$Term_Description[1:100], 
                                                             pathfindR_outputs_stage_1$`GO-BP`$Term_Description[1:100],
                                                             pathfindR_outputs_stage_2$`GO-BP`$Term_Description[1:100], 
                                                             pathfindR_outputs_stage_3$`GO-BP`$Term_Description[1:100],
                                                             pathfindR_outputs_stage_4$`GO-BP`$Term_Description[1:100]))), ncol = 6))
colnames(Top100_GO_BP_map) = c("Term", "Stage_1", "Stage_2", "Stage_3", "Stage_4", "Blood")
Top100_GO_BP_map$Term = sort(unique(c(pathfindR_outputs_blood$`GO-BP`$Term_Description[1:100], 
                                      pathfindR_outputs_stage_1$`GO-BP`$Term_Description[1:100],
                                      pathfindR_outputs_stage_2$`GO-BP`$Term_Description[1:100], 
                                      pathfindR_outputs_stage_3$`GO-BP`$Term_Description[1:100],
                                      pathfindR_outputs_stage_4$`GO-BP`$Term_Description[1:100])))
Top100_GO_BP_map$Stage_1 = ifelse(Top100_GO_BP_map$Term %in% pathfindR_outputs_stage_1$`GO-BP`$Term_Description[1:100],
                                  "Present", "Absent")
Top100_GO_BP_map$Stage_2 = ifelse(Top100_GO_BP_map$Term %in% pathfindR_outputs_stage_2$`GO-BP`$Term_Description[1:100],
                                  "Present", "Absent")
Top100_GO_BP_map$Stage_3 = ifelse(Top100_GO_BP_map$Term %in% pathfindR_outputs_stage_3$`GO-BP`$Term_Description[1:100],
                                  "Present", "Absent")
Top100_GO_BP_map$Stage_4 = ifelse(Top100_GO_BP_map$Term %in% pathfindR_outputs_stage_4$`GO-BP`$Term_Description[1:100],
                                  "Present", "Absent")
Top100_GO_BP_map$Blood = ifelse(Top100_GO_BP_map$Term %in% pathfindR_outputs_blood$`GO-BP`$Term_Description[1:100],
                                "Present", "Absent")
Top100_GO_BP_map[which(Top100_GO_BP_map$Stage_1 == "Present" & Top100_GO_BP_map$Stage_2 == "Absent" &
                         Top100_GO_BP_map$Stage_3 == "Absent" & Top100_GO_BP_map$Stage_4 == "Absent" &
                         Top100_GO_BP_map$Blood == "Present"),]

# GO-CC
Top100_GO_CC_map = as.data.frame(matrix(nrow=length(unique(c(pathfindR_outputs_blood$`GO-CC`$Term_Description[1:100], 
                                                             pathfindR_outputs_stage_1$`GO-CC`$Term_Description[1:100],
                                                             pathfindR_outputs_stage_2$`GO-CC`$Term_Description[1:100], 
                                                             pathfindR_outputs_stage_3$`GO-CC`$Term_Description[1:100],
                                                             pathfindR_outputs_stage_4$`GO-CC`$Term_Description[1:100]))), ncol = 6))
colnames(Top100_GO_CC_map) = c("Term", "Stage_1", "Stage_2", "Stage_3", "Stage_4", "Blood")
Top100_GO_CC_map$Term = sort(unique(c(pathfindR_outputs_blood$`GO-CC`$Term_Description[1:100], 
                                      pathfindR_outputs_stage_1$`GO-CC`$Term_Description[1:100],
                                      pathfindR_outputs_stage_2$`GO-CC`$Term_Description[1:100], 
                                      pathfindR_outputs_stage_3$`GO-CC`$Term_Description[1:100],
                                      pathfindR_outputs_stage_4$`GO-CC`$Term_Description[1:100])))
Top100_GO_CC_map$Stage_1 = ifelse(Top100_GO_CC_map$Term %in% pathfindR_outputs_stage_1$`GO-CC`$Term_Description[1:100],
                                  "Present", "Absent")
Top100_GO_CC_map$Stage_2 = ifelse(Top100_GO_CC_map$Term %in% pathfindR_outputs_stage_2$`GO-CC`$Term_Description[1:100],
                                  "Present", "Absent")
Top100_GO_CC_map$Stage_3 = ifelse(Top100_GO_CC_map$Term %in% pathfindR_outputs_stage_3$`GO-CC`$Term_Description[1:100],
                                  "Present", "Absent")
Top100_GO_CC_map$Stage_4 = ifelse(Top100_GO_CC_map$Term %in% pathfindR_outputs_stage_4$`GO-CC`$Term_Description[1:100],
                                  "Present", "Absent")
Top100_GO_CC_map$Blood = ifelse(Top100_GO_CC_map$Term %in% pathfindR_outputs_blood$`GO-CC`$Term_Description[1:100],
                                "Present", "Absent")
Top100_GO_CC_map[which(Top100_GO_CC_map$Stage_1 == "Present" & Top100_GO_CC_map$Stage_2 == "Absent" &
                         Top100_GO_CC_map$Stage_3 == "Absent" & Top100_GO_CC_map$Stage_4 == "Absent" &
                         Top100_GO_CC_map$Blood == "Present"),]

# GO-MF
Top100_GO_MF_map = as.data.frame(matrix(nrow=length(unique(c(pathfindR_outputs_blood$`GO-MF`$Term_Description[1:100], 
                                                             pathfindR_outputs_stage_1$`GO-MF`$Term_Description[1:100],
                                                             pathfindR_outputs_stage_2$`GO-MF`$Term_Description[1:100], 
                                                             pathfindR_outputs_stage_3$`GO-MF`$Term_Description[1:100],
                                                             pathfindR_outputs_stage_4$`GO-MF`$Term_Description[1:100]))), ncol = 6))
colnames(Top100_GO_MF_map) = c("Term", "Stage_1", "Stage_2", "Stage_3", "Stage_4", "Blood")
Top100_GO_MF_map$Term = sort(unique(c(pathfindR_outputs_blood$`GO-MF`$Term_Description[1:100], 
                                      pathfindR_outputs_stage_1$`GO-MF`$Term_Description[1:100],
                                      pathfindR_outputs_stage_2$`GO-MF`$Term_Description[1:100], 
                                      pathfindR_outputs_stage_3$`GO-MF`$Term_Description[1:100],
                                      pathfindR_outputs_stage_4$`GO-MF`$Term_Description[1:100])))
Top100_GO_MF_map$Stage_1 = ifelse(Top100_GO_MF_map$Term %in% pathfindR_outputs_stage_1$`GO-MF`$Term_Description[1:100],
                                  "Present", "Absent")
Top100_GO_MF_map$Stage_2 = ifelse(Top100_GO_MF_map$Term %in% pathfindR_outputs_stage_2$`GO-MF`$Term_Description[1:100],
                                  "Present", "Absent")
Top100_GO_MF_map$Stage_3 = ifelse(Top100_GO_MF_map$Term %in% pathfindR_outputs_stage_3$`GO-MF`$Term_Description[1:100],
                                  "Present", "Absent")
Top100_GO_MF_map$Stage_4 = ifelse(Top100_GO_MF_map$Term %in% pathfindR_outputs_stage_4$`GO-MF`$Term_Description[1:100],
                                  "Present", "Absent")
Top100_GO_MF_map$Blood = ifelse(Top100_GO_MF_map$Term %in% pathfindR_outputs_blood$`GO-MF`$Term_Description[1:100],
                                "Present", "Absent")
Top100_GO_MF_map[which(Top100_GO_MF_map$Stage_1 == "Present" & Top100_GO_MF_map$Stage_2 == "Absent" &
                         Top100_GO_MF_map$Stage_3 == "Absent" & Top100_GO_MF_map$Stage_4 == "Absent" &
                         Top100_GO_MF_map$Blood == "Present"),]

# KEGG
Top100_KEGG_map = as.data.frame(matrix(nrow=length(unique(c(pathfindR_outputs_blood$KEGG$Term_Description[1:100], 
                                                            pathfindR_outputs_stage_1$KEGG$Term_Description[1:100],
                                                            pathfindR_outputs_stage_2$KEGG$Term_Description[1:100], 
                                                            pathfindR_outputs_stage_3$KEGG$Term_Description[1:100],
                                                            pathfindR_outputs_stage_4$KEGG$Term_Description[1:100]))), ncol = 6))
colnames(Top100_KEGG_map) = c("Term", "Stage_1", "Stage_2", "Stage_3", "Stage_4", "Blood")
Top100_KEGG_map$Term = sort(unique(c(pathfindR_outputs_blood$KEGG$Term_Description[1:100], 
                                     pathfindR_outputs_stage_1$KEGG$Term_Description[1:100],
                                     pathfindR_outputs_stage_2$KEGG$Term_Description[1:100], 
                                     pathfindR_outputs_stage_3$KEGG$Term_Description[1:100],
                                     pathfindR_outputs_stage_4$KEGG$Term_Description[1:100])))
Top100_KEGG_map$Stage_1 = ifelse(Top100_KEGG_map$Term %in% pathfindR_outputs_stage_1$KEGG$Term_Description[1:100],
                                 "Present", "Absent")
Top100_KEGG_map$Stage_2 = ifelse(Top100_KEGG_map$Term %in% pathfindR_outputs_stage_2$KEGG$Term_Description[1:100],
                                 "Present", "Absent")
Top100_KEGG_map$Stage_3 = ifelse(Top100_KEGG_map$Term %in% pathfindR_outputs_stage_3$KEGG$Term_Description[1:100],
                                 "Present", "Absent")
Top100_KEGG_map$Stage_4 = ifelse(Top100_KEGG_map$Term %in% pathfindR_outputs_stage_4$KEGG$Term_Description[1:100],
                                 "Present", "Absent")
Top100_KEGG_map$Blood = ifelse(Top100_KEGG_map$Term %in% pathfindR_outputs_blood$KEGG$Term_Description[1:100],
                               "Present", "Absent")
Top100_KEGG_map[which(Top100_KEGG_map$Stage_1 == "Present" & Top100_KEGG_map$Stage_2 == "Absent" &
                        Top100_KEGG_map$Stage_3 == "Absent" & Top100_KEGG_map$Stage_4 == "Absent" &
                        Top100_KEGG_map$Blood == "Present"),]

# Reactome
Top100_Reactome_map = as.data.frame(matrix(nrow=length(unique(c(pathfindR_outputs_blood$Reactome$Term_Description[1:100], 
                                                                pathfindR_outputs_stage_1$Reactome$Term_Description[1:100],
                                                                pathfindR_outputs_stage_2$Reactome$Term_Description[1:100], 
                                                                pathfindR_outputs_stage_3$Reactome$Term_Description[1:100],
                                                                pathfindR_outputs_stage_4$Reactome$Term_Description[1:100]))), ncol = 6))
colnames(Top100_Reactome_map) = c("Term", "Stage_1", "Stage_2", "Stage_3", "Stage_4", "Blood")
Top100_Reactome_map$Term = sort(unique(c(pathfindR_outputs_blood$Reactome$Term_Description[1:100], 
                                         pathfindR_outputs_stage_1$Reactome$Term_Description[1:100],
                                         pathfindR_outputs_stage_2$Reactome$Term_Description[1:100], 
                                         pathfindR_outputs_stage_3$Reactome$Term_Description[1:100],
                                         pathfindR_outputs_stage_4$Reactome$Term_Description[1:100])))
Top100_Reactome_map$Stage_1 = ifelse(Top100_Reactome_map$Term %in% pathfindR_outputs_stage_1$Reactome$Term_Description[1:100],
                                     "Present", "Absent")
Top100_Reactome_map$Stage_2 = ifelse(Top100_Reactome_map$Term %in% pathfindR_outputs_stage_2$Reactome$Term_Description[1:100],
                                     "Present", "Absent")
Top100_Reactome_map$Stage_3 = ifelse(Top100_Reactome_map$Term %in% pathfindR_outputs_stage_3$Reactome$Term_Description[1:100],
                                     "Present", "Absent")
Top100_Reactome_map$Stage_4 = ifelse(Top100_Reactome_map$Term %in% pathfindR_outputs_stage_4$Reactome$Term_Description[1:100],
                                     "Present", "Absent")
Top100_Reactome_map$Blood = ifelse(Top100_Reactome_map$Term %in% pathfindR_outputs_blood$Reactome$Term_Description[1:100],
                                   "Present", "Absent")
Top100_Reactome_map[which(Top100_Reactome_map$Stage_1 == "Present" & Top100_Reactome_map$Stage_2 == "Absent" &
                            Top100_Reactome_map$Stage_3 == "Absent" & Top100_Reactome_map$Stage_4 == "Absent" &
                            Top100_Reactome_map$Blood == "Present"),]

#####
# Load volcano plots from tumor stage DGEA as an .RData file:
# load("/Your/path/your_volcano_plots.RData")

# All-stage(+blood)-volcano-dotplot (BioCarta)-pairs
tiff("Additional_plots/all_Volcano_CE_Dotplot_multiplot.tif", 
     width = 4320, height = 7680, res = 100)
multiplot(union_one_normal_volcano, union_two_normal_volcano, 
              union_three_normal_volcano, union_four_normal_volcano, 
          TN_z_volcano,
          cluster_enrichment_dotplots_stage_1[["BioCarta"]],
          cluster_enrichment_dotplots_stage_2[["BioCarta"]],
          cluster_enrichment_dotplots_stage_3[["BioCarta"]],
          cluster_enrichment_dotplots_stage_4[["BioCarta"]],
          cluster_enrichment_dotplots_blood[["BioCarta"]], cols = 2)
m = ggplot(multiplot(union_one_normal_volcano, union_two_normal_volcano, 
                     union_three_normal_volcano, union_four_normal_volcano,
                     TN_z_volcano,
                     cluster_enrichment_dotplots_stage_1[["BioCarta"]],
                     cluster_enrichment_dotplots_stage_2[["BioCarta"]],
                     cluster_enrichment_dotplots_stage_3[["BioCarta"]],
                     cluster_enrichment_dotplots_stage_4[["BioCarta"]],
                     cluster_enrichment_dotplots_blood[["BioCarta"]], cols = 2))
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

# Blood
tiff("Additional_plots/Blood_Volcano_CE_Dotplot.tif", 
     width = 1920, height = 1080, res = 150)
multiplot(TN_z_volcano, 
          cluster_enrichment_dotplots_blood[["BioCarta"]], cols = 2)
m = ggplot(multiplot(TN_z_volcano, 
                     cluster_enrichment_dotplots_blood[["BioCarta"]], cols = 2))
dev.off(); rm(m)

# DEmiRNAs #####
DEmiRNAs_stage_1 = read.xlsx("DGEA/DEmiRNAs.xlsx", sheet = 1, cols = c(1:3))$Gene.Symbol
DEmiRNAs_stage_2 = read.xlsx("DGEA/DEmiRNAs.xlsx", sheet = 2, cols = c(1:3))$Gene.Symbol
DEmiRNAs_stage_3 = read.xlsx("DGEA/DEmiRNAs.xlsx", sheet = 3, cols = c(1:3))$Gene.Symbol
DEmiRNAs_stage_4 = read.xlsx("DGEA/DEmiRNAs.xlsx", sheet = 4, cols = c(1:3))$Gene.Symbol
DEmiRNAs_blood = read.xlsx("DGEA/DEmiRNAs.xlsx", sheet = 5, cols = c(1:3))$Gene.Symbol

# exclude host genes and anti-sense
clean_DEmiRNAs_stage_1 = read.xlsx("DGEA/DEmiRNAs.xlsx", sheet = 1, cols = c(8:10))$Gene.Symbol
clean_DEmiRNAs_stage_2 = read.xlsx("DGEA/DEmiRNAs.xlsx", sheet = 2, cols = c(8:10))$Gene.Symbol
clean_DEmiRNAs_stage_3 = read.xlsx("DGEA/DEmiRNAs.xlsx", sheet = 3, cols = c(8:10))$Gene.Symbol
clean_DEmiRNAs_stage_4 = read.xlsx("DGEA/DEmiRNAs.xlsx", sheet = 4, cols = c(8:10))$Gene.Symbol
clean_DEmiRNAs_blood = read.xlsx("DGEA/DEmiRNAs.xlsx", sheet = 5, cols = c(8:10))$Gene.Symbol

miRNA_venn = list(`Stage 1` = DEmiRNAs_stage_1,
                  `Stage 2` = DEmiRNAs_stage_2,
                  `Stage 3` = DEmiRNAs_stage_3,
                  `Stage 4` = DEmiRNAs_stage_4,
                  `Blood samples` = DEmiRNAs_blood)

clean_miRNA_venn = list(`Stage 1` = clean_DEmiRNAs_stage_1,
                        `Stage 2` = clean_DEmiRNAs_stage_2,
                        `Stage 3` = clean_DEmiRNAs_stage_3,
                        `Stage 4` = clean_DEmiRNAs_stage_4,
                        `Blood samples` = clean_DEmiRNAs_blood)

rm(DEmiRNAs_stage_1, DEmiRNAs_stage_2, DEmiRNAs_stage_3, DEmiRNAs_stage_4, DEmiRNAs_blood,
   clean_DEmiRNAs_stage_1, clean_DEmiRNAs_stage_2, clean_DEmiRNAs_stage_3, 
   clean_DEmiRNAs_stage_4, clean_DEmiRNAs_blood); gc()

DEmiRNAs_venn = ggVennDiagram(
  miRNA_venn, label_alpha = 0,
  category.names = names(miRNA_venn)
) +
  scale_fill_gradient(low = "white", high = "darkred") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18))+
  labs(title = "Venn diagram of sig. diff. expressed miRNAs (+ host genes)")
tiff("Additional_plots/Venn/ggVennDiagram_blood_and_Stages_miRNA_Venn.tif", 
     width = 1920, height = 1080, res = 150)
DEmiRNAs_venn
dev.off()

clean_DEmiRNAs_venn = ggVennDiagram(
  clean_miRNA_venn, label_alpha = 0,
  category.names = names(clean_miRNA_venn)
) +
  scale_fill_gradient(low = "white", high = "darkred") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18))+
  labs(title = "Venn diagram of sig. diff.expressed miRNAs")
tiff("Additional_plots/Venn/clean_ggVennDiagram_Stages_miRNA_Venn.tif", 
     width = 1920, height = 1080, res = 150)
clean_DEmiRNAs_venn
dev.off()

tiff("Additional_plots/Venn/ggVennDiagram_DEmiRNA_multiplot.tif", 
     width = 1920, height = 1080, res = 120)
multiplot(DEmiRNAs_venn, clean_DEmiRNAs_venn, cols = 2)
m = ggplot(multiplot(DEmiRNAs_venn, clean_DEmiRNAs_venn, cols = 2))
dev.off(); rm(m)

# Exploring miR21 further
mir21_1 = read.xlsx("Mienturnet/MIR21_predictions_(unique_common_stage_1_blood_up).xlsx",
                    sheet = 1) %>%
  dplyr::select(-X7)
mir21_2 = read.xlsx("Mienturnet/MIR21_predictions_(unique_common_stage_1_blood_up).xlsx",
                    sheet = 2) %>%
  dplyr::select(-X7)
mir21 = rbind(mir21_1, mir21_2); rm(mir21_1, mir21_2)
mir21_DEG_targets_stage_1 = intersect(mir21$Gene.Symbol, sig_DGEA[["Stage_1"]][["Gene.Symbol_pre"]])
mir21_DEG_targets_blood = intersect(mir21$Gene.Symbol, sig_DGEA[["Blood"]][["Gene.Symbol_pre"]])
sig_DGEA[["Stage_1"]][sig_DGEA[["Stage_1"]][["Gene.Symbol_pre"]]
                      %in% mir21_DEG_targets_stage_1,]
sig_DGEA[["Blood"]][sig_DGEA[["Blood"]][["Gene.Symbol_pre"]]
                    %in% mir21_DEG_targets_blood,]