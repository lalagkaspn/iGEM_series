# Identification of early diagnosis markers of pancreatic ductal adenocarcinoma (PDAC) using publicly available transcriptomic tumour and blood sample data

## Aristeidis Sionakidis, Panagiotis Nikolaos Lalagkas, Andigoni Malousi, Ioannis S. Vizirianakis.

## Code authors: Aristeidis Sionakidis, Panagiotis Nikolaos Lalagkas

Link to the [paper](https://doi.org/10.1002/ctd2.248).
DOI: [10.1002/ctd2.248](https://doi.org/10.1002/ctd2.248)

In this work we identify a set of *820 genes* which are consistently deregulated across tumor stages and sample types in PDAC and whose predictive potential is validated on The Cancer Genome Atlas (*TCGA*) data. We use the `pathfindR` package to identify enriched networks and map the genes of these networks to approved drugs (DrugBank) to explore pharmacogenomic relationships. Additionally, we focus on differentially expressed microRNAs (*miRNAs*) and enriched miRNAs as these are predicted by [Mienturnet](http://userver.bio.uniroma1.it/apps/mienturnet/) based on the differentially expressed genes we found. We finally generate Circos plots which illustrate the pharmacogenomic relationships of interest along with miRNA-gene interactions and which also provide additional annotations for genes (e.g. mutational status, cancer driver status).

There are five scripts:

- `GSE_tumor_stage.R` is the script used to download, preprocess and conduct Differential Gene Expression Analysis (*DGEA*) on clinical PDAC samples with tumor staging information available.
- `GSE_blood_samples.R` is the script used to download, preprocess and conduct DGEA on clinical PDAC blood samples (results not adjusted for tumor stage).
- `Pathway_ontology_driver_genes.R` contains the code for performing Active Subnetwork Enrichment analysis with the `pathfindR` package, annotating genes of interest with cancer driver information from the Catalogue Of Somatic Mutations In Cancer (*COSMIC*) and checking the overlaps between significant genes and pathways across stages and sample types.
- `Circos.R` is the script in which the code for generating all the files necessary to produce the Circos plots are generated.
- `TCGA_validations.R` contains the code for downloading relative TCGA PDAC data and validating our signature on them as well as conducting survival analysis.
- `Gut_panc_dgea.R` contains the code for the external validation of our signature in a large cohort of normal, PDAC and chronic pancreatitis samples.

This repository contains all the work conducted for our PDAC work. All scripts are written in **R**. The vast majority of the output is **HTML** plots generated from the `pathfindR` package.
