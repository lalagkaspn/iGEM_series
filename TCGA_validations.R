# TCGA miRNA validation approaches
library(TCGAbiolinks)
library(dplyr)
library(limma)

query.miRNA = GDCquery(
  project = "TCGA-PAAD", 
  experimental.strategy = "miRNA-Seq",
  data.category = "Transcriptome Profiling", 
  data.type = "miRNA Expression Quantification"
)

GDCdownload(query = query.miRNA, files.per.chunk = 183)

dataAssy.miR = GDCprepare(
  query = query.miRNA
)
rownames(dataAssy.miR) = dataAssy.miR$miRNA_ID

# using read_count's data 
read_countData =  colnames(dataAssy.miR)[grep("count", colnames(dataAssy.miR))]
dataAssy.miR = dataAssy.miR[,read_countData]
colnames(dataAssy.miR) = gsub("read_count_","", colnames(dataAssy.miR))

dataFilt = TCGAanalyze_Filtering(
  tabDF = dataAssy.miR,
  method = "quantile", 
  qnt.cut =  0.25
)

# z-score-transformation 
# KBZ transformation method ( https:://www.biostars.org/p/283083/ )
t = as.data.frame(t(dataFilt))
z_TCGA_miRNA_norm = sapply(t, function(t) (t-mean(t, na.rm = T))/sd(t, na.rm = T))
z_TCGA_miRNA_norm = as.matrix(t(z_TCGA_miRNA_norm))
rownames(z_TCGA_miRNA_norm) = rownames(dataFilt)
colnames(z_TCGA_miRNA_norm) = colnames(dataFilt)
rm(t)

# Prepare a design matrix
samplesNT = TCGAquery_SampleTypes(colnames(z_TCGA_miRNA_norm), typesample = c("NT"))
samplesTP = TCGAquery_SampleTypes(colnames(z_TCGA_miRNA_norm), typesample = c("TP"))
pheno = rbind(as.data.frame(list(Sample.ID = samplesNT, Type = "Normal")),
              as.data.frame(list(Sample.ID = samplesTP, Type = "Tumor")))
pheno$Type = factor(pheno$Type, levels = c("Normal", "Tumor"),
                    labels = c("Normal", "Tumor"))

design = model.matrix(~0 + pheno$Type)                            
colnames(design) = c("Normal", "Tumor")                            
cont.matrix = makeContrasts(TumorvsNormal = Tumor - Normal, levels = design)

fit1 = lmFit(z_TCGA_miRNA_norm[,pheno$Sample.ID], design)
fit2 = contrasts.fit(fit1, cont.matrix)
fit3 = eBayes(fit2, robust = TRUE)
summary(decideTests(fit3))
res = topTable(fit3, adjust.method="BH", 
         number = Inf, coef = "TumorvsNormal")

# No miRNA results are replicated here

# This replicates miR-192
dataDEGs = TCGAanalyze_DEA(
  mat1 = dataFilt[,samplesTP],
  mat2 = dataFilt[,samplesNT],
  Cond1type = "Tumor",
  Cond2type = "Normal",
  fdr.cut = 0.1,
  logFC.cut = 0,
  method = "glmLRT"
) %>%
  arrange(., FDR) %>%
  mutate(miRNA = rownames(.)) %>%
  mutate(microRNA = gsub("hsa-mir", "miR", miRNA)) %>%
  dplyr::select(-miRNA)

################################################################################
#                                                                              #
#                  TCGA expression validation approaches                       #
#                                                                              #
################################################################################
query.exp = GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification"
)

GDCdownload(query = query.exp, files.per.chunk = 183)

dataAssy.exp = GDCprepare(
  query = query.exp
)

rownames(dataAssy.exp) = dataAssy.exp@rowRanges@elementMetadata@listData[["gene_name"]]

# TPM data that were log2(TPM +1) transformed
prenorm_exp = log2(dataAssy.exp@assays@data@listData[["tpm_unstrand"]] + 1)
rownames(prenorm_exp) = rownames(dataAssy.exp)
colnames(prenorm_exp) = colnames(dataAssy.exp)

# z-score-transformation 
# KBZ transformation method ( https:://www.biostars.org/p/283083/ )
t = as.data.frame(t(prenorm_exp))
z_TCGA_exp_norm = sapply(t, function(t) (t-mean(t, na.rm = T))/sd(t, na.rm = T))
z_TCGA_exp_norm = as.matrix(t(z_TCGA_exp_norm))
rownames(z_TCGA_exp_norm) = rownames(prenorm_exp)
colnames(z_TCGA_exp_norm) = colnames(prenorm_exp)
rm(t)

# Validate the metascores in TCGA data #####

# Import the 820-gene signature
signature = read.xlsx("DGEA/all_stage_blood_concordant_overlap_set.xlsx")

# Filter TCGA data for outcomes of interest
# Variables: names(dataAssy.exp@colData@listData)
outcomes = data.frame(cbind(dataAssy.exp@colData@listData$barcode,
                            dataAssy.exp@colData@listData$vital_status,
                            dataAssy.exp@colData@listData$days_to_death))
colnames(outcomes) = c("Sample.ID", "Vital Status", "Days to Death")
rownames(outcomes) = outcomes$Sample.ID

metaexp = z_TCGA_exp_norm[intersect(rownames(z_TCGA_exp_norm), signature$Gene.Symbol),]
metagene = data.frame(colMeans(metaexp, na.rm = TRUE)) %>%
  mutate(rownames(.))
colnames(metagene) = c("Metascore", "Sample.ID")
metagene = metagene[metagene$Sample.ID %in% outcomes$Sample.ID,]
metagene$group = factor(outcomes$`Vital Status`, 
                        levels = c("Dead", "Alive"),
                        labels = c("Dead", "Alive"))
metagene$DRS = colMeans(metaexp[intersect(rownames(z_TCGA_exp_norm), signature$Gene.Symbol[signature$logFC_stage_1<0]),], na.rm = TRUE)
metagene$URS = colMeans(metaexp[intersect(rownames(z_TCGA_exp_norm), signature$Gene.Symbol[signature$logFC_stage_1>0]),], na.rm = TRUE)

t.test(metagene$DRS[metagene$group=="Dead"], metagene$DRS[metagene$group=="Alive"])
# t = -2.2014, p = 0.02902

t.test(metagene$URS[metagene$group=="Dead"], metagene$URS[metagene$group=="Alive"])
# t = 2.8918, p = 0.004363 

# Metaplots
metaplot1 = ggplot(metagene, aes(x = group, y = DRS, fill = group)) + 
  geom_boxplot(width=0.35)+
  scale_fill_brewer(palette = "RdPu") +
  theme(panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.line = element_line(),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
  labs(x = "Sample type",
       y = "Mean expression of down-regulated genes",
       title = "Boxplots of mean expression: down-regulated genes",
       fill = "Legend")
metaplot1

metaplot2 = ggplot(metagene, aes(x = group, y = URS, fill = group)) + 
  geom_boxplot(width=0.35)+
  scale_fill_brewer(palette = "RdPu") +
  theme(panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.line = element_line(),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
  labs(x = "Sample type",
       y = "Mean expression of up-regulated genes",
       title = "Boxplots of mean expression: up-regulated genes",
       fill = "Legend")
metaplot2

# Defining the multiplot function

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

# End of multiplot function
tiff("Signatures/TCGA_signature_metaplots.tif", 
     width = 1920, height = 1080, res = 150)
multiplot(metaplot1, metaplot2, cols = 2)
m = ggplot(multiplot(metaplot1, metaplot2, cols = 2))
dev.off(); rm(m)