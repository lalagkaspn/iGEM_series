# TCGA miRNA validation approaches
library(TCGAbiolinks)
library(dplyr)
library(limma)
library(ggplot2)
library(survminer)
library(survival)

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
                            dataAssy.exp@colData@listData$days_to_death,
                            dataAssy.exp@colData@listData$days_to_last_follow_up))
colnames(outcomes) = c("Sample.ID", "Vital Status", "Days to Death", "Days to Last Follow-up")
rownames(outcomes) = outcomes$Sample.ID
outcomes$Deceased = ifelse(outcomes$`Vital Status` == "Alive", FALSE, TRUE)

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
DescTools::CohenD(metagene$DRS[metagene$group=="Dead"], metagene$DRS[metagene$group=="Alive"],
                  correct = TRUE)
# t = -2.2014, p = 0.02902, d = -0.33: small effect size

t.test(metagene$URS[metagene$group=="Dead"], metagene$URS[metagene$group=="Alive"])
DescTools::CohenD(metagene$URS[metagene$group=="Dead"], metagene$URS[metagene$group=="Alive"],
                  correct = TRUE)
# t = 2.8918, p = 0.004363, d = 0.43: small effect size 

# Metaplots
# Without overlaid points
metaplot1 = ggplot(metagene, aes(x = group, y = DRS, fill = group)) + 
  geom_boxplot(width=0.35, outlier.size = 0.5, linewidth = 0.3)+
  scale_fill_brewer(palette = "RdPu") +
  scale_y_continuous(limits = c(-3.5, 3.5), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  ggsignif::geom_signif(comparisons = list(c("Dead", "Alive")), 
                        map_signif_level=TRUE, size = 0.25, textsize = 3) +
  annotate("text", x = 1.5, y = -1.8, size = 1.5, parse = TRUE,
           label = "italic(t) == -2.2\n") +
  annotate("text", x = 1.5, y = -2.0, size = 1.5,
           label = paste0("\n", 
                          "p = 0.03")) +
  annotate("text", x = 1.5, y = -2.4, size = 1.5,
           label = paste0("\n", 
                          "d = -0.33")) +
  geom_rect(aes(xmin = 1.3, xmax = 1.7, ymin = -1.65, ymax = -2.75),
            fill = "transparent", color = "black", linewidth = 0.25) +
  # geom_jitter(color="black", size=0.4, alpha=0.9, width = 0.2) +
  theme(panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 6),
        axis.title = element_text(face = "bold", size = 5),
        axis.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.text = element_text(size = 5),
        legend.key.size = unit(4, "mm")) +
  labs(x = "Sample type",
       y = "Mean DRS expression",
       title = "Boxplots of mean expression: DRS",
       fill = "Legend")
metaplot1

metaplot2 = ggplot(metagene, aes(x = group, y = URS, fill = group)) + 
  geom_boxplot(width=0.35, outlier.size = 0.5, linewidth = 0.3)+
  scale_fill_brewer(palette = "RdPu") +
  scale_y_continuous(limits = c(-3.5, 3.5), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  ggsignif::geom_signif(comparisons = list(c("Dead", "Alive")), 
                        map_signif_level=TRUE, size = 0.25, textsize = 3) +
  annotate("text", x = 1.5, y = -1.8, size = 1.5, parse = TRUE,
           label = "italic(t) == 2.9\n") +
  annotate("text", x = 1.5, y = -2.0, size = 1.5,
           label = paste0("\n", 
                          "p = 0.004")) +
  annotate("text", x = 1.5, y = -2.4, size = 1.5,
           label = paste0("\n", 
                          "d = 0.43")) +
  geom_rect(aes(xmin = 1.3, xmax = 1.7, ymin = -1.65, ymax = -2.75),
            fill = "transparent", color = "black", linewidth = 0.25) +
  # geom_jitter(color="black", size=0.4, alpha=0.9, width = 0.2) +
  theme(panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 6),
        axis.title = element_text(face = "bold", size = 5),
        axis.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.text = element_text(size = 5),
        legend.key.size = unit(4, "mm")) +
  labs(x = "Sample type",
       y = "Mean URS expression",
       title = "Boxplots of mean expression: URS",
       fill = "Legend")
metaplot2

# Metaplots with overlaid points (jitter)
metaplot1_jitter = ggplot(metagene, aes(x = group, y = DRS, fill = group)) + 
  geom_boxplot(width=0.35, outlier.size = 0.5, linewidth = 0.3)+
  scale_fill_brewer(palette = "RdPu") +
  scale_y_continuous(limits = c(-3.5, 3.5), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  geom_jitter(color="black", size=0.4, alpha=0.9, width = 0.2) +
  theme(panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 6),
        axis.title = element_text(face = "bold", size = 5),
        axis.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.text = element_text(size = 5),
        legend.key.size = unit(4, "mm")) +
  labs(x = "Sample type",
       y = "Mean DRS expression",
       title = "Boxplots of mean expression: DRS",
       fill = "Legend")
metaplot1_jitter

metaplot2_jitter = ggplot(metagene, aes(x = group, y = URS, fill = group)) + 
  geom_boxplot(width=0.35, outlier.size = 0.5, linewidth = 0.3)+
  scale_fill_brewer(palette = "RdPu") +
  scale_y_continuous(limits = c(-3.5, 3.5), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  geom_jitter(color="black", size=0.4, alpha=0.9, width = 0.2) +
  theme(panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 6),
        axis.title = element_text(face = "bold", size = 5),
        axis.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.text = element_text(size = 5),
        legend.key.size = unit(4, "mm")) +
  labs(x = "Sample type",
       y = "Mean URS expression",
       title = "Boxplots of mean expression: URS",
       fill = "Legend")
metaplot2_jitter

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

# Multiplot without overlaid points
tiff("Signatures/TCGA_signature_metaplots.tiff", 
     width = 1920*2, height = 1920, res = 700, compression = "lzw")
multiplot(metaplot1, metaplot2, cols = 2)
m = ggplot(multiplot(metaplot1, metaplot2, cols = 2))
dev.off(); rm(m)

# Multiplots with overlaid points
tiff("Signatures/TCGA_signature_metaplots_jitter.tiff", 
     width = 1920*2, height = 1920, res = 700, compression = "lzw")
multiplot(metaplot1_jitter, metaplot2_jitter, cols = 2)
m = ggplot(multiplot(metaplot1_jitter, metaplot2_jitter, cols = 2))
dev.off(); rm(m)

# Survival analysis #####
metagene$DRS_group = ifelse(metagene$DRS < 0, "Low", "High")
metagene$URS_group = ifelse(metagene$URS < 0, "Low", "High")
metagene$URS_group = factor(metagene$URS_group, levels = c("High", "Low"),
                            labels = c("High", "Low"))
metagene$DRS_group = factor(metagene$DRS_group, levels = c("High", "Low"),
                            labels = c("High", "Low"))

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
outcomes$`Days to Death` = as.numeric(outcomes$`Days to Death`)
outcomes$`Days to Last Follow-up` = as.numeric(outcomes$`Days to Last Follow-up`)
outcomes$overall_survival = ifelse(outcomes$`Vital Status` == "Alive",
                                   outcomes$`Days to Last Follow-up`,
                                   outcomes$`Days to Death`)

metagene_surv = metagene %>% inner_join(outcomes, by = "Sample.ID")

# fitting survival curve -----------
# URS
URS_survfit = survfit(Surv(overall_survival, Deceased) ~ URS_group, data = metagene_surv)
URS_survfit

survurs = ggsurvplot(URS_survfit,
                     title = "Survival curves for URS groups",
                     data = metagene_surv,
                     pval = T,
                     size = 0.35,
                     censor.size = 2.5,
                     pval.size = 2,
                     risk.table = F,
                     legend.labs = c("URS > 0", "URS < 0"),
                     legend = c(0.9, 1),
                     font.legend = 4,
                     ggtheme = theme_classic()+
                       theme(plot.title = element_text(face = "bold", size = 6),
                             axis.title = element_text(face = "bold", size = 5),
                             axis.text = element_text(size = 5)))

tiff("Signatures/URS_survival_curve.tiff", width = 1920, height = 1920,
     res = 700, units = "px", compression = "lzw")
survurs
dev.off()

URS_survfit2 = survdiff(Surv(overall_survival, Deceased) ~ URS_group, data = metagene_surv)

# DRS
DRS_survfit = survfit(Surv(overall_survival, Deceased) ~ DRS_group, data = metagene_surv)
DRS_survfit

survdrs = ggsurvplot(DRS_survfit,
                     title = "Survival curves for DRS groups",
                     data = metagene_surv,
                     pval = T,
                     size = 0.35,
                     censor.size = 2.5,
                     pval.size = 2,
                     risk.table = F,
                     legend.labs = c("DRS > 0", "DRS < 0"),
                     legend = c(0.9, 1),
                     font.legend = 4,
                     ggtheme = theme_classic()+
                       theme(plot.title = element_text(face = "bold", size = 6),
                             axis.title = element_text(face = "bold", size = 5),
                             axis.text = element_text(size = 5)))

tiff("Signatures/DRS_survival_curve.tiff", width = 1920, height = 1920,
     res = 700, units = "px", compression = "lzw")
survdrs
dev.off()

DRS_survfit2 = survdiff(Surv(overall_survival, Deceased) ~ DRS_group, data = metagene_surv)

# All plots
ggpubr::ggarrange(metaplot2, metaplot1, survurs$plot, survdrs$plot,
                  ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"),
                  font.label = list(size = 5))
ggsave("TCGA_metaplots_and_survcurves.tiff",
       path = "Signatures/",
       width = 1920*2, height = 1920*2, dpi = 700, compression = "lzw",
       units = "px", device = "tiff")
