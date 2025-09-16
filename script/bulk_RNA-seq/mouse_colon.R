library(DESeq2)
library(apeglm)
#library(IHW)
library(org.Hs.eg.db)
library(biomaRt)
library(curl)
library(ggplot2)
library(Seurat)
# ----------------------- Generating dds object -----------------------------------
meta_data <- read.csv('zenodo/bulk_RNA_mouse_240509.csv',check.names = F)
meta_data$group2 <- paste(meta_data$group, meta_data$treatment, sep = '_')
meta_data <- meta_data[meta_data$is_included=='Y',]
rownames(meta_data) <- meta_data$sampleName
count_matrix <- readRDS('zenodo/bulk_RNA_mouse_240509_Colon_raw_count_matrix.RDS')
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = meta_data, design = ~ group2)
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast = c("group2","GVHD_Vehicle","CTL_Vehicle"), altHypothesis = 'greater',lfcThreshold = 0)
res <- na.omit(res)
res <- res[res$log2FoldChange>0,]
res <- res[res$padj<0.05,]
# calculating gene set score
raw_matrix <- assay(dds)
raw_matrix <- as.data.frame(raw_matrix)
sce <- CreateSeuratObject(counts = raw_matrix)
sce$group <- meta_data$group2
Idents(sce) <- sce$group
sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sce)
sce <- ScaleData(sce, features = all.genes)
sce <- subset(sce, subset = group%in%c('CTL_Vehicle','GVHD_Vehicle'))
temp <- read.csv("zenodo/Top30_Neu_CD177_up_DEG.csv", header = F)
temp <- temp$V1
efer_gene_list <- read.table('zenodo/Ensembl_human_mouse_ID_change.txt', sep = '\t', header = T,check.names = F)
refer_gene_list <- refer_gene_list[refer_gene_list$`Gene name`%in%temp,]
refer_gene_list <- refer_gene_list[!duplicated(refer_gene_list$`Gene name`),]
refer_gene_list <- refer_gene_list[refer_gene_list$`Mouse gene name`!='',]
cd_features <- list(refer_gene_list$`Mouse gene name`)
sce <- AddModuleScore(sce, features = cd_features, name = 'Neutro_score')
VlnPlot(sce, features = 'Neutro_score1')
temp <- sce@meta.data[,c('group','Neutro_score1')]
wilcox.test(temp[temp$group=='CTL_Vehicle','Neutro_score1'],temp[temp$group=='GVHD_Vehicle','Neutro_score1'])
temp %>% mutate(group = factor(group, levels = c('CTL_Vehicle','GVHD_Vehicle'))) %>%
  ggplot(aes(group, Neutro_score1, fill = group)) + geom_boxplot(width = 0.5) +
  geom_jitter(aes(fill  = group),width = 0.3, size = 5, shape = 21)+
  scale_fill_manual(values=c("#91D1C2", "#FEB500")) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = 'black',fill = NA)) +
  labs(title = 'CD177+ Neu Score')

