library(DESeq2)
library(apeglm)
#library(IHW)
library(org.Hs.eg.db)
library(biomaRt)
library(curl)
library(ggplot2)
library(Seurat)
library(dplyr)
dds <- readRDS("zenodo/human_Colon_raw_count_matrix.RDS")
#----------------------- Generating dds object -----------------------------------
meta_data <- read.csv('E://GVHD/example_bulk_RNA/human_Colon_raw_count_matrix_metainfo.csv',check.names = F)
rownames(meta_data) <- meta_data$sampleName
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = meta_data, design = ~ group)
#keep <- rowSums(counts(dds)) >= 1
keep <- apply(counts(dds), 1, max) > 10
dds <- dds[keep,]
dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)
res <- results(dds, contrast = c("group","GVHD","CT"), altHypothesis = 'greater',lfcThreshold = 0)
res <- na.omit(res)
res <- res[res$log2FoldChange>0,]
res <- res[res$padj<0.05,]
# calculating gene set score
raw_matrix <- assay(dds)
raw_matrix <- as.data.frame(raw_matrix)
sce <- CreateSeuratObject(counts = raw_matrix)
sce$group <- meta_data$group
Idents(sce) <- sce$group
sce <- PercentageFeatureSet(sce, pattern = "^MT-", col.name = "percent.mt")
sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sce)
sce <- ScaleData(sce, features = all.genes)
temp <- read.csv("Top30_Neu_CD177_up_DEG.csv", header = F)
cd_features <- list(temp$V1)
sce <- AddModuleScore(sce, features = cd_features, name = 'Neutro_score')
VlnPlot(sce, features = 'Neutro_score1')
temp <- sce@meta.data[,c('group','Neutro_score1')]
temp %>% mutate(group = factor(group, levels = c('CT','GVHD'))) %>%
  ggplot(aes(group,Neutro_score1)) +
  geom_boxplot(width = 0.1) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = 'black',fill = NA)) +
  scale_fill_manual(values=c('#003366','#0B6A9D','#60C4DD','#FCE2DF','#E94C47'))
