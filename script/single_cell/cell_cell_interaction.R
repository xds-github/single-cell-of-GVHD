# convert the h5
library(Seurat)
library(SeuratDisk)
library(CellChat)
library(patchwork)
library(rhdf5)
factor_order_change <- function(new_order, old_factor){
  new_order <- factor(1:length(new_order),labels = new_order)
  new_factor <- factor(old_factor,levels = levels(new_order))
  return(new_factor)
}
set_celltype <- function(sce, new.cluster.ids = c(1,2)){
  names(new.cluster.ids) <- levels(sce)
  sce <- RenameIdents(sce, new.cluster.ids)
  sce$celltype <- Idents(sce)
  return(sce)
}
# Processing blood data
mydata <- h5read("/mnt/e/GVHD/analysis/matrix_blood_raw_merge_231214.h5","mat")
mat <- mydata$block0_values
rownames(mat) <- mydata$axis0
colnames(mat) <- mydata$axis1
mat <- CreateAssayObject(mat)
meta <- read.csv("/mnt/e/GVHD/analysis/meta_subset1k_231214.csv", row.names = 1)
sce <- CreateSeuratObject(mat,assay='RNA',meta.data=meta)
sce <- subset(sce, subset = group%in%c('GE','CT'))
sce_blood <- sce
# Processing IBD data
sce <- ReadMtx("/mnt/e/GVHD/IBD_data/Colitis_V3_combined_raw_matrix_IBD.mtx",
               "/mnt/e/GVHD/IBD_data/Colitis_V3_combined_IBD_barcodes.txt",
               "/mnt/e/GVHD/IBD_data/Colitis_V3_combined_IBD_genes.txt",
               cell.column = 1,
               feature.column = 1)
sce <- CreateAssayObject(sce)
sce <- CreateSeuratObject(sce, assay='RNA')
meta_B <- read.csv("/mnt/e/GVHD/IBD_data/B_meta.csv")
meta_T <- read.csv("/mnt/e/GVHD/IBD_data/T_NK_meta.csv")
meta_FLC <- read.csv("/mnt/e/GVHD/IBD_data/FLC_meta.csv")
meta_Endo <- read.csv("/mnt/e/GVHD/IBD_data/Endo_meta.csv")
meta_Epi <- read.csv("/mnt/e/GVHD/IBD_data/Epi_meta.csv")
meta_DC <- read.csv("/mnt/e/GVHD/IBD_data/DC_Macro_meta.csv")
meta_info <- data.frame()
meta_info <- rbind(meta_B,meta_T,meta_FLC,meta_Endo,meta_Epi,meta_DC)
sce$barcode <- rownames(sce@meta.data)
rownames(meta_info) <- meta_info$barcode
sce <- subset(sce, subset = barcode%in%meta_info$barcode)
sce$celltype <- meta_info[rownames(sce@meta.data),'celltype']
Idents(sce) <- sce$celltype
#saveRDS(sce,"E://GVHD/IBD_data/IBD_subset_raw.RDS")
# merge total
sce_merge <- merge(sce_blood,sce,collapse=F)
sce_merge$ID <- sce_merge$group <- sce_merge$group2 <- sce_merge$barcode <- NULL
saveRDS(sce_merge, "Zenodo/IBD_GVHD_merge_raw_260521.RDS")

#
sce <- sce_merge
sce <- subset(sce,
              subset = celltype%in%c("imNeu-MS4A3","imNeu-CAMP","Neu-IFI6","Neu-TLR2",'Neu-CD177','ins_B','ins_DC_Macro','ins_Epi','ins_Endo','ins_FLC','ins_T_NK'))
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sce)
sce <- ScaleData(sce, features = all.genes)
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
ElbowPlot(sce)
sce <- FindNeighbors(sce, dims = 1:15)
sce <- FindClusters(sce, resolution = 0.5)
sce <- RunUMAP(sce, dims = 1:15)
DimPlot(sce, group.by = 'celltype')
#sce$celltype_ori <- NULL
Idents(sce) <- sce$celltype
# Constructing cellchat object
cellchat <- createCellChat(object = sce, group.by = "celltype", assay = "RNA")
CellChatDB <- CellChatDB.human
options(stringsAsFactors = FALSE)
CellChatDB$interaction <- read.csv(file = 'Zenodo/human_interaction_input_CellChatDB.csv', row.names = 1)
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#future::plan("multicore", workers = 6) # do parallel
# Preprocessing the expression data for cell-cell communication analysis
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat@idents <- factor(cellchat@idents)
#cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
#saveRDS(cellchat, '/mnt/e/GVHD/analysis/Dealing_inter_GVHD_IBD_250825.RDS')
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchat@net$weight
levels(cellchat@idents)
vertex.receiver = c(1,2)
pathways.show.all <- cellchat@netP$pathways
netAnalysis_signalingRole_scatter(cellchat)
