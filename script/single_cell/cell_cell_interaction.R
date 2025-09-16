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
sce <- readRDS('zenodo/sc_IBD_GVHD_merge_raw_250825.RDS')
sce <- subset(sce,
              subset = celltype%in%c("MDSC-MS4A3","MDSC-CAMP","Neu-IFI6","Neu-TLR2",'Neu-CD177','ins_B','ins_DC_Macro','ins_Epi','ins_Endo','ins_FLC','ins_T_NK'))
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
CellChatDB$interaction <- read.csv(file = 'zenodo/human_interaction_input_CellChatDB.csv', row.names = 1)
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
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchat@net$weight
levels(cellchat@idents)
vertex.receiver = c(1,2)
pathways.show.all <- cellchat@netP$pathways
netAnalysis_signalingRole_scatter(cellchat)


