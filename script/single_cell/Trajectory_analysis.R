library(rhdf5)
library(Matrix)
library(Seurat)
library(monocle3)
library(ComplexHeatmap)
library(circlize)
library(SeuratDisk)
library(RColorBrewer)
library(ggplot2)
# import scRNA-seq data
mydata <- h5read("zenodo/sc_matrix_Granucyte_4w_240422_raw.h5","mat")
mat <- mydata$block0_values
rownames(mat) <- mydata$axis0
colnames(mat) <- mydata$axis1
mat <- CreateAssayObject(mat)
meta <- read.csv("zenodo/sc_meta_Granucyte_4w_240422.csv", row.names = 1)
rownames(meta) <- meta$barcode
meta <- meta[colnames(mat),]
sce <- CreateSeuratObject(mat,assay='RNA',meta.data=meta)
#sce <- subset(sce, subset = group == 'GE')
data <- GetAssayData(sce, assay = 'RNA', slot = 'counts')
cell_metadata <- sce@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 100)
cds <- align_cds(cds, alignment_group = "ID")
cds <- reduce_dimension(cds, reduction_method = 'UMAP',umap.min_dist = 0.4)
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- read.csv('Zenodo/sc_Granucyte_4w_umap.csv', row.names = 1)
int.embed <- int.embed[rownames(cds.embed),]
cds.embed[,1] <- int.embed[,'X_X']
cds.embed[,2] <- int.embed[,'X_Y']
cds@int_colData$reducedDims$UMAP <- cds.embed
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "celltype")
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "group")
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "ID")

cds <- cluster_cells(cds)
cds <- learn_graph(cds, learn_graph_control = list(nn.k=25))
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           cell_size = 1)

cds_subset <- choose_cells(cds)
cds_subset <- order_cells(cds_subset)
pt.matrix <- exprs(cds_subset)[match(genes,rownames(rowData(cds_subset))),order(pseudotime(cds_subset))]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=5)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
pt.matrix <- na.omit(pt.matrix)
htkm <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "RdBu"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  row_title_rot                = 0,
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)
htkm = draw(htkm)
cell_pseudotimes <- pseudotime(cds)
df <- data.frame(barcode = colnames(cds), pseudotime = cell_pseudotimes)
